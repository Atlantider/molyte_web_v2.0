from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, List, Tuple
from pathlib import Path
import logging

class BaseQCEngine(ABC):
    """Abstract base class for Quantum Chemistry engines"""
    
    def __init__(self, config: Dict[str, Any], logger: Any):
        self.config = config
        self.logger = logger
        
    @abstractmethod
    def generate_input(self, job: Dict, work_dir: Path) -> Path:
        pass

    @abstractmethod
    def generate_job_script(self, job: Dict, input_path: Path, work_dir: Path) -> Path:
        pass
    
    @abstractmethod
    def submit_job(self, job_script: Path, work_dir: Path) -> Dict[str, Any]:
        pass
        
    @abstractmethod
    def parse_result(self, work_dir: Path) -> Dict[str, Any]:
        pass

    def _sanitize_filename(self, name: str) -> str:
        """Helper to sanitize filenames"""
        import re
        safe = re.sub(r'[^\w\-.]', '_', name)
        safe = re.sub(r'_+', '_', safe)
        safe = safe.strip('_')
        return safe or 'molecule'

    def _parse_xyz_content(self, xyz_content: str) -> Optional[List[Tuple[str, float, float, float]]]:
        """Parse XYZ format content"""
        if not xyz_content or not xyz_content.strip():
            return None

        try:
            lines = xyz_content.strip().split('\n')
            if len(lines) < 3:
                return None

            try:
                num_atoms = int(lines[0].strip())
            except ValueError:
                return None

            coords = []
            for line in lines[2:]:
                line = line.strip()
                if not line: continue
                parts = line.split()
                if len(parts) < 4: continue
                try:
                    coords.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
                except ValueError:
                    continue
            
            return coords if coords else None
        except Exception as e:
            self.logger.error(f"Failed to parse XYZ content: {e}")
            return None

    def _get_3d_coordinates(self, smiles: str, molecule_name: str) -> Optional[List[Tuple[str, float, float, float]]]:
        """Get 3D coordinates from SMILES or PDB"""
        # Try loading from PDB first
        initial_salts_path = Path(self.config['local'].get('initial_salts_path', '/public/home/xiaoji/initial_salts'))
        
        base_name = molecule_name.split('-')[0] if '-' in molecule_name else molecule_name
        base_name = base_name.split('_')[0] if '_' in molecule_name else base_name
        clean_name = molecule_name.replace("+", "").replace("-", "").strip()
        clean_base_name = base_name.replace("+", "").replace("-", "").strip()

        possible_paths = [
            initial_salts_path / f"{base_name}.pdb",
            initial_salts_path / f"{clean_base_name}.pdb",
            initial_salts_path / f"{molecule_name}.pdb",
            initial_salts_path / f"{clean_name}.pdb",
        ]

        for pdb_path in possible_paths:
            if pdb_path.exists():
                coords = self._parse_pdb_coordinates(pdb_path)
                if coords:
                    self.logger.info(f"Loaded coordinates from PDB: {pdb_path}")
                    return coords

        if not smiles:
            return None

        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mol = Chem.AddHs(mol)
                result = AllChem.EmbedMolecule(mol, randomSeed=42)
                if result == -1:
                    result = AllChem.EmbedMolecule(mol, useRandomCoords=True, maxAttempts=100, randomSeed=42)
                if result == -1:
                    params = AllChem.ETKDGv3()
                    params.randomSeed = 42
                    result = AllChem.EmbedMolecule(mol, params)
                
                if result != -1:
                    try:
                        AllChem.MMFFOptimizeMolecule(mol)
                    except: pass
                    
                    coords = []
                    conf = mol.GetConformer()
                    for i, atom in enumerate(mol.GetAtoms()):
                        pos = conf.GetAtomPosition(i)
                        coords.append((atom.GetSymbol(), pos.x, pos.y, pos.z))
                    return coords
        except ImportError:
            self.logger.warning("RDKit not installed, cannot generate from SMILES")
        except Exception as e:
            self.logger.warning(f"Failed to generate coordinates from SMILES: {e}")
            
        return None

    def _parse_pdb_coordinates(self, pdb_path: Path) -> Optional[List[Tuple[str, float, float, float]]]:
        """Parse PDB file"""
        coords = []
        encodings = ['utf-8', 'latin1', 'gbk', 'gb2312']
        
        for encoding in encodings:
            try:
                with open(pdb_path, 'r', encoding=encoding) as f:
                    for line in f:
                        if line.startswith(('ATOM', 'HETATM')):
                            try:
                                atom_name = line[12:16].strip()
                                element = ''
                                if len(line) > 76:
                                    element = line[76:78].strip()
                                    element = ''.join(c for c in element if c.isalpha())
                                if not element:
                                    element = ''.join(c for c in atom_name if c.isalpha())
                                if not element:
                                    element = atom_name[0] if atom_name else 'C'
                                    
                                x = float(line[30:38])
                                y = float(line[38:46])
                                z = float(line[46:54])
                                coords.append((element, x, y, z))
                            except: continue
                if coords: return coords
            except UnicodeDecodeError: continue
            except Exception: continue
            
        return None

    def _validate_and_correct_spin(self, smiles: str, charge: int, spin: int, coords: list = None) -> Tuple[int, int]:
        """Validate spin multiplicity"""
        try:
            from rdkit import Chem
            
            if coords:
                total_electrons = 0
                atomic_numbers = {
                    'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'P': 15, 'S': 16, 'Cl': 17,
                    'Li': 3, 'Na': 11, 'K': 19, 'Mg': 12, 'Ca': 20, 'Al': 13, 'Si': 14,
                    'B': 5, 'Br': 35, 'I': 53
                }
                for atom_symbol, _, _, _ in coords:
                    total_electrons += atomic_numbers.get(atom_symbol, 6) # Default to C if unknown
                total_electrons -= charge
                
                correct_spin = 1 if total_electrons % 2 == 0 else 2
                if correct_spin != spin:
                    self.logger.warning(f"Correcting spin from {spin} to {correct_spin} based on coordinates")
                    return charge, correct_spin
            
            elif smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    mol = Chem.AddHs(mol)
                    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
                    total_electrons = sum(atom.GetAtomicNum() for atom in mol.GetAtoms())
                    total_electrons -= total_charge
                    
                    num_radicals = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
                    if num_radicals > 0:
                        correct_spin = num_radicals + 1
                    else:
                        correct_spin = 1 if total_electrons % 2 == 0 else 2
                        
                    if correct_spin != spin or total_charge != charge:
                        self.logger.warning(f"Correcting charge/spin: {charge}/{spin} -> {total_charge}/{correct_spin}")
                        return total_charge, correct_spin
                        
        except Exception as e:
            self.logger.warning(f"Failed to validate spin: {e}")
            
        return charge, spin
