"""
Molecular structure tagging and feature detection for RSNet.

This module provides functions to identify structural features in molecules
that are relevant for reaction operator activation and site identification.
"""

from typing import Dict, List, Set, Tuple, Optional
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
import numpy as np


class StructureAnnotator:
    """结构标注器 - 分析分子结构特征"""

    def __init__(self):
        """初始化结构标注器"""
        pass

    def annotate_molecule(self, molecule) -> Dict[str, any]:
        """
        对分子进行结构标注

        Args:
            molecule: RSNet Molecule对象

        Returns:
            结构标签字典
        """
        return get_structure_tags(molecule)


def get_structure_tags(molecules) -> Dict[str, any]:
    """
    Get comprehensive structure tags for molecules.

    Args:
        molecules: RSNet Molecule object or list of Molecule objects

    Returns:
        Dictionary of structure tags and their values
    """
    # Handle both single molecule and list of molecules
    if isinstance(molecules, list):
        if len(molecules) == 0:
            return {}
        molecule = molecules[0]  # Use first molecule for now
    else:
        molecule = molecules

    mol = molecule.rdkit_mol
    tags = {}
    
    # Basic molecular properties
    tags['molecular_weight'] = Descriptors.MolWt(mol)
    tags['num_atoms'] = mol.GetNumAtoms()
    tags['num_bonds'] = mol.GetNumBonds()
    tags['num_rings'] = rdMolDescriptors.CalcNumRings(mol)
    
    # Ring analysis
    tags['small_rings'] = detect_small_rings(mol)
    tags['aromatic_rings'] = detect_aromatic_rings(mol)
    tags['saturated_rings'] = detect_saturated_rings(mol)
    
    # Bond analysis
    tags['polar_bonds'] = detect_polar_bonds(mol)
    tags['weak_bonds'] = detect_weak_bonds(mol)
    tags['multiple_bonds'] = detect_multiple_bonds(mol)
    
    # Atom analysis
    tags['heteroatoms'] = detect_heteroatoms(mol)
    tags['lewis_acid_sites'] = detect_lewis_acid_sites(mol)
    tags['lewis_base_sites'] = detect_lewis_base_sites(mol)
    
    # Electronic structure
    tags['pi_systems'] = detect_pi_systems(mol)
    tags['conjugated_systems'] = detect_conjugated_systems(mol)
    
    # Functional groups
    tags['functional_groups'] = detect_functional_groups(mol)
    
    # Reactivity indicators
    tags['electrophilic_sites'] = detect_electrophilic_sites(mol)
    tags['nucleophilic_sites'] = detect_nucleophilic_sites(mol)
    tags['radical_sites'] = detect_radical_sites(mol)
    
    return tags


def detect_small_rings(mol: Chem.Mol) -> List[Dict]:
    """Detect small rings (3-6 membered) in the molecule."""
    rings = []
    ring_info = mol.GetRingInfo()
    
    for ring in ring_info.AtomRings():
        if len(ring) <= 6:  # Small rings
            ring_data = {
                'size': len(ring),
                'atoms': list(ring),
                'is_aromatic': all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring),
                'strain_energy': estimate_ring_strain(len(ring))
            }
            rings.append(ring_data)
    
    return rings


def detect_polar_bonds(mol: Chem.Mol) -> List[Dict]:
    """Detect polar bonds based on electronegativity differences."""
    polar_bonds = []
    
    # Electronegativity values (Pauling scale)
    electronegativity = {
        'H': 2.20, 'C': 2.55, 'N': 3.04, 'O': 3.44, 'F': 3.98,
        'P': 2.19, 'S': 2.58, 'Cl': 3.16, 'Br': 2.96, 'I': 2.66,
        'Li': 0.98, 'Na': 0.93, 'K': 0.82, 'Mg': 1.31, 'Ca': 1.00
    }
    
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        
        symbol1 = atom1.GetSymbol()
        symbol2 = atom2.GetSymbol()
        
        if symbol1 in electronegativity and symbol2 in electronegativity:
            en_diff = abs(electronegativity[symbol1] - electronegativity[symbol2])
            
            if en_diff > 0.5:  # Threshold for polar bond
                polar_bonds.append({
                    'bond_idx': bond.GetIdx(),
                    'atoms': (atom1.GetIdx(), atom2.GetIdx()),
                    'symbols': (symbol1, symbol2),
                    'electronegativity_diff': en_diff,
                    'bond_type': bond.GetBondType(),
                    'polarity': 'high' if en_diff > 1.5 else 'medium'
                })
    
    return polar_bonds


def detect_heteroatoms(mol: Chem.Mol) -> List[Dict]:
    """Detect heteroatoms (non-carbon, non-hydrogen atoms)."""
    heteroatoms = []
    
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ['C', 'H']:
            heteroatoms.append({
                'atom_idx': atom.GetIdx(),
                'symbol': atom.GetSymbol(),
                'atomic_number': atom.GetAtomicNum(),
                'formal_charge': atom.GetFormalCharge(),
                'hybridization': atom.GetHybridization(),
                'num_neighbors': atom.GetDegree(),
                'is_aromatic': atom.GetIsAromatic(),
                'lone_pairs': estimate_lone_pairs(atom)
            })
    
    return heteroatoms


def detect_pi_systems(mol: Chem.Mol) -> List[Dict]:
    """Detect π-electron systems in the molecule."""
    pi_systems = []
    
    # Find aromatic systems
    aromatic_atoms = set()
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            aromatic_atoms.add(atom.GetIdx())
    
    if aromatic_atoms:
        pi_systems.append({
            'type': 'aromatic',
            'atoms': list(aromatic_atoms),
            'size': len(aromatic_atoms)
        })
    
    # Find double bonds
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            pi_systems.append({
                'type': 'double_bond',
                'atoms': [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()],
                'size': 2
            })
    
    # Find triple bonds
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.TRIPLE:
            pi_systems.append({
                'type': 'triple_bond',
                'atoms': [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()],
                'size': 2
            })
    
    return pi_systems


def detect_weak_bonds(mol: Chem.Mol) -> List[Dict]:
    """Detect potentially weak bonds that might break easily."""
    weak_bonds = []
    
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        
        # Check for weak bond types
        is_weak = False
        weakness_reason = []
        
        # Long bonds (C-C single bonds > 1.6 Å are considered weak)
        if bond.GetBondType() == Chem.BondType.SINGLE:
            # Estimate bond length based on atom types
            if (atom1.GetSymbol() == 'C' and atom2.GetSymbol() == 'C' and
                (atom1.GetDegree() > 3 or atom2.GetDegree() > 3)):
                is_weak = True
                weakness_reason.append('steric_strain')
        
        # Bonds to heteroatoms with lone pairs
        if atom1.GetSymbol() in ['O', 'N', 'S'] or atom2.GetSymbol() in ['O', 'N', 'S']:
            if bond.GetBondType() == Chem.BondType.SINGLE:
                is_weak = True
                weakness_reason.append('heteroatom_single_bond')
        
        # Bonds in small rings (high strain)
        ring_info = mol.GetRingInfo()
        if ring_info.NumBondRings(bond.GetIdx()) > 0:
            for ring in ring_info.BondRings():
                if bond.GetIdx() in ring and len(ring) <= 4:
                    is_weak = True
                    weakness_reason.append('ring_strain')
        
        if is_weak:
            weak_bonds.append({
                'bond_idx': bond.GetIdx(),
                'atoms': (atom1.GetIdx(), atom2.GetIdx()),
                'bond_type': bond.GetBondType(),
                'weakness_reasons': weakness_reason
            })
    
    return weak_bonds


def detect_lewis_acid_sites(mol: Chem.Mol) -> List[Dict]:
    """Detect potential Lewis acid sites (electron-deficient atoms)."""
    lewis_acids = []
    
    for atom in mol.GetAtoms():
        is_lewis_acid = False
        reasons = []
        
        # Atoms with positive formal charge
        if atom.GetFormalCharge() > 0:
            is_lewis_acid = True
            reasons.append('positive_charge')
        
        # Electron-deficient atoms (less than octet)
        if atom.GetSymbol() in ['B', 'Al']:  # Typical Lewis acids
            is_lewis_acid = True
            reasons.append('electron_deficient_element')
        
        # Transition metals (can act as Lewis acids)
        if atom.GetAtomicNum() in range(21, 31):  # 3d transition metals
            is_lewis_acid = True
            reasons.append('transition_metal')
        
        if is_lewis_acid:
            lewis_acids.append({
                'atom_idx': atom.GetIdx(),
                'symbol': atom.GetSymbol(),
                'formal_charge': atom.GetFormalCharge(),
                'reasons': reasons
            })
    
    return lewis_acids


def detect_functional_groups(mol: Chem.Mol) -> List[Dict]:
    """Detect common functional groups."""
    functional_groups = []
    
    # Define SMARTS patterns for common functional groups
    patterns = {
        'alcohol': '[OH]',
        'aldehyde': '[CX3H1](=O)',
        'ketone': '[CX3](=O)[C,c]',
        'carboxylic_acid': '[CX3](=O)[OX2H1]',
        'ester': '[CX3](=O)[OX2][C,c]',
        'ether': '[OX2]([C,c])[C,c]',
        'amine': '[NX3;H2,H1,H0]',
        'amide': '[NX3][CX3](=O)',
        'nitrile': '[CX2]#N',
        'nitro': '[NX3+](=O)[O-]',
        'sulfoxide': '[SX3](=O)',
        'sulfone': '[SX4](=O)(=O)',
        'phosphate': '[PX4](=O)([O,o])([O,o])[O,o]'
    }
    
    for group_name, smarts in patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern:
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                functional_groups.append({
                    'name': group_name,
                    'atoms': list(match),
                    'smarts': smarts
                })
    
    return functional_groups


# Helper functions
def estimate_ring_strain(ring_size: int) -> float:
    """Estimate ring strain energy in kcal/mol."""
    strain_energies = {3: 27.5, 4: 26.3, 5: 6.2, 6: 0.1, 7: 6.2, 8: 9.7}
    return strain_energies.get(ring_size, 0.0)


def estimate_lone_pairs(atom: Chem.Atom) -> int:
    """Estimate number of lone pairs on an atom."""
    valence_electrons = {
        'N': 5, 'O': 6, 'F': 7, 'P': 5, 'S': 6, 'Cl': 7, 'Br': 7, 'I': 7
    }
    
    symbol = atom.GetSymbol()
    if symbol not in valence_electrons:
        return 0
    
    total_electrons = valence_electrons[symbol]
    bonding_electrons = atom.GetTotalValence()
    formal_charge = atom.GetFormalCharge()
    
    lone_pair_electrons = total_electrons - bonding_electrons + formal_charge
    return max(0, lone_pair_electrons // 2)


def detect_aromatic_rings(mol: Chem.Mol) -> List[Dict]:
    """Detect aromatic ring systems."""
    aromatic_rings = []
    ring_info = mol.GetRingInfo()
    
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            aromatic_rings.append({
                'size': len(ring),
                'atoms': list(ring)
            })
    
    return aromatic_rings


def detect_saturated_rings(mol: Chem.Mol) -> List[Dict]:
    """Detect saturated (non-aromatic) ring systems."""
    saturated_rings = []
    ring_info = mol.GetRingInfo()
    
    for ring in ring_info.AtomRings():
        if not any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            saturated_rings.append({
                'size': len(ring),
                'atoms': list(ring)
            })
    
    return saturated_rings


def detect_multiple_bonds(mol: Chem.Mol) -> List[Dict]:
    """Detect double and triple bonds."""
    multiple_bonds = []
    
    for bond in mol.GetBonds():
        if bond.GetBondType() in [Chem.BondType.DOUBLE, Chem.BondType.TRIPLE]:
            multiple_bonds.append({
                'bond_idx': bond.GetIdx(),
                'atoms': (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
                'bond_type': str(bond.GetBondType())
            })
    
    return multiple_bonds


def detect_lewis_base_sites(mol: Chem.Mol) -> List[Dict]:
    """Detect potential Lewis base sites (electron-rich atoms)."""
    lewis_bases = []
    
    for atom in mol.GetAtoms():
        is_lewis_base = False
        reasons = []
        
        # Atoms with negative formal charge
        if atom.GetFormalCharge() < 0:
            is_lewis_base = True
            reasons.append('negative_charge')
        
        # Atoms with lone pairs
        if atom.GetSymbol() in ['N', 'O', 'S', 'P'] and estimate_lone_pairs(atom) > 0:
            is_lewis_base = True
            reasons.append('lone_pairs')
        
        if is_lewis_base:
            lewis_bases.append({
                'atom_idx': atom.GetIdx(),
                'symbol': atom.GetSymbol(),
                'formal_charge': atom.GetFormalCharge(),
                'lone_pairs': estimate_lone_pairs(atom),
                'reasons': reasons
            })
    
    return lewis_bases


def detect_conjugated_systems(mol: Chem.Mol) -> List[Dict]:
    """Detect conjugated π-electron systems."""
    # This is a simplified implementation
    # A more sophisticated version would use graph algorithms
    conjugated_systems = []
    
    # Find aromatic systems (automatically conjugated)
    aromatic_atoms = set()
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            aromatic_atoms.add(atom.GetIdx())
    
    if aromatic_atoms:
        conjugated_systems.append({
            'type': 'aromatic',
            'atoms': list(aromatic_atoms),
            'size': len(aromatic_atoms)
        })
    
    return conjugated_systems


def detect_electrophilic_sites(mol: Chem.Mol) -> List[Dict]:
    """Detect potential electrophilic sites."""
    electrophilic_sites = []
    
    for atom in mol.GetAtoms():
        is_electrophilic = False
        reasons = []
        
        # Positive formal charge
        if atom.GetFormalCharge() > 0:
            is_electrophilic = True
            reasons.append('positive_charge')
        
        # Carbonyl carbon
        if atom.GetSymbol() == 'C':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O':
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        is_electrophilic = True
                        reasons.append('carbonyl_carbon')
        
        if is_electrophilic:
            electrophilic_sites.append({
                'atom_idx': atom.GetIdx(),
                'symbol': atom.GetSymbol(),
                'reasons': reasons
            })
    
    return electrophilic_sites


def detect_nucleophilic_sites(mol: Chem.Mol) -> List[Dict]:
    """Detect potential nucleophilic sites."""
    nucleophilic_sites = []
    
    for atom in mol.GetAtoms():
        is_nucleophilic = False
        reasons = []
        
        # Negative formal charge
        if atom.GetFormalCharge() < 0:
            is_nucleophilic = True
            reasons.append('negative_charge')
        
        # Atoms with lone pairs
        if estimate_lone_pairs(atom) > 0:
            is_nucleophilic = True
            reasons.append('lone_pairs')
        
        if is_nucleophilic:
            nucleophilic_sites.append({
                'atom_idx': atom.GetIdx(),
                'symbol': atom.GetSymbol(),
                'lone_pairs': estimate_lone_pairs(atom),
                'reasons': reasons
            })
    
    return nucleophilic_sites


def detect_radical_sites(mol: Chem.Mol) -> List[Dict]:
    """Detect potential radical sites (atoms with unpaired electrons)."""
    radical_sites = []
    
    for atom in mol.GetAtoms():
        # Check for radical character based on valence
        num_radical_electrons = atom.GetNumRadicalElectrons()
        
        if num_radical_electrons > 0:
            radical_sites.append({
                'atom_idx': atom.GetIdx(),
                'symbol': atom.GetSymbol(),
                'num_radical_electrons': num_radical_electrons
            })
    
    return radical_sites
