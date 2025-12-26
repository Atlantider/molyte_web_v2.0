"""
Molecule class for RSNet.

This module provides the core Molecule class that wraps RDKit functionality
and adds additional features needed for reaction network generation.
"""

import hashlib
from typing import Optional, List, Dict, Any
import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
    from rdkit.Chem.rdchem import Mol as RDKitMol
except ImportError:
    raise ImportError("RDKit is required but not installed. Please install RDKit.")


class Molecule:
    """
    A molecule representation that wraps RDKit Mol objects with additional functionality.
    
    This class provides methods for molecular structure analysis, property calculation,
    and integration with reaction operators.
    """
    
    def __init__(self, rdkit_mol: RDKitMol, name: Optional[str] = None):
        """
        Initialize a Molecule from an RDKit Mol object.
        
        Args:
            rdkit_mol: RDKit Mol object
            name: Optional name for the molecule
        """
        if rdkit_mol is None:
            raise ValueError("RDKit molecule cannot be None")
        
        self._mol = rdkit_mol
        self._name = name
        self._properties = {}
        self._hash = None
        
        # Ensure molecule has explicit hydrogens for accurate analysis
        self._mol = Chem.AddHs(self._mol)
    
    @classmethod
    def from_smiles(cls, smiles: str, name: Optional[str] = None) -> 'Molecule':
        """
        Create a Molecule from a SMILES string.
        
        Args:
            smiles: SMILES representation of the molecule
            name: Optional name for the molecule
            
        Returns:
            Molecule object
            
        Raises:
            ValueError: If SMILES string is invalid
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        
        return cls(mol, name or smiles)
    
    @classmethod
    def from_mol_file(cls, filepath: str, name: Optional[str] = None) -> 'Molecule':
        """
        Create a Molecule from a MOL file.
        
        Args:
            filepath: Path to the MOL file
            name: Optional name for the molecule
            
        Returns:
            Molecule object
        """
        mol = Chem.MolFromMolFile(filepath)
        if mol is None:
            raise ValueError(f"Could not read molecule from file: {filepath}")
        
        return cls(mol, name)
    
    @property
    def rdkit_mol(self) -> RDKitMol:
        """Get the underlying RDKit Mol object."""
        return self._mol
    
    @property
    def name(self) -> str:
        """Get the molecule name."""
        return self._name or self.smiles
    
    @name.setter
    def name(self, value: str):
        """Set the molecule name."""
        self._name = value
    
    @property
    def smiles(self) -> str:
        """Get the canonical SMILES representation (without explicit hydrogens)."""
        # Remove explicit hydrogens for canonical SMILES representation
        mol_no_hs = Chem.RemoveHs(self._mol)
        return Chem.MolToSmiles(mol_no_hs)
    
    @property
    def formula(self) -> str:
        """Get the molecular formula."""
        return rdMolDescriptors.CalcMolFormula(self._mol)
    
    @property
    def molecular_weight(self) -> float:
        """Get the molecular weight in g/mol."""
        return Descriptors.MolWt(self._mol)
    
    @property
    def num_atoms(self) -> int:
        """Get the number of atoms (including hydrogens)."""
        return self._mol.GetNumAtoms()
    
    @property
    def num_heavy_atoms(self) -> int:
        """Get the number of heavy atoms (excluding hydrogens)."""
        return self._mol.GetNumHeavyAtoms()
    
    @property
    def num_bonds(self) -> int:
        """Get the number of bonds."""
        return self._mol.GetNumBonds()
    
    def get_hash(self) -> str:
        """
        Get a unique hash for this molecule based on its structure.
        
        Returns:
            MD5 hash string
        """
        if self._hash is None:
            # Use canonical SMILES for hashing to ensure consistency
            canonical_smiles = Chem.MolToSmiles(self._mol, canonical=True)
            self._hash = hashlib.md5(canonical_smiles.encode()).hexdigest()
        return self._hash
    
    def get_atoms(self) -> List[Dict[str, Any]]:
        """
        Get information about all atoms in the molecule.
        
        Returns:
            List of dictionaries containing atom information
        """
        atoms = []
        for atom in self._mol.GetAtoms():
            atoms.append({
                'idx': atom.GetIdx(),
                'symbol': atom.GetSymbol(),
                'atomic_num': atom.GetAtomicNum(),
                'formal_charge': atom.GetFormalCharge(),
                'hybridization': str(atom.GetHybridization()),
                'is_aromatic': atom.GetIsAromatic(),
                'num_explicit_hs': atom.GetNumExplicitHs(),
                'num_implicit_hs': atom.GetNumImplicitHs(),
            })
        return atoms
    
    def get_bonds(self) -> List[Dict[str, Any]]:
        """
        Get information about all bonds in the molecule.
        
        Returns:
            List of dictionaries containing bond information
        """
        bonds = []
        for bond in self._mol.GetBonds():
            bonds.append({
                'idx': bond.GetIdx(),
                'begin_atom': bond.GetBeginAtomIdx(),
                'end_atom': bond.GetEndAtomIdx(),
                'bond_type': str(bond.GetBondType()),
                'is_aromatic': bond.GetIsAromatic(),
                'is_conjugated': bond.GetIsConjugated(),
                'is_in_ring': bond.IsInRing(),
            })
        return bonds

    def has_conformer(self) -> bool:
        """
        Check if the molecule has a 3D conformer.

        Returns:
            True if molecule has at least one conformer, False otherwise
        """
        return self._mol.GetNumConformers() > 0

    def generate_conformer(self, num_confs: int = 1, random_seed: int = 42) -> bool:
        """
        Generate 3D conformer(s) for the molecule.
        
        Args:
            num_confs: Number of conformers to generate
            random_seed: Random seed for reproducibility
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Remove existing conformers
            self._mol.RemoveAllConformers()
            
            # Generate conformers
            conf_ids = AllChem.EmbedMultipleConfs(
                self._mol, 
                numConfs=num_confs,
                randomSeed=random_seed,
                useExpTorsionAnglePrefs=True,
                useBasicKnowledge=True
            )
            
            if len(conf_ids) == 0:
                return False
            
            # Optimize conformers with MMFF
            for conf_id in conf_ids:
                AllChem.MMFFOptimizeMolecule(self._mol, confId=conf_id)
            
            return True
            
        except Exception:
            return False
    
    def get_coordinates(self, conf_id: int = 0) -> Optional[np.ndarray]:
        """
        Get 3D coordinates of atoms.
        
        Args:
            conf_id: Conformer ID
            
        Returns:
            Numpy array of shape (n_atoms, 3) or None if no conformer
        """
        if self._mol.GetNumConformers() == 0:
            return None
        
        conf = self._mol.GetConformer(conf_id)
        coords = []
        for i in range(self._mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            coords.append([pos.x, pos.y, pos.z])
        
        return np.array(coords)
    
    def copy(self) -> 'Molecule':
        """
        Create a copy of this molecule.
        
        Returns:
            New Molecule object
        """
        mol_copy = Chem.Mol(self._mol)
        return Molecule(mol_copy, self._name)
    
    def __eq__(self, other) -> bool:
        """Check if two molecules are the same based on their hash."""
        if not isinstance(other, Molecule):
            return False
        return self.get_hash() == other.get_hash()
    
    def __hash__(self) -> int:
        """Hash function for use in sets and dictionaries."""
        return hash(self.get_hash())
    
    def __str__(self) -> str:
        """String representation of the molecule."""
        return f"Molecule(name='{self.name}', formula='{self.formula}', smiles='{self.smiles}')"
    
    def __repr__(self) -> str:
        """Detailed string representation."""
        return self.__str__()
