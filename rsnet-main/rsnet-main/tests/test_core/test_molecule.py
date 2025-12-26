"""
Tests for the Molecule class.
"""

import pytest
import numpy as np
from rdkit import Chem

from rsnet.core.molecule import Molecule


class TestMolecule:
    """Test cases for the Molecule class."""
    
    def test_from_smiles_valid(self):
        """Test creating molecule from valid SMILES."""
        mol = Molecule.from_smiles("CCO")  # ethanol
        assert mol.smiles == "CCO"
        assert mol.formula == "C2H6O"
        assert mol.name == "CCO"
    
    def test_from_smiles_invalid(self):
        """Test creating molecule from invalid SMILES."""
        with pytest.raises(ValueError, match="Invalid SMILES string"):
            Molecule.from_smiles("invalid_smiles")
    
    def test_from_smiles_with_name(self):
        """Test creating molecule with custom name."""
        mol = Molecule.from_smiles("CCO", name="ethanol")
        assert mol.name == "ethanol"
        assert mol.smiles == "CCO"
    
    def test_basic_properties(self):
        """Test basic molecular properties."""
        mol = Molecule.from_smiles("CCO")  # ethanol
        
        assert mol.formula == "C2H6O"
        assert mol.num_heavy_atoms == 3
        assert mol.num_atoms == 9  # including hydrogens
        assert mol.num_bonds == 8
        assert mol.molecular_weight > 46.0  # approximately 46.07
    
    def test_methane_properties(self):
        """Test properties of methane."""
        mol = Molecule.from_smiles("C")  # methane
        
        assert mol.formula == "CH4"
        assert mol.num_heavy_atoms == 1
        assert mol.num_atoms == 5  # 1 C + 4 H
        assert mol.num_bonds == 4
    
    def test_benzene_properties(self):
        """Test properties of benzene."""
        mol = Molecule.from_smiles("c1ccccc1")  # benzene
        
        assert mol.formula == "C6H6"
        assert mol.num_heavy_atoms == 6
        assert mol.num_atoms == 12  # 6 C + 6 H
    
    def test_hash_consistency(self):
        """Test that hash is consistent for same molecule."""
        mol1 = Molecule.from_smiles("CCO")
        mol2 = Molecule.from_smiles("CCO")
        
        assert mol1.get_hash() == mol2.get_hash()
        assert mol1 == mol2
        assert hash(mol1) == hash(mol2)
    
    def test_hash_different_molecules(self):
        """Test that different molecules have different hashes."""
        mol1 = Molecule.from_smiles("CCO")  # ethanol
        mol2 = Molecule.from_smiles("CCC")  # propane
        
        assert mol1.get_hash() != mol2.get_hash()
        assert mol1 != mol2
        assert hash(mol1) != hash(mol2)
    
    def test_get_atoms(self):
        """Test getting atom information."""
        mol = Molecule.from_smiles("CO")  # methanol
        atoms = mol.get_atoms()
        
        assert len(atoms) == 6  # C, O, 4H
        
        # Check carbon atom
        carbon_atoms = [a for a in atoms if a['symbol'] == 'C']
        assert len(carbon_atoms) == 1
        assert carbon_atoms[0]['atomic_num'] == 6
        
        # Check oxygen atom
        oxygen_atoms = [a for a in atoms if a['symbol'] == 'O']
        assert len(oxygen_atoms) == 1
        assert oxygen_atoms[0]['atomic_num'] == 8
    
    def test_get_bonds(self):
        """Test getting bond information."""
        mol = Molecule.from_smiles("CCO")  # ethanol
        bonds = mol.get_bonds()
        
        assert len(bonds) == 8  # C-C, C-O, and 6 C-H/O-H bonds
        
        # Check that we have single bonds
        single_bonds = [b for b in bonds if b['bond_type'] == 'SINGLE']
        assert len(single_bonds) == 8
    
    def test_conformer_generation(self):
        """Test 3D conformer generation."""
        mol = Molecule.from_smiles("CCO")  # ethanol
        
        # Initially no conformers
        assert mol.get_coordinates() is None
        
        # Generate conformer
        success = mol.generate_conformer()
        assert success is True
        
        # Now should have coordinates
        coords = mol.get_coordinates()
        assert coords is not None
        assert coords.shape == (9, 3)  # 9 atoms, 3D coordinates
        assert isinstance(coords, np.ndarray)
    
    def test_copy(self):
        """Test molecule copying."""
        mol1 = Molecule.from_smiles("CCO", name="original")
        mol2 = mol1.copy()
        
        # Should be equal but not the same object
        assert mol1 == mol2
        assert mol1 is not mol2
        assert mol1.rdkit_mol is not mol2.rdkit_mol
        
        # Names should be the same
        assert mol1.name == mol2.name
    
    def test_string_representations(self):
        """Test string representations."""
        mol = Molecule.from_smiles("CCO", name="ethanol")
        
        str_repr = str(mol)
        assert "ethanol" in str_repr
        assert "C2H6O" in str_repr
        assert "CCO" in str_repr
        
        repr_str = repr(mol)
        assert repr_str == str_repr
    
    def test_name_property(self):
        """Test name property getter and setter."""
        mol = Molecule.from_smiles("CCO")
        
        # Default name should be SMILES
        assert mol.name == "CCO"
        
        # Set custom name
        mol.name = "ethanol"
        assert mol.name == "ethanol"
    
    def test_rdkit_mol_property(self):
        """Test access to underlying RDKit molecule."""
        mol = Molecule.from_smiles("CCO")
        rdkit_mol = mol.rdkit_mol

        assert isinstance(rdkit_mol, Chem.Mol)
        # The underlying RDKit mol has explicit hydrogens, so we need to remove them for comparison
        mol_no_hs = Chem.RemoveHs(rdkit_mol)
        assert Chem.MolToSmiles(mol_no_hs) == "CCO"


if __name__ == "__main__":
    pytest.main([__file__])
