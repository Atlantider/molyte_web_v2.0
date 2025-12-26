"""
Tests for hydrogen transfer operator.
"""

import pytest
from rdkit import Chem

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.operators.hydrogen_transfer import HydrogenTransferOperator, HydrogenTransferTemplate
from rsnet.utils.structure_analysis import find_hydrogen_transfer_sites, get_structure_tags


class TestHydrogenTransferOperator:
    """Test cases for HydrogenTransferOperator."""
    
    def test_operator_initialization(self):
        """Test operator initialization."""
        op = HydrogenTransferOperator()
        assert op.name == "HydrogenTransfer"
        assert op.enabled is True
        assert isinstance(op.get_templates(), list)
    
    def test_operator_with_config(self):
        """Test operator initialization with config."""
        config = {
            'enabled': True,
            'max_distance': 3.5,
            'min_distance': 1.2,
            'energy_threshold': 40.0
        }
        op = HydrogenTransferOperator(config=config)
        assert op.max_distance == 3.5
        assert op.min_distance == 1.2
        assert op.energy_threshold == 40.0
    
    def test_is_applicable_simple_alcohol(self):
        """Test applicability check with simple alcohol."""
        mol = Molecule.from_smiles("CCO", name="ethanol")  # Has O-H
        op = HydrogenTransferOperator()
        
        assert op.is_applicable([mol]) is True
    
    def test_is_applicable_hydrocarbon(self):
        """Test applicability check with hydrocarbon (should be False)."""
        mol = Molecule.from_smiles("CC", name="ethane")  # No heteroatoms
        op = HydrogenTransferOperator()
        
        assert op.is_applicable([mol]) is False
    
    def test_is_applicable_multiple_molecules(self):
        """Test applicability check with multiple molecules (should be False for now)."""
        mol1 = Molecule.from_smiles("CCO", name="ethanol")
        mol2 = Molecule.from_smiles("CO", name="methanol")
        op = HydrogenTransferOperator()
        
        assert op.is_applicable([mol1, mol2]) is False
    
    def test_structure_analysis_ethanol(self):
        """Test structure analysis on ethanol."""
        mol = Molecule.from_smiles("CCO", name="ethanol")
        tags = get_structure_tags(mol)
        
        assert tags['has_heteroatoms'] is True
        assert tags['has_acidic_hydrogens'] is True
        assert tags['has_polar_bonds'] is True
        assert len(tags['heteroatoms']) == 1  # One oxygen
        assert len(tags['acidic_hydrogens']) == 1  # One O-H hydrogen
    
    def test_find_hydrogen_transfer_sites_ethanol(self):
        """Test finding H-transfer sites in ethanol."""
        mol = Molecule.from_smiles("CCO", name="ethanol")
        sites = find_hydrogen_transfer_sites(mol)
        
        # Ethanol should have at least one potential site
        # (though intramolecular transfer might not be favorable)
        assert isinstance(sites, list)
        # Note: The exact number depends on the molecule's 3D structure
    
    def test_apply_operator_ethanol(self):
        """Test applying operator to ethanol."""
        mol = Molecule.from_smiles("CCO", name="ethanol")
        env = Environment(temperature=298.15)
        op = HydrogenTransferOperator()
        
        reactions = op.apply([mol], env)
        
        # Should return a list (might be empty if no favorable transfers)
        assert isinstance(reactions, list)
    
    def test_apply_operator_carboxylic_acid(self):
        """Test applying operator to carboxylic acid."""
        mol = Molecule.from_smiles("CC(=O)O", name="acetic_acid")
        env = Environment(temperature=298.15)
        op = HydrogenTransferOperator()
        
        reactions = op.apply([mol], env)
        
        # Should return a list
        assert isinstance(reactions, list)
        
        # Print results for debugging
        print(f"Found {len(reactions)} reactions for acetic acid")
        for i, rxn in enumerate(reactions):
            print(f"  Reaction {i+1}: {rxn.name}")
    
    def test_template_initialization(self):
        """Test hydrogen transfer template initialization."""
        template = HydrogenTransferTemplate(
            donor_heavy_idx=1,
            donor_hydrogen_idx=2,
            acceptor_idx=3,
            name="test_template"
        )
        
        assert template.name == "test_template"
        assert template.donor_heavy_idx == 1
        assert template.donor_hydrogen_idx == 2
        assert template.acceptor_idx == 3
    
    def test_template_auto_naming(self):
        """Test automatic template naming."""
        template = HydrogenTransferTemplate(
            donor_heavy_idx=1,
            donor_hydrogen_idx=2,
            acceptor_idx=3
        )
        
        assert "H_transfer_1_3" in template.name
    
    def test_required_drives(self):
        """Test getting required driving forces."""
        op = HydrogenTransferOperator()
        drives = op.get_required_drives()
        
        assert isinstance(drives, list)
        assert 'thermal' in drives or 'polar_environment' in drives
    
    def test_structure_conditions(self):
        """Test structure condition checking."""
        op = HydrogenTransferOperator()
        
        # Test with ethanol tags
        mol = Molecule.from_smiles("CCO", name="ethanol")
        tags = get_structure_tags(mol)
        assert op.check_structure_conditions(tags) is True
        
        # Test with ethane tags (no heteroatoms)
        mol = Molecule.from_smiles("CC", name="ethane")
        tags = get_structure_tags(mol)
        assert op.check_structure_conditions(tags) is False


class TestStructureAnalysis:
    """Test cases for structure analysis functions."""
    
    def test_detect_heteroatoms_ethanol(self):
        """Test heteroatom detection in ethanol."""
        mol = Molecule.from_smiles("CCO", name="ethanol")
        tags = get_structure_tags(mol)
        
        heteroatoms = tags['heteroatoms']
        assert len(heteroatoms) == 1
        
        # Check that the heteroatom is oxygen
        rdkit_mol = mol.rdkit_mol
        oxygen_atom = rdkit_mol.GetAtomWithIdx(heteroatoms[0])
        assert oxygen_atom.GetSymbol() == 'O'
    
    def test_detect_acidic_hydrogens_ethanol(self):
        """Test acidic hydrogen detection in ethanol."""
        mol = Molecule.from_smiles("CCO", name="ethanol")
        tags = get_structure_tags(mol)
        
        acidic_h = tags['acidic_hydrogens']
        assert len(acidic_h) == 1
        
        # Check that the hydrogen is bonded to oxygen
        rdkit_mol = mol.rdkit_mol
        h_atom = rdkit_mol.GetAtomWithIdx(acidic_h[0])
        assert h_atom.GetSymbol() == 'H'
        
        # Check its neighbor is oxygen
        neighbors = [n for n in h_atom.GetNeighbors()]
        assert len(neighbors) == 1
        assert neighbors[0].GetSymbol() == 'O'
    
    def test_detect_polar_bonds_ethanol(self):
        """Test polar bond detection in ethanol."""
        mol = Molecule.from_smiles("CCO", name="ethanol")
        tags = get_structure_tags(mol)
        
        polar_bonds = tags['polar_bonds']
        assert len(polar_bonds) >= 1  # At least C-O and O-H bonds
    
    def test_structure_tags_comprehensive(self):
        """Test comprehensive structure tagging."""
        mol = Molecule.from_smiles("CC(=O)O", name="acetic_acid")
        tags = get_structure_tags(mol)
        
        # Check all expected keys are present
        expected_keys = [
            'heteroatoms', 'polar_bonds', 'acidic_hydrogens', 'h_bond_acceptors',
            'weak_bonds', 'small_rings', 'pi_systems', 'molecular_weight',
            'num_rotatable_bonds', 'num_hbd', 'num_hba', 'logp', 'tpsa',
            'has_heteroatoms', 'has_polar_bonds', 'has_acidic_hydrogens',
            'has_small_rings', 'has_pi_systems', 'has_weak_bonds'
        ]
        
        for key in expected_keys:
            assert key in tags
        
        # Check boolean flags
        assert tags['has_heteroatoms'] is True
        assert tags['has_polar_bonds'] is True
        assert tags['has_acidic_hydrogens'] is True
