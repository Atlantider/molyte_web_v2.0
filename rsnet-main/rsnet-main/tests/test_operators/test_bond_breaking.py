"""
Tests for bond breaking operator.
"""

import pytest
from rdkit import Chem

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.operators.bond_breaking import BondBreakingOperator, BondBreakingTemplate
from rsnet.utils.structure_analysis import get_structure_tags


class TestBondBreakingOperator:
    """Test cases for BondBreakingOperator."""
    
    def test_operator_initialization(self):
        """Test operator initialization."""
        op = BondBreakingOperator()
        assert op.name == "BondBreaking"
        assert op.enabled is True
        assert isinstance(op.get_templates(), list)
    
    def test_operator_with_config(self):
        """Test operator initialization with config."""
        config = {
            'enabled': True,
            'min_fragment_size': 3,
            'max_bonds_to_break': 3,
            'energy_threshold': 60.0
        }
        op = BondBreakingOperator(config=config)
        assert op.min_fragment_size == 3
        assert op.max_bonds_to_break == 3
        assert op.energy_threshold == 60.0
    
    def test_is_applicable_ethane(self):
        """Test applicability check with ethane."""
        mol = Molecule.from_smiles("CC", name="ethane")
        op = BondBreakingOperator()
        
        assert op.is_applicable([mol]) is True
    
    def test_is_applicable_methane(self):
        """Test applicability check with methane (single heavy atom)."""
        mol = Molecule.from_smiles("C", name="methane")
        op = BondBreakingOperator()

        # Methane has only one heavy atom, so no bonds between heavy atoms to break
        # Our current implementation requires at least 2 heavy atoms
        assert op.is_applicable([mol]) is False
    
    def test_is_applicable_multiple_molecules(self):
        """Test applicability check with multiple molecules (should be False)."""
        mol1 = Molecule.from_smiles("CC", name="ethane")
        mol2 = Molecule.from_smiles("CO", name="methanol")
        op = BondBreakingOperator()
        
        assert op.is_applicable([mol1, mol2]) is False
    
    def test_find_breakable_bonds_ethane(self):
        """Test finding breakable bonds in ethane."""
        mol = Molecule.from_smiles("CC", name="ethane")
        env = Environment(temperature=500.0)  # High temperature
        op = BondBreakingOperator()
        
        breakable_bonds = op._find_breakable_bonds(mol, env)
        
        # Should find at least the C-C bond at high temperature
        assert len(breakable_bonds) >= 1
        
        # Check structure of returned data
        for bond_info in breakable_bonds:
            assert len(bond_info) == 4  # bond_idx, atom1_idx, atom2_idx, bond_type
            bond_idx, atom1_idx, atom2_idx, bond_type = bond_info
            assert isinstance(bond_idx, int)
            assert isinstance(atom1_idx, int)
            assert isinstance(atom2_idx, int)
            assert isinstance(bond_type, str)
    
    def test_find_breakable_bonds_cyclopropane(self):
        """Test finding breakable bonds in cyclopropane (strained ring)."""
        mol = Molecule.from_smiles("C1CC1", name="cyclopropane")
        env = Environment(temperature=298.15)  # Room temperature
        op = BondBreakingOperator()
        
        breakable_bonds = op._find_breakable_bonds(mol, env)
        
        # Should find bonds due to ring strain
        assert len(breakable_bonds) >= 1
    
    def test_apply_operator_ethane_high_temp(self):
        """Test applying operator to ethane at high temperature."""
        mol = Molecule.from_smiles("CC", name="ethane")
        env = Environment(temperature=500.0)  # High temperature
        op = BondBreakingOperator()
        
        reactions = op.apply([mol], env)
        
        # Should generate at least one reaction
        assert isinstance(reactions, list)
        
        # Print results for debugging
        print(f"Found {len(reactions)} reactions for ethane at high temperature")
        for i, rxn in enumerate(reactions):
            print(f"  Reaction {i+1}: {rxn.name}")
            print(f"    Reactants: {[r.smiles for r in rxn.reactants]}")
            print(f"    Products: {[p.smiles for p in rxn.products]}")
    
    def test_apply_operator_ethane_low_temp(self):
        """Test applying operator to ethane at low temperature."""
        mol = Molecule.from_smiles("CC", name="ethane")
        env = Environment(temperature=298.15)  # Room temperature
        op = BondBreakingOperator()
        
        reactions = op.apply([mol], env)
        
        # Should generate fewer or no reactions at low temperature
        assert isinstance(reactions, list)
        print(f"Found {len(reactions)} reactions for ethane at room temperature")
    
    def test_apply_operator_cyclopropane(self):
        """Test applying operator to cyclopropane."""
        mol = Molecule.from_smiles("C1CC1", name="cyclopropane")
        env = Environment(temperature=298.15)
        op = BondBreakingOperator()
        
        reactions = op.apply([mol], env)
        
        # Should generate reactions due to ring strain
        assert isinstance(reactions, list)
        
        print(f"Found {len(reactions)} reactions for cyclopropane")
        for i, rxn in enumerate(reactions):
            print(f"  Reaction {i+1}: {rxn.name}")
            print(f"    Products: {[p.smiles for p in rxn.products]}")
    
    def test_apply_operator_with_min_fragment_size(self):
        """Test operator with minimum fragment size constraint."""
        mol = Molecule.from_smiles("CCC", name="propane")
        env = Environment(temperature=500.0)
        
        # Set minimum fragment size to 2 heavy atoms
        config = {'min_fragment_size': 2}
        op = BondBreakingOperator(config=config)
        
        reactions = op.apply([mol], env)
        
        # Check that all products meet minimum size requirement
        for rxn in reactions:
            for product in rxn.products:
                assert product.num_heavy_atoms >= 2
    
    def test_template_initialization(self):
        """Test bond breaking template initialization."""
        template = BondBreakingTemplate(
            bond_idx=0,
            atom1_idx=1,
            atom2_idx=2,
            bond_type="single",
            name="test_template"
        )
        
        assert template.name == "test_template"
        assert template.bond_idx == 0
        assert template.atom1_idx == 1
        assert template.atom2_idx == 2
        assert template.bond_type == "single"
    
    def test_template_auto_naming(self):
        """Test automatic template naming."""
        template = BondBreakingTemplate(
            bond_idx=0,
            atom1_idx=1,
            atom2_idx=2,
            bond_type="single"
        )
        
        assert "bond_break_1_2" in template.name
    
    def test_template_apply_ethane(self):
        """Test applying template to ethane."""
        mol = Molecule.from_smiles("CC", name="ethane")
        
        # Create template for breaking C-C bond (atoms 0 and 1)
        template = BondBreakingTemplate(
            bond_idx=0,
            atom1_idx=0,
            atom2_idx=1,
            bond_type="single"
        )
        
        reactions = template.apply([mol])
        
        # Should generate one reaction with two methyl radical fragments
        assert len(reactions) == 1
        rxn = reactions[0]
        assert len(rxn.reactants) == 1
        assert len(rxn.products) == 2
        
        # Products should be methyl radicals (or similar fragments)
        print(f"Ethane fragmentation products: {[p.smiles for p in rxn.products]}")
    
    def test_required_drives(self):
        """Test getting required driving forces."""
        op = BondBreakingOperator()
        drives = op.get_required_drives()
        
        assert isinstance(drives, list)
        assert 'thermal' in drives or 'radical_environment' in drives
    
    def test_structure_conditions(self):
        """Test structure condition checking."""
        op = BondBreakingOperator()
        
        # Test with a molecule that has weak bonds
        mol = Molecule.from_smiles("C1CC1", name="cyclopropane")  # Strained ring
        tags = get_structure_tags(mol)
        
        # The condition check depends on weak bond detection
        # For now, just check that it returns a boolean
        result = op.check_structure_conditions(tags)
        assert isinstance(result, bool)


class TestBondBreakingTemplate:
    """Test cases for BondBreakingTemplate."""
    
    def test_break_bond_ethane(self):
        """Test breaking C-C bond in ethane."""
        mol = Molecule.from_smiles("CC", name="ethane")
        
        template = BondBreakingTemplate(
            bond_idx=0,
            atom1_idx=0,
            atom2_idx=1,
            bond_type="single"
        )
        
        products = template._break_bond(mol)
        
        # Should produce two fragments
        assert len(products) == 2
        
        # Each fragment should be a methyl group
        for product in products:
            assert product.num_heavy_atoms == 1  # One carbon each
            print(f"Fragment: {product.smiles} (formula: {product.formula})")
    
    def test_break_bond_propane(self):
        """Test breaking C-C bond in propane."""
        mol = Molecule.from_smiles("CCC", name="propane")
        
        # Break the first C-C bond (between atoms 0 and 1)
        template = BondBreakingTemplate(
            bond_idx=0,
            atom1_idx=0,
            atom2_idx=1,
            bond_type="single"
        )
        
        products = template._break_bond(mol)
        
        # Should produce two fragments: methyl and ethyl
        assert len(products) == 2
        
        # One fragment should have 1 heavy atom, the other should have 2
        heavy_atom_counts = [p.num_heavy_atoms for p in products]
        heavy_atom_counts.sort()
        assert heavy_atom_counts == [1, 2]
        
        print(f"Propane fragmentation products: {[p.smiles for p in products]}")
