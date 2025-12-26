"""
Tests for enhanced RSNet features including physical driving forces,
structure annotation, and expanded operators.
"""

import pytest
from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.features.structure_tags import get_structure_tags
from rsnet.features.driving_forces import get_driving_forces
from rsnet.operators.registry import OPERATOR_REGISTRY
from rsnet.operators.redox import RedoxOperator
from rsnet.operators.addition import AdditionOperator
from rsnet.operators.rearrangement import RearrangementOperator


class TestEnhancedEnvironment:
    """Test enhanced environment with physical driving forces."""
    
    def test_basic_environment(self):
        """Test basic environment creation."""
        env = Environment(temperature=298.15, pressure=1.0)
        assert env.temperature == 298.15
        assert env.pressure == 1.0
    
    def test_electrochemical_environment(self):
        """Test electrochemical environment parameters."""
        env = Environment(
            temperature=298.15,
            electrode_type='cathode',
            voltage=4.2,
            li_activity=0.5,
            interface_type='SEI'
        )
        
        assert env.electrode_type == 'cathode'
        assert env.voltage == 4.2
        assert env.li_activity == 0.5
        assert env.interface_type == 'SEI'
    
    def test_active_drives(self):
        """Test active driving forces detection."""
        # High temperature environment
        env1 = Environment(temperature=800.0)
        drives1 = env1.get_active_drives()
        assert drives1['thermal'] == True
        assert drives1['high_temperature'] == True
        
        # Electrochemical environment
        env2 = Environment(electrode_type='cathode', voltage=4.5)
        drives2 = env2.get_active_drives()
        assert drives2['oxidation'] == True
        assert drives2['electrochemical'] == True
        
        # Reducing environment
        env3 = Environment(electrode_type='anode', voltage=0.1)
        drives3 = env3.get_active_drives()
        assert drives3['reduction'] == True


class TestStructureAnnotation:
    """Test structure annotation system."""
    
    def test_basic_structure_tags(self):
        """Test basic structure tag detection."""
        mol = Molecule.from_smiles('CCO')  # Ethanol
        tags = get_structure_tags([mol])
        
        assert 'molecular_weight' in tags
        assert 'num_atoms' in tags
        assert 'num_bonds' in tags
        assert tags['num_atoms'] == 9  # 3C + 6H + 1O
        assert tags['num_bonds'] == 8
    
    def test_functional_group_detection(self):
        """Test functional group detection."""
        # Alcohol
        mol1 = Molecule.from_smiles('CCO')
        tags1 = get_structure_tags([mol1])
        assert tags1.get('functional_groups', False)
        
        # Ketone
        mol2 = Molecule.from_smiles('CC(=O)C')
        tags2 = get_structure_tags([mol2])
        assert tags2.get('functional_groups', False)
    
    def test_ring_detection(self):
        """Test ring structure detection."""
        # Benzene
        mol1 = Molecule.from_smiles('c1ccccc1')
        tags1 = get_structure_tags([mol1])
        assert tags1.get('num_rings', 0) > 0
        assert tags1.get('aromatic_rings', False)
        
        # Cyclopropane (small ring)
        mol2 = Molecule.from_smiles('C1CC1')
        tags2 = get_structure_tags([mol2])
        assert tags2.get('small_rings', False)


class TestDrivingForces:
    """Test driving force evaluation system."""
    
    def test_thermal_driving_force(self):
        """Test thermal driving force evaluation."""
        mol = Molecule.from_smiles('CCO')
        
        # Room temperature
        env1 = Environment(temperature=298.15)
        forces1 = get_driving_forces([mol], env1)
        assert forces1['thermal'] == 0.0
        
        # High temperature
        env2 = Environment(temperature=500.0)
        forces2 = get_driving_forces([mol], env2)
        assert forces2['thermal'] > 0.0
    
    def test_electrochemical_driving_force(self):
        """Test electrochemical driving force evaluation."""
        mol = Molecule.from_smiles('CCO')
        
        # High voltage cathode
        env1 = Environment(electrode_type='cathode', voltage=4.5)
        forces1 = get_driving_forces([mol], env1)
        assert forces1['oxidation'] > 0.5
        assert forces1['electrochemical'] > 0.5
        
        # Low voltage anode
        env2 = Environment(electrode_type='anode', voltage=0.1)
        forces2 = get_driving_forces([mol], env2)
        assert forces2['reduction'] > 0.5


class TestOperatorRegistry:
    """Test operator registry system."""
    
    def test_registry_initialization(self):
        """Test operator registry initialization."""
        operators = list(OPERATOR_REGISTRY._operators.keys())
        expected_operators = [
            'hydrogen_transfer', 'bond_breaking', 'cyclization',
            'addition', 'rearrangement', 'redox'
        ]
        
        for op in expected_operators:
            assert op in operators
    
    def test_get_operator(self):
        """Test getting operator instances."""
        redox_op = OPERATOR_REGISTRY.get_operator('redox')
        assert redox_op is not None
        assert redox_op.name == 'Redox'
        
        # Test caching
        redox_op2 = OPERATOR_REGISTRY.get_operator('redox')
        assert redox_op is redox_op2  # Same instance
    
    def test_active_operators(self):
        """Test active operator selection."""
        mol = Molecule.from_smiles('CCO')
        env = Environment(temperature=500.0, electrode_type='cathode', voltage=4.2)
        
        active_ops = OPERATOR_REGISTRY.get_active_operators([mol], env)
        assert len(active_ops) > 0
        
        op_names = [op.name for op in active_ops]
        assert 'Redox' in op_names  # Should be active in oxidizing environment
    
    def test_recommended_operators(self):
        """Test recommended operator selection."""
        mol = Molecule.from_smiles('CCO')
        env = Environment(temperature=500.0, electrode_type='cathode', voltage=4.2)
        
        recommended = OPERATOR_REGISTRY.get_recommended_operators([mol], env, max_operators=3)
        assert len(recommended) <= 3
        assert len(recommended) > 0


class TestNewOperators:
    """Test new operator implementations."""
    
    def test_redox_operator(self):
        """Test redox operator functionality."""
        redox_op = RedoxOperator()
        
        # Test can_apply
        mol = Molecule.from_smiles('CCO')
        env_oxidizing = Environment(electrode_type='cathode', voltage=4.2)
        env_neutral = Environment(temperature=298.15)
        
        assert redox_op.can_apply([mol], env_oxidizing) == True
        assert redox_op.can_apply([mol], env_neutral) == False
        
        # Test reaction site finding
        sites = redox_op.find_reaction_sites([mol])
        assert isinstance(sites, list)
    
    def test_addition_operator(self):
        """Test addition operator functionality."""
        add_op = AdditionOperator()

        # Test with alkene in solution phase
        mol = Molecule.from_smiles('C=C')
        env = Environment(temperature=298.15, solvent='water')  # Solution phase

        assert add_op.can_apply([mol], env) == True
        
        sites = add_op.find_reaction_sites([mol])
        assert isinstance(sites, list)
    
    def test_rearrangement_operator(self):
        """Test rearrangement operator functionality."""
        rearr_op = RearrangementOperator()
        
        mol = Molecule.from_smiles('CCCC')
        env = Environment(temperature=500.0)
        
        assert rearr_op.can_apply([mol], env) == True
        
        sites = rearr_op.find_reaction_sites([mol])
        assert isinstance(sites, list)


class TestIntegration:
    """Integration tests for enhanced features."""
    
    def test_full_workflow(self):
        """Test complete workflow with enhanced features."""
        # Create molecule and environment
        mol = Molecule.from_smiles('CCO')
        env = Environment(
            temperature=500.0,
            electrode_type='cathode',
            voltage=4.2,
            li_activity=0.5
        )
        
        # Get structure tags
        tags = get_structure_tags([mol])
        assert len(tags) > 0
        
        # Get driving forces
        forces = get_driving_forces([mol], env)
        assert len(forces) > 0
        assert forces['oxidation'] > 0.0
        
        # Get active operators
        active_ops = OPERATOR_REGISTRY.get_active_operators([mol], env)
        assert len(active_ops) > 0
        
        # Test operator application
        for op in active_ops[:2]:  # Test first 2 operators
            sites = op.find_reaction_sites([mol])
            assert isinstance(sites, list)
    
    def test_different_molecules(self):
        """Test with different types of molecules."""
        molecules = [
            ('CCO', 'alcohol'),
            ('C=C', 'alkene'),
            ('C=O', 'carbonyl'),
            ('c1ccccc1', 'aromatic')
        ]
        
        env = Environment(temperature=400.0)
        
        for smiles, mol_type in molecules:
            mol = Molecule.from_smiles(smiles)
            
            # Should be able to get structure tags
            tags = get_structure_tags([mol])
            assert len(tags) > 0
            
            # Should be able to get driving forces
            forces = get_driving_forces([mol], env)
            assert len(forces) > 0
            
            # Should find some active operators
            active_ops = OPERATOR_REGISTRY.get_active_operators([mol], env)
            # Note: Not all molecules will have active operators in all environments
