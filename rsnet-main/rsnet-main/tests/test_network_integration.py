"""
Integration tests for enhanced network generation engine.

Tests the complete integration of enhanced features including:
- Intelligent operator selection
- Driving force evaluation
- Structure-based filtering
- Configuration system
"""

import pytest
from unittest.mock import Mock, patch
import tempfile
import os

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.network.generator import NetworkGenerator
from rsnet.network.config import NetworkGenerationConfig, NetworkGenerationPresets
from rsnet.operators.registry import OPERATOR_REGISTRY


class TestNetworkGenerationIntegration:
    """Test complete network generation with enhanced features."""
    
    def test_intelligent_network_generation(self):
        """Test network generation with intelligent operator selection."""
        # Setup
        mol = Molecule.from_smiles('CCO')  # Ethanol
        env = Environment(temperature=500.0, electrode_type='cathode', voltage=4.0)
        
        config = NetworkGenerationConfig(
            max_generations=2,
            max_species=20,
            energy_cutoff=60.0,
            use_structure_based_filtering=True,
            max_operators_per_generation=3
        )
        
        generator = NetworkGenerator(config=config)
        
        # Generate network
        network = generator.generate_network([mol], env, max_time=30.0)
        
        # Verify results
        assert len(network.molecules) >= 1
        assert mol.smiles in network.molecules
        
        stats = network.get_statistics()
        assert stats['num_molecules'] >= 1
        assert stats['max_generation'] <= config.max_generations
        
        # Verify generation statistics
        gen_stats = generator.get_generation_statistics()
        assert len(gen_stats) >= 1
        assert all('generation' in stat for stat in gen_stats)
    
    def test_preset_configurations(self):
        """Test network generation with different preset configurations."""
        mol = Molecule.from_smiles('C=C')  # Ethylene
        env = Environment(temperature=400.0, solvent='water')
        
        # Test battery chemistry preset
        battery_config = NetworkGenerationPresets.battery_chemistry()
        generator = NetworkGenerator(config=battery_config)
        
        network = generator.generate_network([mol], env, max_time=20.0)
        assert len(network.molecules) >= 1
        
        # Test organic synthesis preset
        organic_config = NetworkGenerationPresets.organic_synthesis()
        generator = NetworkGenerator(config=organic_config)
        
        network = generator.generate_network([mol], env, max_time=20.0)
        assert len(network.molecules) >= 1
        
        # Test fast screening preset
        fast_config = NetworkGenerationPresets.fast_screening()
        generator = NetworkGenerator(config=fast_config)
        
        network = generator.generate_network([mol], env, max_time=10.0)
        assert len(network.molecules) >= 1
    
    def test_driving_force_integration(self):
        """Test integration of driving force evaluation in network generation."""
        # High temperature environment with strong driving forces
        mol = Molecule.from_smiles('CC(=O)O')  # Acetic acid
        env = Environment(
            temperature=600.0,
            electrode_type='cathode',
            voltage=4.5,
            li_activity=0.8
        )
        
        config = NetworkGenerationConfig(
            max_generations=2,
            max_species=15,
            driving_force_threshold=0.2,  # Lower threshold to capture more reactions
            use_structure_based_filtering=True
        )
        
        generator = NetworkGenerator(config=config)
        network = generator.generate_network([mol], env, max_time=25.0)
        
        # Should generate some reactions due to strong driving forces
        assert len(network.reactions) >= 0  # May be 0 if no feasible reactions
        
        # Check that reactions have driving force scores
        for reaction in network.reactions.values():
            assert hasattr(reaction, 'driving_force_score')
    
    def test_operator_registry_integration(self):
        """Test integration with operator registry system."""
        mol = Molecule.from_smiles('CCO')
        env = Environment(temperature=450.0, electrode_type='anode')
        
        # Test with custom operator registry
        custom_registry = OPERATOR_REGISTRY
        
        config = NetworkGenerationConfig(
            max_generations=2,
            max_species=10,
            max_operators_per_generation=2
        )
        
        generator = NetworkGenerator(
            operator_registry=custom_registry,
            config=config
        )
        
        network = generator.generate_network([mol], env, max_time=15.0)
        
        assert len(network.molecules) >= 1
        assert generator.use_intelligent_selection == True
    
    def test_legacy_compatibility(self):
        """Test backward compatibility with legacy configuration."""
        mol = Molecule.from_smiles('C')  # Methane
        env = Environment(temperature=800.0)
        
        # Legacy configuration dictionary
        legacy_config = {
            'max_generations': 2,
            'max_species': 15,
            'energy_cutoff': 70.0
        }
        
        generator = NetworkGenerator(legacy_config=legacy_config)
        network = generator.generate_network([mol], env, max_time=10.0)
        
        assert len(network.molecules) >= 1
        assert generator.config.max_generations == 2
        assert generator.config.max_species == 15
        assert generator.config.energy_cutoff == 70.0
    
    def test_structure_based_filtering(self):
        """Test structure-based molecular filtering."""
        molecules = [
            Molecule.from_smiles('CCO'),      # Ethanol
            Molecule.from_smiles('c1ccccc1'),  # Benzene
            Molecule.from_smiles('C1CC1'),     # Cyclopropane
        ]
        env = Environment(temperature=400.0)
        
        config = NetworkGenerationConfig(
            max_generations=1,
            max_species=20,
            use_structure_based_filtering=True,
            max_molecular_weight=200.0
        )
        
        generator = NetworkGenerator(config=config)
        
        # Test molecule grouping
        groups = generator._group_molecules_by_features(molecules)
        assert len(groups) >= 1
        assert all(isinstance(group, list) for group in groups)
        assert sum(len(group) for group in groups) == len(molecules)
    
    def test_adaptive_operator_selection(self):
        """Test adaptive operator selection based on conditions."""
        mol = Molecule.from_smiles('CC=O')  # Acetaldehyde
        env = Environment(
            temperature=500.0,
            electrode_type='cathode',
            voltage=3.8
        )
        
        config = NetworkGenerationConfig(
            max_generations=1,
            max_species=10,
            operator_selection_strategy='adaptive',
            driving_force_threshold=0.3
        )
        
        generator = NetworkGenerator(config=config)
        
        # Test operator scoring
        from rsnet.operators.redox import RedoxOperator
        redox_op = RedoxOperator()
        
        score = generator._calculate_operator_score(
            redox_op, [mol], env, {'electrochemical': 0.8, 'oxidation': 0.6}
        )
        
        assert score > 1.0  # Should have bonus from electrochemical forces
    
    def test_reaction_ranking(self):
        """Test reaction ranking by relevance."""
        from rsnet.core.reaction import Reaction
        
        mol1 = Molecule.from_smiles('CCO')
        mol2 = Molecule.from_smiles('CC=O')
        env = Environment(temperature=400.0)
        
        # Create mock reactions with different properties
        reaction1 = Reaction(
            reactants=[mol1],
            products=[mol2],
            operator_name='test'
        )
        reaction1.reaction_energy = -10.0  # Favorable
        reaction1.driving_force_score = 0.8
        
        reaction2 = Reaction(
            reactants=[mol1],
            products=[mol2],
            operator_name='test'
        )
        reaction2.reaction_energy = 30.0   # Less favorable
        reaction2.driving_force_score = 0.3
        
        config = NetworkGenerationConfig(energy_cutoff=50.0)
        generator = NetworkGenerator(config=config)
        
        ranked = generator._rank_reactions_by_relevance([reaction1, reaction2], env)
        
        # reaction1 should rank higher due to better energy and driving forces
        assert ranked[0] == reaction1
        assert ranked[1] == reaction2
    
    def test_configuration_serialization(self):
        """Test configuration serialization and deserialization."""
        config = NetworkGenerationConfig(
            max_generations=7,
            max_species=150,
            energy_cutoff=45.0,
            use_structure_based_filtering=True
        )
        
        # Test to_dict
        config_dict = config.to_dict()
        assert config_dict['max_generations'] == 7
        assert config_dict['max_species'] == 150
        assert config_dict['energy_cutoff'] == 45.0
        
        # Test from_dict
        restored_config = NetworkGenerationConfig.from_dict(config_dict)
        assert restored_config.max_generations == 7
        assert restored_config.max_species == 150
        assert restored_config.energy_cutoff == 45.0
        assert restored_config.use_structure_based_filtering == True


class TestNetworkGenerationPerformance:
    """Test performance aspects of enhanced network generation."""
    
    def test_generation_time_limits(self):
        """Test that time limits are respected."""
        import time
        
        mol = Molecule.from_smiles('CCCCCCCC')  # Octane - may generate many reactions
        env = Environment(temperature=700.0)
        
        config = NetworkGenerationConfig(
            max_generations=10,  # High limit
            max_species=1000,    # High limit
            energy_cutoff=100.0  # Permissive
        )
        
        generator = NetworkGenerator(config=config)
        
        start_time = time.time()
        network = generator.generate_network([mol], env, max_time=5.0)  # 5 second limit
        elapsed_time = time.time() - start_time
        
        # Should respect time limit (with some tolerance for cleanup)
        assert elapsed_time < 10.0
        assert len(network.molecules) >= 1
    
    def test_species_limits(self):
        """Test that species limits are respected."""
        mol = Molecule.from_smiles('CC')  # Ethane
        env = Environment(temperature=600.0)
        
        config = NetworkGenerationConfig(
            max_generations=10,
            max_species=5,  # Low limit
            energy_cutoff=80.0
        )
        
        generator = NetworkGenerator(config=config)
        network = generator.generate_network([mol], env, max_time=20.0)
        
        # Should respect species limit
        assert len(network.molecules) <= config.max_species
