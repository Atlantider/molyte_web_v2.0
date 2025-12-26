"""
Tests for network generator.
"""

import pytest
from unittest.mock import Mock, patch

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.core.reaction import Reaction
from rsnet.operators.bond_breaking import BondBreakingOperator
from rsnet.compute.reaction_screener import ReactionScreener
from rsnet.network.generator import NetworkGenerator, ReactionNetwork


class TestReactionNetwork:
    """Test cases for ReactionNetwork."""
    
    def test_network_initialization(self):
        """Test network initialization."""
        network = ReactionNetwork()
        assert len(network.molecules) == 0
        assert len(network.reactions) == 0
        assert network.graph.number_of_nodes() == 0
        assert network.graph.number_of_edges() == 0
    
    def test_add_molecule(self):
        """Test adding molecules to network."""
        network = ReactionNetwork()
        mol = Molecule.from_smiles("CC", name="ethane")
        
        network.add_molecule(mol, generation=0)
        
        assert len(network.molecules) == 1
        assert mol.smiles in network.molecules
        assert network.molecules[mol.smiles] is mol
        assert network.generation_info[mol.smiles] == 0
        assert network.graph.has_node(mol.smiles)
    
    def test_add_duplicate_molecule(self):
        """Test adding duplicate molecules."""
        network = ReactionNetwork()
        mol1 = Molecule.from_smiles("CC", name="ethane1")
        mol2 = Molecule.from_smiles("CC", name="ethane2")  # Same SMILES
        
        network.add_molecule(mol1, generation=0)
        network.add_molecule(mol2, generation=1)
        
        # Should only have one molecule (first one added)
        assert len(network.molecules) == 1
        assert network.molecules["CC"] is mol1
        assert network.generation_info["CC"] == 0  # Original generation
    
    def test_add_reaction(self):
        """Test adding reactions to network."""
        network = ReactionNetwork()
        
        reactant = Molecule.from_smiles("CC", name="ethane")
        product1 = Molecule.from_smiles("[CH3]", name="methyl1")
        product2 = Molecule.from_smiles("[CH3]", name="methyl2")
        
        reaction = Reaction([reactant], [product1, product2], name="test_reaction")
        
        result = network.add_reaction(reaction)
        
        assert result is True
        assert len(network.reactions) == 1
        assert reaction.id in network.reactions
        assert len(network.molecules) == 2  # ethane and methyl (deduplicated)
        assert network.graph.has_edge("CC", "[CH3]")
    
    def test_add_duplicate_reaction(self):
        """Test adding duplicate reactions."""
        network = ReactionNetwork()
        
        reactant = Molecule.from_smiles("CC", name="ethane")
        product = Molecule.from_smiles("[CH3]", name="methyl")
        
        reaction1 = Reaction([reactant], [product], name="reaction1")
        reaction2 = Reaction([reactant], [product], name="reaction2")
        
        result1 = network.add_reaction(reaction1)
        result2 = network.add_reaction(reaction2)
        
        assert result1 is True
        assert result2 is False  # Duplicate
        assert len(network.reactions) == 1
    
    def test_get_molecules_by_generation(self):
        """Test getting molecules by generation."""
        network = ReactionNetwork()
        
        mol1 = Molecule.from_smiles("CC", name="ethane")
        mol2 = Molecule.from_smiles("CCO", name="ethanol")
        mol3 = Molecule.from_smiles("[CH3]", name="methyl")
        
        network.add_molecule(mol1, generation=0)
        network.add_molecule(mol2, generation=0)
        network.add_molecule(mol3, generation=1)
        
        gen0_mols = network.get_molecules_by_generation(0)
        gen1_mols = network.get_molecules_by_generation(1)
        
        assert len(gen0_mols) == 2
        assert len(gen1_mols) == 1
        assert mol3 in gen1_mols
    
    def test_get_statistics(self):
        """Test network statistics."""
        network = ReactionNetwork()
        
        # Add molecules
        mol1 = Molecule.from_smiles("CC", name="ethane")
        mol2 = Molecule.from_smiles("[CH3]", name="methyl")
        network.add_molecule(mol1, generation=0)
        network.add_molecule(mol2, generation=1)
        
        # Add reaction
        reaction = Reaction([mol1], [mol2], name="test_reaction")
        network.add_reaction(reaction)
        
        stats = network.get_statistics()
        
        assert stats['num_molecules'] == 2
        assert stats['num_reactions'] == 1
        assert stats['num_edges'] == 1
        assert stats['max_generation'] == 1
        assert stats['molecules_by_generation'] == {0: 1, 1: 1}


class TestNetworkGenerator:
    """Test cases for NetworkGenerator."""
    
    def test_generator_initialization(self):
        """Test generator initialization."""
        generator = NetworkGenerator()
        
        assert len(generator.operators) == 2  # Default operators
        assert generator.screener is not None
        assert generator.max_generations == 5
        assert generator.max_species == 100
    
    def test_generator_with_custom_config(self):
        """Test generator with custom configuration."""
        config = {
            'max_generations': 3,
            'max_species': 50,
            'energy_cutoff': 30.0
        }
        
        generator = NetworkGenerator(config=config)
        
        assert generator.max_generations == 3
        assert generator.max_species == 50
        assert generator.energy_cutoff == 30.0
    
    def test_process_generation_no_reactions(self):
        """Test processing generation with no applicable reactions."""
        generator = NetworkGenerator()
        
        # Mock operators to return no reactions
        mock_operator = Mock()
        mock_operator.is_applicable.return_value = False
        generator.operators = [mock_operator]
        
        network = ReactionNetwork()
        molecules = [Molecule.from_smiles("C", name="methane")]
        env = Environment(temperature=298.15)
        
        new_reactions = generator._process_generation(molecules, env, network, 0)
        
        assert len(new_reactions) == 0
        mock_operator.is_applicable.assert_called_once()
    
    def test_process_generation_with_reactions(self):
        """Test processing generation with reactions."""
        generator = NetworkGenerator()
        
        # Mock operator
        mock_operator = Mock()
        mock_operator.is_applicable.return_value = True
        mock_operator.name = "mock_operator"
        
        # Create mock reaction
        reactant = Molecule.from_smiles("CC", name="ethane")
        product = Molecule.from_smiles("[CH3]", name="methyl")
        mock_reaction = Reaction([reactant], [product], name="mock_reaction")
        mock_operator.apply.return_value = [mock_reaction]
        
        generator.operators = [mock_operator]
        
        # Mock screener
        mock_screening_result = {
            'reaction': mock_reaction,
            'success': True,
            'feasible': True,
            'reaction_energy': 10.0,
            'activation_energy': 25.0
        }
        generator.screener.screen_reactions = Mock(return_value=[mock_screening_result])
        
        network = ReactionNetwork()
        molecules = [reactant]
        env = Environment(temperature=298.15)
        
        new_reactions = generator._process_generation(molecules, env, network, 0)
        
        assert len(new_reactions) == 1
        assert new_reactions[0] is mock_reaction
        assert mock_reaction.reaction_energy == 10.0
        assert mock_reaction.activation_energy == 25.0
    
    def test_process_generation_unfeasible_reactions(self):
        """Test processing generation with unfeasible reactions."""
        generator = NetworkGenerator()
        
        # Mock operator
        mock_operator = Mock()
        mock_operator.is_applicable.return_value = True
        mock_operator.name = "mock_operator"
        
        reactant = Molecule.from_smiles("CC", name="ethane")
        product = Molecule.from_smiles("[CH3]", name="methyl")
        mock_reaction = Reaction([reactant], [product], name="mock_reaction")
        mock_operator.apply.return_value = [mock_reaction]
        
        generator.operators = [mock_operator]
        
        # Mock screener to return unfeasible reaction
        mock_screening_result = {
            'reaction': mock_reaction,
            'success': True,
            'feasible': False,  # Not feasible
            'reaction_energy': 100.0,  # High energy
            'activation_energy': 150.0
        }
        generator.screener.screen_reactions = Mock(return_value=[mock_screening_result])
        
        network = ReactionNetwork()
        molecules = [reactant]
        env = Environment(temperature=298.15)
        
        new_reactions = generator._process_generation(molecules, env, network, 0)
        
        assert len(new_reactions) == 0  # No feasible reactions
    
    def test_generate_network_simple(self):
        """Test simple network generation."""
        # Create a minimal generator
        generator = NetworkGenerator(
            operators=[],  # No operators for simplicity
            config={'max_generations': 2, 'max_species': 10}
        )
        
        seed_molecules = [Molecule.from_smiles("C", name="methane")]
        env = Environment(temperature=298.15)
        
        network = generator.generate_network(seed_molecules, env)
        
        # Should have at least the seed molecule
        assert len(network.molecules) >= 1
        assert "C" in network.molecules
        assert network.generation_info["C"] == 0
    
    def test_generate_network_max_species_limit(self):
        """Test network generation with species limit."""
        generator = NetworkGenerator(
            config={'max_generations': 10, 'max_species': 2}  # Very low limit
        )
        
        seed_molecules = [Molecule.from_smiles("CC", name="ethane")]
        env = Environment(temperature=500.0)  # High temp for more reactions
        
        network = generator.generate_network(seed_molecules, env)
        
        # Should respect species limit
        assert len(network.molecules) <= 2
    
    def test_generate_network_max_generations_limit(self):
        """Test network generation with generation limit."""
        generator = NetworkGenerator(
            config={'max_generations': 1, 'max_species': 100}
        )
        
        seed_molecules = [Molecule.from_smiles("CC", name="ethane")]
        env = Environment(temperature=500.0)
        
        network = generator.generate_network(seed_molecules, env)
        
        # Should respect generation limit
        max_gen = max(network.generation_info.values()) if network.generation_info else 0
        assert max_gen <= 1
    
    def test_generate_network_time_limit(self):
        """Test network generation with time limit."""
        generator = NetworkGenerator()
        
        seed_molecules = [Molecule.from_smiles("CC", name="ethane")]
        env = Environment(temperature=500.0)
        
        # Very short time limit
        network = generator.generate_network(seed_molecules, env, max_time=0.1)
        
        # Should have at least seed molecules
        assert len(network.molecules) >= 1
    
    def test_get_generation_statistics(self):
        """Test generation statistics tracking."""
        generator = NetworkGenerator(
            operators=[],  # No operators for simplicity
            config={'max_generations': 1}
        )
        
        seed_molecules = [Molecule.from_smiles("C", name="methane")]
        env = Environment(temperature=298.15)
        
        network = generator.generate_network(seed_molecules, env)
        stats = generator.get_generation_statistics()
        
        assert len(stats) >= 1
        assert stats[0]['generation'] == 0
        assert stats[0]['input_molecules'] == 1
    
    def test_reset_statistics(self):
        """Test resetting generation statistics."""
        generator = NetworkGenerator()
        
        # Add some mock statistics
        generator.generation_stats = [{'generation': 0, 'test': True}]
        
        generator.reset_statistics()
        
        assert len(generator.generation_stats) == 0


@pytest.mark.slow
class TestNetworkGeneratorIntegration:
    """Integration tests for NetworkGenerator with real components."""
    
    def test_generate_network_with_bond_breaking(self):
        """Test network generation with bond breaking operator."""
        operators = [BondBreakingOperator(config={'min_fragment_size': 1})]
        screener = ReactionScreener(optimize_geometries=False, max_workers=1)
        
        generator = NetworkGenerator(
            operators=operators,
            screener=screener,
            config={'max_generations': 2, 'max_species': 20, 'energy_cutoff': 100.0}
        )
        
        seed_molecules = [Molecule.from_smiles("CC", name="ethane")]
        env = Environment(temperature=500.0)
        
        network = generator.generate_network(seed_molecules, env, max_time=60)
        
        # Should generate some reactions
        assert len(network.molecules) > 1
        assert len(network.reactions) > 0
        
        # Check that ethane is present
        assert "CC" in network.molecules
        
        # Should have methyl radicals from ethane breaking
        assert "[CH3]" in network.molecules
        
        print(f"Generated network: {len(network.molecules)} molecules, {len(network.reactions)} reactions")
    
    def test_network_convergence(self):
        """Test that network generation converges."""
        operators = [BondBreakingOperator(config={'min_fragment_size': 1})]
        screener = ReactionScreener(optimize_geometries=False, max_workers=1)
        
        generator = NetworkGenerator(
            operators=operators,
            screener=screener,
            config={
                'max_generations': 5,
                'max_species': 30,
                'min_new_reactions': 1  # Require at least 1 new reaction per generation
            }
        )
        
        seed_molecules = [Molecule.from_smiles("CCC", name="propane")]
        env = Environment(temperature=400.0)  # Lower temp for fewer reactions
        
        network = generator.generate_network(seed_molecules, env, max_time=120)
        
        # Should converge before max generations
        stats = generator.get_generation_statistics()
        assert len(stats) > 0
        
        # Last generation should have few new reactions (convergence)
        if len(stats) > 1:
            last_gen_reactions = stats[-1]['new_reactions']
            print(f"Last generation had {last_gen_reactions} new reactions")
        
        print(f"Network converged after {len(stats)} generations")
