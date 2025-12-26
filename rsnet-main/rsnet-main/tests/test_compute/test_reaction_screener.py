"""
Tests for reaction screener.
"""

import pytest
from unittest.mock import Mock, patch

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.core.reaction import Reaction
from rsnet.compute.xtb_calculator import XTBCalculator
from rsnet.compute.reaction_screener import ReactionScreener


class TestReactionScreener:
    """Test cases for ReactionScreener."""
    
    def test_screener_initialization(self):
        """Test screener initialization."""
        screener = ReactionScreener()
        assert screener.calculator is not None
        assert screener.max_workers == 4
        assert screener.energy_threshold == 50.0
        assert screener.optimize_geometries is True
    
    def test_screener_with_custom_settings(self):
        """Test screener with custom settings."""
        calc = XTBCalculator()
        screener = ReactionScreener(
            calculator=calc,
            max_workers=2,
            energy_threshold=30.0,
            optimize_geometries=False
        )
        assert screener.calculator is calc
        assert screener.max_workers == 2
        assert screener.energy_threshold == 30.0
        assert screener.optimize_geometries is False
    
    def test_molecule_key_generation(self):
        """Test molecule cache key generation."""
        screener = ReactionScreener()
        mol = Molecule.from_smiles("CC", name="ethane")
        
        key = screener._get_molecule_key(mol)
        assert isinstance(key, str)
        assert "CC" in key
        assert "C2H6" in key
    
    def test_energy_caching(self):
        """Test energy calculation caching."""
        screener = ReactionScreener()
        mol = Molecule.from_smiles("C", name="methane")
        env = Environment(temperature=298.15)
        
        # Mock the calculator to return a fixed energy
        mock_result = {
            'success': True,
            'energy_kcal_mol': -100.0
        }
        
        with patch.object(screener.calculator, 'single_point', return_value=mock_result):
            # First call should hit the calculator
            energy1 = screener._calculate_molecule_energy(mol, env)
            assert energy1 == -100.0
            
            # Second call should use cache
            energy2 = screener._calculate_molecule_energy(mol, env)
            assert energy2 == -100.0
            
            # Verify cache was used
            assert len(screener._energy_cache) == 1
    
    def test_calculate_reaction_energy_simple(self):
        """Test reaction energy calculation for simple reaction."""
        screener = ReactionScreener()
        
        # Create simple reaction: A -> B
        reactant = Molecule.from_smiles("C", name="methane")
        product = Molecule.from_smiles("CC", name="ethane")
        reaction = Reaction([reactant], [product], name="test_reaction")
        env = Environment(temperature=298.15)
        
        # Mock energies
        def mock_energy(mol, environment):
            if mol.name == "methane":
                return -50.0
            elif mol.name == "ethane":
                return -80.0
            return None
        
        with patch.object(screener, '_calculate_molecule_energy', side_effect=mock_energy):
            result = screener.calculate_reaction_energy(reaction, env)
            
            assert result['success'] is True
            assert result['reaction_energy'] == -30.0  # -80 - (-50)
            assert result['exothermic'] is True
            assert result['endothermic'] is False
            assert result['total_reactant_energy'] == -50.0
            assert result['total_product_energy'] == -80.0
    
    def test_calculate_reaction_energy_multiple_molecules(self):
        """Test reaction energy calculation with multiple reactants/products."""
        screener = ReactionScreener()
        
        # Create reaction: A + B -> C + D
        reactant1 = Molecule.from_smiles("C", name="methane")
        reactant2 = Molecule.from_smiles("O", name="water")
        product1 = Molecule.from_smiles("CO", name="methanol")
        product2 = Molecule.from_smiles("[H][H]", name="hydrogen")
        
        reaction = Reaction([reactant1, reactant2], [product1, product2], name="test_reaction")
        env = Environment(temperature=298.15)
        
        # Mock energies
        energy_map = {
            "methane": -40.0,
            "water": -60.0,
            "methanol": -90.0,
            "hydrogen": -10.0
        }
        
        def mock_energy(mol, environment):
            return energy_map.get(mol.name)
        
        with patch.object(screener, '_calculate_molecule_energy', side_effect=mock_energy):
            result = screener.calculate_reaction_energy(reaction, env)
            
            assert result['success'] is True
            # Reaction energy = (-90 + -10) - (-40 + -60) = -100 - (-100) = 0
            assert result['reaction_energy'] == 0.0
            assert result['exothermic'] is False
            assert result['endothermic'] is False
    
    def test_calculate_reaction_energy_failure(self):
        """Test reaction energy calculation when molecule energy fails."""
        screener = ReactionScreener()
        
        reactant = Molecule.from_smiles("C", name="methane")
        product = Molecule.from_smiles("CC", name="ethane")
        reaction = Reaction([reactant], [product], name="test_reaction")
        env = Environment(temperature=298.15)
        
        # Mock energy calculation to fail for product
        def mock_energy(mol, environment):
            if mol.name == "methane":
                return -50.0
            return None  # Fail for ethane
        
        with patch.object(screener, '_calculate_molecule_energy', side_effect=mock_energy):
            result = screener.calculate_reaction_energy(reaction, env)
            
            assert result['success'] is False
            assert 'error' in result
            assert 'ethane' in result['error']
    
    def test_estimate_activation_energy_endothermic(self):
        """Test activation energy estimation for endothermic reaction."""
        screener = ReactionScreener()
        
        reactant = Molecule.from_smiles("C", name="methane")
        product = Molecule.from_smiles("CC", name="ethane")
        reaction = Reaction([reactant], [product], name="test_reaction")
        env = Environment(temperature=298.15)
        
        # Mock endothermic reaction (ΔE = +20 kcal/mol)
        mock_energetics = {
            'success': True,
            'reaction_energy': 20.0,
            'exothermic': False,
            'endothermic': True
        }
        
        with patch.object(screener, 'calculate_reaction_energy', return_value=mock_energetics):
            ea = screener.estimate_activation_energy(reaction, env)
            
            assert ea is not None
            assert ea > 20.0  # Should be higher than reaction energy
            # Expected: 0.5 * 20 + 15 = 25 kcal/mol
            assert abs(ea - 25.0) < 1.0
    
    def test_estimate_activation_energy_exothermic(self):
        """Test activation energy estimation for exothermic reaction."""
        screener = ReactionScreener()
        
        reactant = Molecule.from_smiles("CC", name="ethane")
        product = Molecule.from_smiles("C", name="methane")
        reaction = Reaction([reactant], [product], name="test_reaction")
        env = Environment(temperature=298.15)
        
        # Mock exothermic reaction (ΔE = -20 kcal/mol)
        mock_energetics = {
            'success': True,
            'reaction_energy': -20.0,
            'exothermic': True,
            'endothermic': False
        }
        
        with patch.object(screener, 'calculate_reaction_energy', return_value=mock_energetics):
            ea = screener.estimate_activation_energy(reaction, env)
            
            assert ea is not None
            # Expected: 10 kcal/mol (intrinsic barrier for exothermic)
            assert abs(ea - 10.0) < 1.0
    
    def test_screen_reaction_feasible(self):
        """Test screening a feasible reaction."""
        screener = ReactionScreener(energy_threshold=50.0)
        
        reactant = Molecule.from_smiles("C", name="methane")
        product = Molecule.from_smiles("CC", name="ethane")
        reaction = Reaction([reactant], [product], name="test_reaction")
        env = Environment(temperature=298.15)
        
        # Mock feasible reaction
        mock_energetics = {
            'success': True,
            'reaction_energy': 10.0,  # Within threshold
            'exothermic': False,
            'endothermic': True,
            'reactant_energies': [-50.0],
            'product_energies': [-40.0]
        }
        
        with patch.object(screener, 'calculate_reaction_energy', return_value=mock_energetics):
            with patch.object(screener, 'estimate_activation_energy', return_value=25.0):
                result = screener.screen_reaction(reaction, env)
                
                assert result['success'] is True
                assert result['feasible'] is True
                assert result['thermodynamically_feasible'] is True
                assert result['kinetically_accessible'] is True
                assert result['reaction_energy'] == 10.0
                assert result['activation_energy'] == 25.0
                assert result['rate_constant'] is not None
    
    def test_screen_reaction_unfeasible_thermodynamics(self):
        """Test screening a thermodynamically unfeasible reaction."""
        screener = ReactionScreener(energy_threshold=30.0)
        
        reactant = Molecule.from_smiles("C", name="methane")
        product = Molecule.from_smiles("CC", name="ethane")
        reaction = Reaction([reactant], [product], name="test_reaction")
        env = Environment(temperature=298.15)
        
        # Mock unfeasible reaction (high energy)
        mock_energetics = {
            'success': True,
            'reaction_energy': 50.0,  # Above threshold
            'exothermic': False,
            'endothermic': True,
            'reactant_energies': [-50.0],
            'product_energies': [0.0]
        }
        
        with patch.object(screener, 'calculate_reaction_energy', return_value=mock_energetics):
            with patch.object(screener, 'estimate_activation_energy', return_value=40.0):
                result = screener.screen_reaction(reaction, env)
                
                assert result['success'] is True
                assert result['feasible'] is False
                assert result['thermodynamically_feasible'] is False
                assert result['kinetically_accessible'] is False
    
    def test_screen_reactions_multiple(self):
        """Test screening multiple reactions."""
        screener = ReactionScreener(max_workers=1)  # Sequential for predictable testing
        
        # Create multiple reactions
        reactions = []
        for i in range(3):
            reactant = Molecule.from_smiles("C", name=f"reactant_{i}")
            product = Molecule.from_smiles("CC", name=f"product_{i}")
            reaction = Reaction([reactant], [product], name=f"reaction_{i}")
            reactions.append(reaction)
        
        env = Environment(temperature=298.15)
        
        # Mock screening results
        def mock_screen(reaction, environment):
            return {
                'reaction': reaction,
                'success': True,
                'feasible': True,
                'reaction_energy': 10.0,
                'calculation_time': 1.0
            }
        
        with patch.object(screener, 'screen_reaction', side_effect=mock_screen):
            results = screener.screen_reactions(reactions, env)
            
            assert len(results) == 3
            for result in results:
                assert result['success'] is True
                assert result['feasible'] is True
    
    def test_filter_feasible_reactions(self):
        """Test filtering feasible reactions."""
        screener = ReactionScreener()
        
        # Create test reactions
        reactions = []
        for i in range(3):
            reactant = Molecule.from_smiles("C", name=f"reactant_{i}")
            product = Molecule.from_smiles("CC", name=f"product_{i}")
            reaction = Reaction([reactant], [product], name=f"reaction_{i}")
            reactions.append(reaction)
        
        env = Environment(temperature=298.15)
        
        # Mock screening results - only first two are feasible
        mock_results = [
            {'reaction': reactions[0], 'success': True, 'feasible': True},
            {'reaction': reactions[1], 'success': True, 'feasible': True},
            {'reaction': reactions[2], 'success': True, 'feasible': False}
        ]
        
        with patch.object(screener, 'screen_reactions', return_value=mock_results):
            feasible = screener.filter_feasible_reactions(reactions, env)
            
            assert len(feasible) == 2
            assert feasible[0] is reactions[0]
            assert feasible[1] is reactions[1]
    
    def test_get_energy_statistics_empty(self):
        """Test energy statistics with empty cache."""
        screener = ReactionScreener()
        stats = screener.get_energy_statistics()
        
        assert stats['cache_size'] == 0
    
    def test_get_energy_statistics_with_data(self):
        """Test energy statistics with cached data."""
        screener = ReactionScreener()
        
        # Add some mock data to cache
        screener._energy_cache = {
            'mol1': -50.0,
            'mol2': -100.0,
            'mol3': -75.0
        }
        
        stats = screener.get_energy_statistics()
        
        assert stats['cache_size'] == 3
        assert stats['min_energy'] == -100.0
        assert stats['max_energy'] == -50.0
        assert stats['mean_energy'] == -75.0


@pytest.mark.slow
class TestReactionScreenerIntegration:
    """Integration tests for ReactionScreener with real xTB calculations."""
    
    def test_screen_simple_reaction(self):
        """Test screening a simple reaction with real calculations."""
        screener = ReactionScreener(optimize_geometries=False)  # Faster
        
        # Simple reaction: H2 -> 2H
        reactant = Molecule.from_smiles("[H][H]", name="hydrogen")
        product1 = Molecule.from_smiles("[H]", name="hydrogen_atom_1")
        product2 = Molecule.from_smiles("[H]", name="hydrogen_atom_2")
        
        reaction = Reaction([reactant], [product1, product2], name="H2_dissociation")
        env = Environment(temperature=298.15)
        
        result = screener.screen_reaction(reaction, env)
        
        assert result['success'] is True
        assert result['reaction_energy'] is not None
        assert result['activation_energy'] is not None
        
        # H2 dissociation should be highly endothermic
        assert result['reaction_energy'] > 50.0
        assert result['endothermic'] is True
        
        print(f"H2 dissociation:")
        print(f"  ΔE = {result['reaction_energy']:.2f} kcal/mol")
        print(f"  Ea = {result['activation_energy']:.2f} kcal/mol")
        print(f"  Feasible: {result['feasible']}")
    
    def test_energy_caching_real(self):
        """Test energy caching with real calculations."""
        screener = ReactionScreener(optimize_geometries=False)
        
        mol = Molecule.from_smiles("C", name="methane")
        env = Environment(temperature=298.15)
        
        # First calculation
        energy1 = screener._calculate_molecule_energy(mol, env)
        assert energy1 is not None
        
        # Second calculation should use cache
        energy2 = screener._calculate_molecule_energy(mol, env)
        assert energy2 == energy1
        
        # Check cache
        assert len(screener._energy_cache) == 1
        
        print(f"Methane energy: {energy1:.2f} kcal/mol (cached)")
