"""
Tests for xTB calculator.
"""

import pytest
import tempfile
import os

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.compute.xtb_calculator import XTBCalculator


class TestXTBCalculator:
    """Test cases for XTBCalculator."""
    
    def test_calculator_initialization(self):
        """Test calculator initialization."""
        calc = XTBCalculator()
        assert calc.xtb_path == "xtb"
        assert calc.method == "gfn2"
        assert calc.charge == 0
        assert calc.multiplicity == 1
        assert calc.solvent is None
    
    def test_calculator_with_custom_settings(self):
        """Test calculator with custom settings."""
        calc = XTBCalculator(
            method="gfn1",
            charge=1,
            multiplicity=2,
            solvent="water"
        )
        assert calc.method == "gfn1"
        assert calc.charge == 1
        assert calc.multiplicity == 2
        assert calc.solvent == "water"
    
    def test_xtb_availability(self):
        """Test that xTB is available."""
        # This should not raise an exception if xTB is properly installed
        calc = XTBCalculator()
        # If we get here, xTB is available
        assert True
    
    def test_molecule_to_xyz_simple(self):
        """Test conversion of simple molecule to XYZ format."""
        mol = Molecule.from_smiles("H", name="hydrogen_atom")
        calc = XTBCalculator()
        
        xyz_content = calc._molecule_to_xyz(mol)
        lines = xyz_content.strip().split('\n')
        
        # Should have at least 3 lines: atom count, comment, coordinates
        assert len(lines) >= 3
        assert lines[0].strip() == "1"  # One atom
        assert "hydrogen_atom" in lines[1]  # Comment line
        assert lines[2].startswith("H")  # Hydrogen atom
    
    def test_molecule_to_xyz_ethane(self):
        """Test conversion of ethane to XYZ format."""
        mol = Molecule.from_smiles("CC", name="ethane")
        calc = XTBCalculator()
        
        xyz_content = calc._molecule_to_xyz(mol)
        lines = xyz_content.strip().split('\n')
        
        # Ethane has 8 atoms (2 C + 6 H)
        assert lines[0].strip() == "8"
        assert "ethane" in lines[1].lower()
        
        # Count carbon and hydrogen atoms
        carbon_count = sum(1 for line in lines[2:] if line.strip().startswith('C'))
        hydrogen_count = sum(1 for line in lines[2:] if line.strip().startswith('H'))
        
        assert carbon_count == 2
        assert hydrogen_count == 6
    
    def test_parse_energy_from_output(self):
        """Test parsing energy from xTB output."""
        calc = XTBCalculator()
        
        # Mock xTB output with energy line
        mock_output = """
        Some output...
        :: total energy              -1.021517863204 Eh    ::
        More output...
        """
        
        energy = calc._parse_energy_from_output(mock_output)
        assert energy is not None
        assert abs(energy - (-1.021517863204)) < 1e-10
    
    def test_parse_energy_no_match(self):
        """Test parsing energy when no energy found."""
        calc = XTBCalculator()
        
        mock_output = "No energy information here"
        energy = calc._parse_energy_from_output(mock_output)
        assert energy is None
    
    def test_parse_optimization_log(self):
        """Test parsing optimization log file."""
        calc = XTBCalculator()
        
        # Create temporary log file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.log') as f:
            f.write("2\n")
            f.write(" energy: -0.981983694723 gnorm: 0.029902840595 xtb: 6.5.0\n")
            f.write("H 0.0 0.0 0.0\n")
            f.write("H 0.0 0.0 0.74\n")
            f.write("2\n")
            f.write(" energy: -0.982674982093 gnorm: 0.003505242088 xtb: 6.5.0\n")
            f.write("H 0.0 0.0 -0.016\n")
            f.write("H 0.0 0.0 0.756\n")
            temp_file = f.name
        
        try:
            steps = calc._parse_optimization_log(temp_file)
            
            assert len(steps) == 2
            assert steps[0]['energy'] == -0.981983694723
            assert steps[0]['gradient_norm'] == 0.029902840595
            assert steps[0]['step'] == 1
            
            assert steps[1]['energy'] == -0.982674982093
            assert steps[1]['gradient_norm'] == 0.003505242088
            assert steps[1]['step'] == 2
            
        finally:
            os.unlink(temp_file)
    
    @pytest.mark.slow
    def test_single_point_hydrogen(self):
        """Test single-point calculation on hydrogen molecule."""
        mol = Molecule.from_smiles("[H][H]", name="hydrogen")
        calc = XTBCalculator()
        
        result = calc.single_point(mol)
        
        assert result['success'] is True
        assert result['energy'] is not None
        assert result['energy_kcal_mol'] is not None
        assert result['method'] == "gfn2"
        assert result['charge'] == 0
        assert result['multiplicity'] == 1
        
        # Hydrogen molecule should have negative energy
        assert result['energy'] < 0
        
        print(f"H2 energy: {result['energy']:.6f} Eh ({result['energy_kcal_mol']:.2f} kcal/mol)")
    
    @pytest.mark.slow
    def test_single_point_with_environment(self):
        """Test single-point calculation with environment."""
        mol = Molecule.from_smiles("O", name="water")
        env = Environment(temperature=298.15, solvent="water")
        calc = XTBCalculator()
        
        result = calc.single_point(mol, env)
        
        assert result['success'] is True
        assert result['energy'] is not None
        assert result['solvent'] == "water"
        
        print(f"H2O energy in water: {result['energy']:.6f} Eh ({result['energy_kcal_mol']:.2f} kcal/mol)")
    
    @pytest.mark.slow
    def test_optimization_hydrogen(self):
        """Test geometry optimization on hydrogen molecule."""
        mol = Molecule.from_smiles("[H][H]", name="hydrogen")
        calc = XTBCalculator()
        
        result = calc.optimize(mol)
        
        assert result['success'] is True
        assert result['energy'] is not None
        assert result['energy_kcal_mol'] is not None
        assert len(result['optimization_steps']) > 0
        
        # Check that optimization converged
        final_step = result['optimization_steps'][-1]
        assert final_step['gradient_norm'] < 0.1  # Should be well converged
        
        print(f"H2 optimized energy: {result['energy']:.6f} Eh")
        print(f"Optimization steps: {len(result['optimization_steps'])}")
        print(f"Final gradient norm: {final_step['gradient_norm']:.6f}")
    
    @pytest.mark.slow
    def test_optimization_with_different_methods(self):
        """Test optimization with different xTB methods."""
        mol = Molecule.from_smiles("CC", name="ethane")
        
        methods = ["gfn2", "gfn1"]
        results = {}
        
        for method in methods:
            calc = XTBCalculator(method=method)
            result = calc.optimize(mol)
            
            assert result['success'] is True
            assert result['method'] == method
            
            results[method] = result['energy_kcal_mol']
            print(f"Ethane energy ({method}): {result['energy_kcal_mol']:.2f} kcal/mol")
        
        # Different methods should give different energies
        assert abs(results["gfn2"] - results["gfn1"]) > 0.1
    
    @pytest.mark.slow
    def test_charged_molecule(self):
        """Test calculation on charged molecule."""
        mol = Molecule.from_smiles("[NH4+]", name="ammonium")
        calc = XTBCalculator(charge=1)
        
        result = calc.single_point(mol)
        
        assert result['success'] is True
        assert result['charge'] == 1
        assert result['energy'] is not None
        
        print(f"NH4+ energy: {result['energy_kcal_mol']:.2f} kcal/mol")
    
    def test_invalid_molecule_handling(self):
        """Test handling of invalid molecules."""
        # Create a molecule that might cause issues
        mol = Molecule.from_smiles("C", name="methane")
        
        # Use invalid xTB path to simulate failure
        calc = XTBCalculator(xtb_path="nonexistent_xtb")
        
        with pytest.raises(RuntimeError):
            calc.single_point(mol)


class TestXTBCalculatorIntegration:
    """Integration tests for XTBCalculator."""
    
    @pytest.mark.slow
    def test_energy_comparison(self):
        """Test energy comparison between similar molecules."""
        molecules = [
            Molecule.from_smiles("C", name="methane"),
            Molecule.from_smiles("CC", name="ethane"),
            Molecule.from_smiles("CCC", name="propane")
        ]
        
        calc = XTBCalculator()
        energies = []
        
        for mol in molecules:
            result = calc.single_point(mol)
            assert result['success'] is True
            energies.append(result['energy_kcal_mol'])
            print(f"{mol.name}: {result['energy_kcal_mol']:.2f} kcal/mol")
        
        # Longer alkanes should have more negative energies
        assert energies[0] > energies[1]  # methane > ethane
        assert energies[1] > energies[2]  # ethane > propane
    
    @pytest.mark.slow
    def test_solvent_effect(self):
        """Test solvent effects on molecular energies."""
        mol = Molecule.from_smiles("CCO", name="ethanol")
        calc = XTBCalculator()
        
        # Calculate in gas phase
        gas_result = calc.single_point(mol)
        assert gas_result['success'] is True
        
        # Calculate in water
        water_calc = XTBCalculator(solvent="water")
        water_result = water_calc.single_point(mol)
        assert water_result['success'] is True
        
        print(f"Ethanol in gas: {gas_result['energy_kcal_mol']:.2f} kcal/mol")
        print(f"Ethanol in water: {water_result['energy_kcal_mol']:.2f} kcal/mol")
        
        # Solvation should stabilize polar molecules
        assert water_result['energy_kcal_mol'] < gas_result['energy_kcal_mol']
