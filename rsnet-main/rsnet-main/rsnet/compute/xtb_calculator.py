"""
xTB calculator for quantum chemical calculations.

This module provides an interface to the xTB (extended Tight Binding) program
for fast semi-empirical quantum chemical calculations.
"""

import os
import tempfile
import subprocess
import shutil
from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path
import re

from rdkit import Chem
from rdkit.Chem import AllChem

from ..core.molecule import Molecule
from ..core.environment import Environment


class XTBCalculator:
    """
    Interface to xTB (extended Tight Binding) calculations.
    
    This class provides methods to perform single-point energy calculations,
    geometry optimizations, and frequency calculations using xTB.
    """
    
    def __init__(self, 
                 xtb_path: str = "xtb",
                 method: str = "gfn2",
                 charge: int = 0,
                 multiplicity: int = 1,
                 solvent: Optional[str] = None,
                 temp_dir: Optional[str] = None):
        """
        Initialize xTB calculator.
        
        Args:
            xtb_path: Path to xTB executable
            method: xTB method (gfn0, gfn1, gfn2, gfnff)
            charge: Molecular charge
            multiplicity: Spin multiplicity (2S+1)
            solvent: Solvent for implicit solvation (water, acetone, etc.)
            temp_dir: Directory for temporary files
        """
        self.xtb_path = xtb_path
        self.method = method
        self.charge = charge
        self.multiplicity = multiplicity
        self.solvent = solvent
        self.temp_dir = temp_dir or tempfile.gettempdir()
        
        # Check if xTB is available
        self._check_xtb_availability()
    
    def _check_xtb_availability(self):
        """Check if xTB executable is available."""
        try:
            result = subprocess.run([self.xtb_path, "--help"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode != 0:
                raise RuntimeError(f"xTB executable not working: {self.xtb_path}")
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            raise RuntimeError(f"xTB executable not found or not working: {self.xtb_path}") from e
    
    def _molecule_to_xyz(self, molecule: Molecule) -> str:
        """
        Convert molecule to XYZ format string.
        
        Args:
            molecule: Molecule object
            
        Returns:
            XYZ format string
        """
        # Get 3D coordinates (generate if not available)
        if not molecule.has_conformer():
            molecule.generate_conformer()
        
        rdkit_mol = molecule.rdkit_mol
        conf = rdkit_mol.GetConformer()
        
        # Build XYZ string
        xyz_lines = [str(rdkit_mol.GetNumAtoms())]
        xyz_lines.append(f"{molecule.name} - {molecule.smiles}")
        
        for i in range(rdkit_mol.GetNumAtoms()):
            atom = rdkit_mol.GetAtomWithIdx(i)
            pos = conf.GetAtomPosition(i)
            symbol = atom.GetSymbol()
            xyz_lines.append(f"{symbol:2s} {pos.x:12.8f} {pos.y:12.8f} {pos.z:12.8f}")
        
        return "\n".join(xyz_lines)
    
    def _parse_energy_from_output(self, output: str) -> Optional[float]:
        """
        Parse total energy from xTB output.
        
        Args:
            output: xTB output text
            
        Returns:
            Total energy in Hartree, or None if not found
        """
        # Look for the final total energy line
        energy_pattern = r":: total energy\s+(-?\d+\.\d+)\s+Eh"
        matches = re.findall(energy_pattern, output)
        
        if matches:
            return float(matches[-1])  # Return the last (final) energy
        
        return None
    
    def _parse_optimization_log(self, log_file: str) -> List[Dict]:
        """
        Parse optimization trajectory from xtbopt.log file.
        
        Args:
            log_file: Path to xtbopt.log file
            
        Returns:
            List of optimization steps with energy and gradient norm
        """
        if not os.path.exists(log_file):
            return []
        
        steps = []
        with open(log_file, 'r') as f:
            lines = f.readlines()
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if line.startswith('energy:'):
                # Parse energy line: "energy: -0.981983694723 gnorm: 0.029902840595 xtb: 6.5.0"
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        energy = float(parts[1])
                        gnorm = float(parts[3])
                        steps.append({
                            'energy': energy,
                            'gradient_norm': gnorm,
                            'step': len(steps) + 1
                        })
                    except ValueError:
                        pass
            i += 1
        
        return steps
    
    def single_point(self, molecule: Molecule, environment: Optional[Environment] = None) -> Dict:
        """
        Perform single-point energy calculation.
        
        Args:
            molecule: Molecule to calculate
            environment: Environmental conditions
            
        Returns:
            Dictionary with calculation results
        """
        with tempfile.TemporaryDirectory(dir=self.temp_dir) as temp_dir:
            # Write XYZ file
            xyz_content = self._molecule_to_xyz(molecule)
            xyz_file = os.path.join(temp_dir, "input.xyz")
            with open(xyz_file, 'w') as f:
                f.write(xyz_content)
            
            # Build xTB command
            cmd = [self.xtb_path, "input.xyz"]
            
            # Add method
            if self.method != "gfn2":  # gfn2 is default
                if self.method == "gfnff":
                    cmd.append("--gfnff")
                else:
                    cmd.extend(["--gfn", self.method[-1]])  # Extract number from gfn1, gfn2, etc.
            
            # Add charge and multiplicity
            if self.charge != 0:
                cmd.extend(["--chrg", str(self.charge)])
            
            if self.multiplicity != 1:
                uhf = self.multiplicity - 1
                cmd.extend(["--uhf", str(uhf)])
            
            # Add solvent
            if self.solvent and self.solvent != "gas":
                cmd.extend(["--alpb", self.solvent])
            elif environment and environment.solvent != "gas":
                cmd.extend(["--alpb", environment.solvent])
            
            # Run calculation
            try:
                result = subprocess.run(cmd, cwd=temp_dir, capture_output=True, 
                                      text=True, timeout=300)
                
                if result.returncode != 0:
                    raise RuntimeError(f"xTB calculation failed: {result.stderr}")
                
                # Parse results
                energy = self._parse_energy_from_output(result.stdout)
                
                return {
                    'energy': energy,  # Hartree
                    'energy_kcal_mol': energy * 627.509 if energy else None,  # kcal/mol
                    'success': energy is not None,
                    'method': self.method,
                    'charge': self.charge,
                    'multiplicity': self.multiplicity,
                    'solvent': self.solvent or (environment.solvent if environment else "gas"),
                    'stdout': result.stdout,
                    'stderr': result.stderr
                }
                
            except subprocess.TimeoutExpired:
                raise RuntimeError("xTB calculation timed out")
    
    def optimize(self, molecule: Molecule, environment: Optional[Environment] = None) -> Dict:
        """
        Perform geometry optimization.
        
        Args:
            molecule: Molecule to optimize
            environment: Environmental conditions
            
        Returns:
            Dictionary with optimization results including optimized molecule
        """
        with tempfile.TemporaryDirectory(dir=self.temp_dir) as temp_dir:
            # Write XYZ file
            xyz_content = self._molecule_to_xyz(molecule)
            xyz_file = os.path.join(temp_dir, "input.xyz")
            with open(xyz_file, 'w') as f:
                f.write(xyz_content)
            
            # Build xTB command for optimization
            cmd = [self.xtb_path, "input.xyz", "--opt"]
            
            # Add method
            if self.method != "gfn2":
                if self.method == "gfnff":
                    cmd.append("--gfnff")
                else:
                    cmd.extend(["--gfn", self.method[-1]])
            
            # Add charge and multiplicity
            if self.charge != 0:
                cmd.extend(["--chrg", str(self.charge)])
            
            if self.multiplicity != 1:
                uhf = self.multiplicity - 1
                cmd.extend(["--uhf", str(uhf)])
            
            # Add solvent
            if self.solvent and self.solvent != "gas":
                cmd.extend(["--alpb", self.solvent])
            elif environment and environment.solvent != "gas":
                cmd.extend(["--alpb", environment.solvent])
            
            # Run optimization
            try:
                result = subprocess.run(cmd, cwd=temp_dir, capture_output=True, 
                                      text=True, timeout=600)
                
                if result.returncode != 0:
                    raise RuntimeError(f"xTB optimization failed: {result.stderr}")
                
                # Parse final energy
                energy = self._parse_energy_from_output(result.stdout)
                
                # Parse optimization trajectory
                opt_log = os.path.join(temp_dir, "xtbopt.log")
                opt_steps = self._parse_optimization_log(opt_log)
                
                # Read optimized geometry
                opt_xyz_file = os.path.join(temp_dir, "xtbopt.xyz")
                optimized_molecule = None
                
                if os.path.exists(opt_xyz_file):
                    # Read the last geometry from xtbopt.xyz
                    with open(opt_xyz_file, 'r') as f:
                        lines = f.readlines()
                    
                    # Find the last geometry block
                    last_geom_start = -1
                    for i in range(len(lines)):
                        if lines[i].strip().isdigit():
                            last_geom_start = i
                    
                    if last_geom_start >= 0:
                        # Create optimized molecule with new coordinates
                        optimized_molecule = self._create_molecule_from_xyz_block(
                            lines[last_geom_start:], molecule)
                
                return {
                    'energy': energy,  # Hartree
                    'energy_kcal_mol': energy * 627.509 if energy else None,  # kcal/mol
                    'optimized_molecule': optimized_molecule,
                    'optimization_steps': opt_steps,
                    'converged': len(opt_steps) > 0 and opt_steps[-1]['gradient_norm'] < 0.01,
                    'success': energy is not None,
                    'method': self.method,
                    'charge': self.charge,
                    'multiplicity': self.multiplicity,
                    'solvent': self.solvent or (environment.solvent if environment else "gas"),
                    'stdout': result.stdout,
                    'stderr': result.stderr
                }
                
            except subprocess.TimeoutExpired:
                raise RuntimeError("xTB optimization timed out")
    
    def _create_molecule_from_xyz_block(self, xyz_lines: List[str], template_molecule: Molecule) -> Molecule:
        """
        Create a molecule from XYZ coordinate block using template molecule.
        
        Args:
            xyz_lines: Lines from XYZ file starting with atom count
            template_molecule: Template molecule for connectivity
            
        Returns:
            New molecule with updated coordinates
        """
        try:
            # Parse XYZ block
            num_atoms = int(xyz_lines[0].strip())
            
            # Create new conformer for the template molecule
            new_mol = Chem.Mol(template_molecule.rdkit_mol)
            new_mol = Chem.AddHs(new_mol)
            
            # Create new conformer
            conf = Chem.Conformer(new_mol.GetNumAtoms())
            
            # Read coordinates (skip first two lines: atom count and comment)
            coord_lines = xyz_lines[2:2+num_atoms]
            
            for i, line in enumerate(coord_lines):
                parts = line.strip().split()
                if len(parts) >= 4:
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                    conf.SetAtomPosition(i, (x, y, z))
            
            new_mol.AddConformer(conf)
            
            # Create new Molecule object
            return Molecule(new_mol, name=f"{template_molecule.name}_optimized")
            
        except Exception as e:
            print(f"Warning: Could not create optimized molecule: {e}")
            return template_molecule
