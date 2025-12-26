#!/usr/bin/env python3
"""
Basic usage example for RSNet.

This example demonstrates how to create molecules and environments
using the RSNet package.
"""

import sys
import os

# Add the parent directory to the path so we can import rsnet
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment


def main():
    """Demonstrate basic RSNet functionality."""
    
    print("=== RSNet Basic Usage Example ===\n")
    
    # Create some molecules
    print("1. Creating molecules from SMILES:")
    
    molecules = [
        ("CCO", "ethanol"),
        ("CC", "ethane"), 
        ("C", "methane"),
        ("c1ccccc1", "benzene"),
        ("CC(=O)O", "acetic acid")
    ]
    
    mol_objects = []
    for smiles, name in molecules:
        mol = Molecule.from_smiles(smiles, name=name)
        mol_objects.append(mol)
        print(f"  {name}: {mol.formula} (MW: {mol.molecular_weight:.2f} g/mol)")
        print(f"    SMILES: {mol.smiles}")
        print(f"    Atoms: {mol.num_atoms}, Heavy atoms: {mol.num_heavy_atoms}")
        print(f"    Hash: {mol.get_hash()[:8]}...")
        print()
    
    # Create reaction environments
    print("2. Creating reaction environments:")
    
    environments = [
        Environment(temperature=298.15, pressure=1.0, solvent="gas"),
        Environment(temperature=373.15, pressure=1.0, solvent="water"),
        Environment(temperature=500.0, pressure=10.0, solvent="gas"),
    ]
    
    for i, env in enumerate(environments, 1):
        print(f"  Environment {i}: {env}")
        print(f"    Thermal energy: {env.get_thermal_energy():.4f} kcal/mol")
        print(f"    Gas phase: {env.is_gas_phase}")
        print()
    
    # Demonstrate 3D conformer generation
    print("3. Generating 3D conformers:")
    
    ethanol = mol_objects[0]  # ethanol
    print(f"  Molecule: {ethanol.name}")
    
    # Generate conformer
    success = ethanol.generate_conformer()
    if success:
        coords = ethanol.get_coordinates()
        print(f"  Conformer generated successfully!")
        print(f"  Coordinates shape: {coords.shape}")
        print(f"  First atom coordinates: {coords[0]}")
    else:
        print("  Failed to generate conformer")
    
    print()
    
    # Demonstrate molecule comparison
    print("4. Molecule comparison:")
    
    ethanol1 = Molecule.from_smiles("CCO")
    ethanol2 = Molecule.from_smiles("CCO", name="ethanol_copy")
    methanol = Molecule.from_smiles("CO")
    
    print(f"  ethanol1 == ethanol2: {ethanol1 == ethanol2}")
    print(f"  ethanol1 == methanol: {ethanol1 == methanol}")
    print(f"  ethanol1 hash: {ethanol1.get_hash()[:8]}...")
    print(f"  ethanol2 hash: {ethanol2.get_hash()[:8]}...")
    print(f"  methanol hash: {methanol.get_hash()[:8]}...")
    
    print("\n=== Example completed successfully! ===")


if __name__ == "__main__":
    main()
