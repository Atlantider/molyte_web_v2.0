#!/usr/bin/env python3
"""
Bond breaking operator demonstration.

This example shows how to use the bond breaking operator
to generate fragmentation reactions from molecules.
"""

import sys
import os

# Add the parent directory to the path so we can import rsnet
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.operators.bond_breaking import BondBreakingOperator
from rsnet.utils.structure_analysis import get_structure_tags


def analyze_molecule_for_breaking(mol: Molecule):
    """Analyze a molecule's structure for bond breaking potential."""
    print(f"\n=== Analyzing {mol.name} ({mol.smiles}) ===")
    
    # Get structure tags
    tags = get_structure_tags(mol)
    
    print(f"Formula: {mol.formula}")
    print(f"Molecular weight: {mol.molecular_weight:.2f} g/mol")
    print(f"Number of atoms: {mol.num_atoms}")
    print(f"Number of heavy atoms: {mol.num_heavy_atoms}")
    
    rdkit_mol = mol.rdkit_mol
    print(f"Number of bonds: {rdkit_mol.GetNumBonds()}")
    
    print("\nStructural features:")
    print(f"  Has weak bonds: {tags['has_weak_bonds']}")
    print(f"  Has small rings: {tags['has_small_rings']}")
    print(f"  Number of weak bonds: {len(tags['weak_bonds'])}")
    
    if tags['weak_bonds']:
        print("  Weak bond pairs:")
        for bond_pair in tags['weak_bonds']:
            atom1 = rdkit_mol.GetAtomWithIdx(bond_pair[0])
            atom2 = rdkit_mol.GetAtomWithIdx(bond_pair[1])
            print(f"    {atom1.GetSymbol()}{bond_pair[0]}-{atom2.GetSymbol()}{bond_pair[1]}")
    
    if tags['small_rings']:
        print(f"  Small rings: {len(tags['small_rings'])}")
        for i, ring in enumerate(tags['small_rings']):
            print(f"    Ring {i+1}: {len(ring)} atoms")
    
    return tags


def test_bond_breaking_operator(mol: Molecule, env: Environment, config=None):
    """Test the bond breaking operator on a molecule."""
    print(f"\n=== Testing bond breaking on {mol.name} ===")
    print(f"Environment: T={env.temperature}K, P={env.pressure}atm")
    
    # Create operator
    op = BondBreakingOperator(config=config)
    
    # Check if applicable
    applicable = op.is_applicable([mol])
    print(f"Operator applicable: {applicable}")
    
    if not applicable:
        print("Skipping - operator not applicable")
        return []
    
    # Find breakable bonds
    breakable_bonds = op._find_breakable_bonds(mol, env)
    print(f"Breakable bonds found: {len(breakable_bonds)}")
    
    # Apply operator
    reactions = op.apply([mol], env)
    print(f"Generated reactions: {len(reactions)}")
    
    for i, rxn in enumerate(reactions, 1):
        print(f"\nReaction {i}: {rxn.name}")
        print(f"  Reactant: {rxn.reactants[0].smiles}")
        print(f"  Products: {len(rxn.products)} fragments")
        for j, product in enumerate(rxn.products, 1):
            print(f"    Fragment {j}: {product.smiles} ({product.formula}, {product.num_heavy_atoms} heavy atoms)")
        if rxn.reaction_energy is not None:
            print(f"  ΔE: {rxn.reaction_energy:.2f} kcal/mol")
        if rxn.activation_energy is not None:
            print(f"  Ea: {rxn.activation_energy:.2f} kcal/mol")
    
    return reactions


def main():
    """Demonstrate bond breaking operator."""
    
    print("=== RSNet Bond Breaking Operator Demo ===")
    
    # Test molecules with different bond breaking potential
    test_cases = [
        # (smiles, name, temperature, config)
        ("CC", "ethane", 500.0, {'min_fragment_size': 1}),  # Simple C-C bond
        ("CCC", "propane", 500.0, {'min_fragment_size': 1}),  # Multiple C-C bonds
        ("C1CC1", "cyclopropane", 298.15, {'min_fragment_size': 1}),  # Ring strain
        ("C1CCC1", "cyclobutane", 298.15, {'min_fragment_size': 1}),  # Ring strain
        ("CCCCCC", "hexane", 600.0, {'min_fragment_size': 2}),  # Long chain
        ("CC(C)C", "isobutane", 500.0, {'min_fragment_size': 1}),  # Branched
        ("CCO", "ethanol", 400.0, {'min_fragment_size': 1}),  # With heteroatom
        ("c1ccccc1", "benzene", 800.0, {'min_fragment_size': 1}),  # Aromatic (stable)
    ]
    
    all_reactions = []
    
    for smiles, name, temperature, config in test_cases:
        try:
            # Create molecule and environment
            mol = Molecule.from_smiles(smiles, name=name)
            env = Environment(temperature=temperature, pressure=1.0, solvent="gas")
            
            # Analyze structure
            tags = analyze_molecule_for_breaking(mol)
            
            # Test operator
            reactions = test_bond_breaking_operator(mol, env, config)
            all_reactions.extend(reactions)
            
        except Exception as e:
            print(f"Error processing {name}: {e}")
            continue
    
    # Summary
    print(f"\n=== Summary ===")
    print(f"Total molecules tested: {len(test_cases)}")
    print(f"Total reactions generated: {len(all_reactions)}")
    
    if all_reactions:
        print("\nAll generated reactions:")
        for i, rxn in enumerate(all_reactions, 1):
            reactant_smiles = " + ".join([r.smiles for r in rxn.reactants])
            product_smiles = " + ".join([p.smiles for p in rxn.products])
            print(f"  {i}. {reactant_smiles} → {product_smiles}")
    
    # Demonstrate temperature effect
    print(f"\n=== Temperature Effect Demo ===")
    mol = Molecule.from_smiles("CC", name="ethane")
    config = {'min_fragment_size': 1}
    
    temperatures = [298.15, 400.0, 500.0, 600.0, 800.0]
    
    for temp in temperatures:
        env = Environment(temperature=temp)
        op = BondBreakingOperator(config=config)
        breakable_bonds = op._find_breakable_bonds(mol, env)
        reactions = op.apply([mol], env)
        print(f"T={temp}K: {len(breakable_bonds)} breakable bonds, {len(reactions)} reactions")
    
    print("\n=== Demo completed! ===")


if __name__ == "__main__":
    main()
