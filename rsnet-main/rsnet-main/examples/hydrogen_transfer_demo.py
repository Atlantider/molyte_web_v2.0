#!/usr/bin/env python3
"""
Hydrogen transfer operator demonstration.

This example shows how to use the hydrogen transfer operator
to generate reactions from molecules.
"""

import sys
import os

# Add the parent directory to the path so we can import rsnet
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.operators.hydrogen_transfer import HydrogenTransferOperator
from rsnet.utils.structure_analysis import get_structure_tags, find_hydrogen_transfer_sites


def analyze_molecule(mol: Molecule):
    """Analyze a molecule's structure for hydrogen transfer potential."""
    print(f"\n=== Analyzing {mol.name} ({mol.smiles}) ===")
    
    # Get structure tags
    tags = get_structure_tags(mol)
    
    print(f"Formula: {mol.formula}")
    print(f"Molecular weight: {mol.molecular_weight:.2f} g/mol")
    print(f"Number of atoms: {mol.num_atoms}")
    print(f"Number of heavy atoms: {mol.num_heavy_atoms}")
    
    print("\nStructural features:")
    print(f"  Has heteroatoms: {tags['has_heteroatoms']}")
    print(f"  Has acidic hydrogens: {tags['has_acidic_hydrogens']}")
    print(f"  Has polar bonds: {tags['has_polar_bonds']}")
    print(f"  Has H-bond acceptors: {len(tags['h_bond_acceptors'])} acceptors")
    print(f"  Has H-bond donors: {tags['num_hbd']} donors")
    
    if tags['heteroatoms']:
        rdkit_mol = mol.rdkit_mol
        hetero_symbols = [rdkit_mol.GetAtomWithIdx(idx).GetSymbol() 
                         for idx in tags['heteroatoms']]
        print(f"  Heteroatoms: {hetero_symbols}")
    
    # Find hydrogen transfer sites
    sites = find_hydrogen_transfer_sites(mol)
    print(f"\nPotential H-transfer sites: {len(sites)}")
    
    for i, site in enumerate(sites, 1):
        print(f"  Site {i}: {site['donor_element']}-H → {site['acceptor_element']}")
        if site['distance'] is not None:
            print(f"    Distance: {site['distance']:.2f} Å")
    
    return tags, sites


def test_hydrogen_transfer_operator(mol: Molecule, env: Environment):
    """Test the hydrogen transfer operator on a molecule."""
    print(f"\n=== Testing H-transfer operator on {mol.name} ===")
    
    # Create operator
    op = HydrogenTransferOperator()
    
    # Check if applicable
    applicable = op.is_applicable([mol])
    print(f"Operator applicable: {applicable}")
    
    if not applicable:
        print("Skipping - operator not applicable")
        return []
    
    # Apply operator
    reactions = op.apply([mol], env)
    print(f"Generated reactions: {len(reactions)}")
    
    for i, rxn in enumerate(reactions, 1):
        print(f"\nReaction {i}: {rxn.name}")
        print(f"  Reactants: {[r.smiles for r in rxn.reactants]}")
        print(f"  Products: {[p.smiles for p in rxn.products]}")
        if rxn.reaction_energy is not None:
            print(f"  ΔE: {rxn.reaction_energy:.2f} kcal/mol")
        if rxn.activation_energy is not None:
            print(f"  Ea: {rxn.activation_energy:.2f} kcal/mol")
    
    return reactions


def main():
    """Demonstrate hydrogen transfer operator."""
    
    print("=== RSNet Hydrogen Transfer Operator Demo ===")
    
    # Create reaction environment
    env = Environment(temperature=298.15, pressure=1.0, solvent="gas")
    print(f"Environment: {env}")
    print(f"Thermal energy: {env.get_thermal_energy():.4f} kcal/mol")
    
    # Test molecules with different H-transfer potential
    test_molecules = [
        ("CCO", "ethanol"),
        ("CO", "methanol"), 
        ("CC(=O)O", "acetic_acid"),
        ("CC(C)O", "isopropanol"),
        ("c1ccc(O)cc1", "phenol"),
        ("CC", "ethane"),  # No heteroatoms - should not be applicable
        ("CCN", "ethylamine"),
        ("CC(=O)N", "acetamide"),
    ]
    
    all_reactions = []
    
    for smiles, name in test_molecules:
        try:
            # Create molecule
            mol = Molecule.from_smiles(smiles, name=name)
            
            # Analyze structure
            tags, sites = analyze_molecule(mol)
            
            # Test operator
            reactions = test_hydrogen_transfer_operator(mol, env)
            all_reactions.extend(reactions)
            
        except Exception as e:
            print(f"Error processing {name}: {e}")
            continue
    
    # Summary
    print(f"\n=== Summary ===")
    print(f"Total molecules tested: {len(test_molecules)}")
    print(f"Total reactions generated: {len(all_reactions)}")
    
    if all_reactions:
        print("\nAll generated reactions:")
        for i, rxn in enumerate(all_reactions, 1):
            reactant_smiles = " + ".join([r.smiles for r in rxn.reactants])
            product_smiles = " + ".join([p.smiles for p in rxn.products])
            print(f"  {i}. {reactant_smiles} → {product_smiles}")
    
    print("\n=== Demo completed! ===")


if __name__ == "__main__":
    main()
