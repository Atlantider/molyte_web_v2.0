#!/usr/bin/env python3
"""
Comprehensive operators demonstration.

This example shows how to use multiple reaction operators
to generate reactions from molecules.
"""

import sys
import os

# Add the parent directory to the path so we can import rsnet
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.operators.hydrogen_transfer import HydrogenTransferOperator
from rsnet.operators.bond_breaking import BondBreakingOperator
from rsnet.utils.structure_analysis import get_structure_tags


def test_all_operators(mol: Molecule, env: Environment):
    """Test all available operators on a molecule."""
    print(f"\n=== Testing all operators on {mol.name} ({mol.smiles}) ===")
    print(f"Environment: T={env.temperature}K, P={env.pressure}atm")
    
    # Get structure tags
    tags = get_structure_tags(mol)
    print(f"Molecular formula: {mol.formula}")
    print(f"Heavy atoms: {mol.num_heavy_atoms}")
    print(f"Has heteroatoms: {tags['has_heteroatoms']}")
    print(f"Has acidic hydrogens: {tags['has_acidic_hydrogens']}")
    print(f"Has weak bonds: {tags['has_weak_bonds']}")
    
    all_reactions = []
    
    # Test Hydrogen Transfer Operator
    print(f"\n--- Hydrogen Transfer Operator ---")
    h_transfer_op = HydrogenTransferOperator()
    
    if h_transfer_op.is_applicable([mol]):
        h_reactions = h_transfer_op.apply([mol], env)
        print(f"Generated {len(h_reactions)} hydrogen transfer reactions")
        
        for i, rxn in enumerate(h_reactions, 1):
            print(f"  H-Transfer {i}: {rxn.reactants[0].smiles} → {' + '.join([p.smiles for p in rxn.products])}")
        
        all_reactions.extend(h_reactions)
    else:
        print("Hydrogen transfer not applicable")
    
    # Test Bond Breaking Operator
    print(f"\n--- Bond Breaking Operator ---")
    bond_break_config = {'min_fragment_size': 1}  # Allow small fragments
    bond_break_op = BondBreakingOperator(config=bond_break_config)
    
    if bond_break_op.is_applicable([mol]):
        bond_reactions = bond_break_op.apply([mol], env)
        print(f"Generated {len(bond_reactions)} bond breaking reactions")
        
        for i, rxn in enumerate(bond_reactions, 1):
            products_str = ' + '.join([p.smiles for p in rxn.products])
            print(f"  Bond-Break {i}: {rxn.reactants[0].smiles} → {products_str}")
        
        all_reactions.extend(bond_reactions)
    else:
        print("Bond breaking not applicable")
    
    return all_reactions


def main():
    """Demonstrate all available operators."""
    
    print("=== RSNet Comprehensive Operators Demo ===")
    
    # Test different types of molecules
    test_molecules = [
        # Simple hydrocarbons
        ("CC", "ethane"),
        ("CCC", "propane"),
        ("C1CC1", "cyclopropane"),
        
        # Molecules with heteroatoms
        ("CCO", "ethanol"),
        ("CO", "methanol"),
        ("CC(=O)O", "acetic_acid"),
        ("CCN", "ethylamine"),
        
        # More complex molecules
        ("CC(C)C", "isobutane"),
        ("CCCCCC", "hexane"),
        ("c1ccccc1", "benzene"),
    ]
    
    # Test at different conditions
    environments = [
        Environment(temperature=298.15, pressure=1.0, solvent="gas"),  # Room temperature
        Environment(temperature=500.0, pressure=1.0, solvent="gas"),   # High temperature
        Environment(temperature=298.15, pressure=1.0, solvent="water"), # Aqueous
    ]
    
    all_results = []
    
    for smiles, name in test_molecules:
        try:
            mol = Molecule.from_smiles(smiles, name=name)
            
            for env in environments:
                env_name = f"T={env.temperature}K"
                if env.solvent != "gas":
                    env_name += f", {env.solvent}"
                
                print(f"\n{'='*60}")
                print(f"Testing {name} at {env_name}")
                print(f"{'='*60}")
                
                reactions = test_all_operators(mol, env)
                
                result = {
                    'molecule': name,
                    'smiles': smiles,
                    'environment': env_name,
                    'reactions': len(reactions),
                    'reaction_list': reactions
                }
                all_results.append(result)
                
        except Exception as e:
            print(f"Error processing {name}: {e}")
            continue
    
    # Summary
    print(f"\n{'='*60}")
    print("COMPREHENSIVE SUMMARY")
    print(f"{'='*60}")
    
    total_reactions = sum(r['reactions'] for r in all_results)
    print(f"Total molecules tested: {len(test_molecules)}")
    print(f"Total conditions tested: {len(test_molecules) * len(environments)}")
    print(f"Total reactions generated: {total_reactions}")
    
    # Group by molecule type
    print(f"\n--- Results by Molecule ---")
    molecule_summary = {}
    for result in all_results:
        mol_name = result['molecule']
        if mol_name not in molecule_summary:
            molecule_summary[mol_name] = {'total_reactions': 0, 'conditions': []}
        
        molecule_summary[mol_name]['total_reactions'] += result['reactions']
        molecule_summary[mol_name]['conditions'].append({
            'env': result['environment'],
            'reactions': result['reactions']
        })
    
    for mol_name, summary in molecule_summary.items():
        print(f"\n{mol_name}: {summary['total_reactions']} total reactions")
        for condition in summary['conditions']:
            if condition['reactions'] > 0:
                print(f"  {condition['env']}: {condition['reactions']} reactions")
    
    # Show most reactive molecules
    print(f"\n--- Most Reactive Molecules ---")
    reactive_molecules = [(name, data['total_reactions']) 
                         for name, data in molecule_summary.items() 
                         if data['total_reactions'] > 0]
    reactive_molecules.sort(key=lambda x: x[1], reverse=True)
    
    for mol_name, reaction_count in reactive_molecules[:5]:
        print(f"  {mol_name}: {reaction_count} reactions")
    
    # Show example reactions
    print(f"\n--- Example Reactions ---")
    reaction_examples = []
    for result in all_results:
        for rxn in result['reaction_list']:
            reactant = ' + '.join([r.smiles for r in rxn.reactants])
            product = ' + '.join([p.smiles for p in rxn.products])
            reaction_examples.append(f"{reactant} → {product}")
    
    # Remove duplicates and show first 10
    unique_reactions = list(set(reaction_examples))
    for i, rxn in enumerate(unique_reactions[:10], 1):
        print(f"  {i}. {rxn}")
    
    if len(unique_reactions) > 10:
        print(f"  ... and {len(unique_reactions) - 10} more unique reactions")
    
    print(f"\n=== Demo completed! ===")
    print(f"Generated {len(unique_reactions)} unique reaction types")


if __name__ == "__main__":
    main()
