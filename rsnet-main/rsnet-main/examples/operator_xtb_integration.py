#!/usr/bin/env python3
"""
Operator + xTB Integration Demo

This example demonstrates the complete workflow:
1. Generate reactions using operators
2. Screen reactions using xTB calculations
3. Filter feasible reactions
4. Analyze results
"""

import sys
import os
import time

# Add the parent directory to the path so we can import rsnet
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.operators.bond_breaking import BondBreakingOperator
from rsnet.compute.reaction_screener import ReactionScreener


def main():
    """Demonstrate operator + xTB integration."""
    print("=== RSNet Operator + xTB Integration Demo ===")
    print("This demo shows the complete workflow from reaction generation")
    print("to thermodynamic screening using quantum calculations.\n")
    
    # Test molecules
    molecules = [
        Molecule.from_smiles("CC", name="ethane"),
        Molecule.from_smiles("CCO", name="ethanol"),
        Molecule.from_smiles("CCC", name="propane"),
    ]
    
    # Environment conditions
    env = Environment(temperature=500.0, pressure=1.0, solvent="gas")
    print(f"Environment: T={env.temperature}K, P={env.pressure}atm, solvent={env.solvent}")
    
    # Step 1: Generate reactions using bond breaking operator
    print(f"\n--- Step 1: Generating Reactions ---")
    bond_break_op = BondBreakingOperator(config={'min_fragment_size': 1})
    
    all_reactions = []
    for mol in molecules:
        print(f"Processing {mol.name} ({mol.smiles})...")
        
        if bond_break_op.is_applicable([mol]):
            reactions = bond_break_op.apply([mol], env)
            all_reactions.extend(reactions)
            print(f"  Generated {len(reactions)} reactions")
            
            # Show first reaction as example
            if reactions:
                reaction = reactions[0]
                reactant_smiles = " + ".join([r.smiles for r in reaction.reactants])
                product_smiles = " + ".join([p.smiles for p in reaction.products])
                print(f"  Example: {reactant_smiles} → {product_smiles}")
        else:
            print(f"  No applicable reactions")
    
    print(f"\nTotal reactions generated: {len(all_reactions)}")
    
    # Step 2: Screen reactions using xTB
    print(f"\n--- Step 2: Screening Reactions with xTB ---")
    screener = ReactionScreener(
        max_workers=2,
        energy_threshold=80.0,  # kcal/mol
        optimize_geometries=False  # Faster for demo
    )
    
    # Limit to first 3 reactions for demo speed
    reactions_to_screen = all_reactions[:3]
    print(f"Screening {len(reactions_to_screen)} reactions...")
    
    start_time = time.time()
    screening_results = screener.screen_reactions(reactions_to_screen, env)
    screening_time = time.time() - start_time
    
    print(f"Screening completed in {screening_time:.2f} seconds")
    
    # Step 3: Analyze results
    print(f"\n--- Step 3: Analysis Results ---")
    feasible_reactions = []
    
    for i, result in enumerate(screening_results, 1):
        if result['success']:
            reaction = result['reaction']
            reactant_smiles = " + ".join([r.smiles for r in reaction.reactants])
            product_smiles = " + ".join([p.smiles for p in reaction.products])
            
            print(f"\nReaction {i}: {reactant_smiles} → {product_smiles}")
            print(f"  ΔE: {result['reaction_energy']:8.2f} kcal/mol")
            print(f"  Ea: {result['activation_energy']:8.2f} kcal/mol")
            print(f"  Feasible: {result['feasible']}")
            print(f"  Rate constant: {result['rate_constant']:.2e} s⁻¹")
            
            if result['feasible']:
                feasible_reactions.append(reaction)
                print(f"  ✓ This reaction is feasible!")
            else:
                print(f"  ✗ This reaction is not feasible")
                
        else:
            print(f"\nReaction {i}: Failed - {result.get('error', 'Unknown error')}")
    
    # Step 4: Summary
    print(f"\n--- Step 4: Summary ---")
    print(f"Molecules processed: {len(molecules)}")
    print(f"Total reactions generated: {len(all_reactions)}")
    print(f"Reactions screened: {len(screening_results)}")
    print(f"Feasible reactions: {len(feasible_reactions)}")
    
    if len(screening_results) > 0:
        success_rate = len([r for r in screening_results if r['success']]) / len(screening_results)
        feasibility_rate = len(feasible_reactions) / len(screening_results)
        print(f"Calculation success rate: {success_rate*100:.1f}%")
        print(f"Feasibility rate: {feasibility_rate*100:.1f}%")
    
    # Energy cache statistics
    stats = screener.get_energy_statistics()
    if stats['cache_size'] > 0:
        print(f"\n--- Energy Cache Statistics ---")
        print(f"Unique molecules calculated: {stats['cache_size']}")
        print(f"Energy range: {stats['min_energy']:.2f} to {stats['max_energy']:.2f} kcal/mol")
        print(f"Mean energy: {stats['mean_energy']:.2f} kcal/mol")
    
    # Step 5: Demonstrate filtering
    if feasible_reactions:
        print(f"\n--- Step 5: Feasible Reaction Details ---")
        for i, reaction in enumerate(feasible_reactions, 1):
            reactant_smiles = " + ".join([r.smiles for r in reaction.reactants])
            product_smiles = " + ".join([p.smiles for p in reaction.products])
            print(f"Feasible Reaction {i}: {reactant_smiles} → {product_smiles}")
            print(f"  Template: {reaction.template}")
            print(f"  Operator: {reaction.operator}")
    
    print(f"\n=== Integration Demo Completed Successfully! ===")
    print("The RSNet system can now:")
    print("✓ Generate reactions using chemical operators")
    print("✓ Calculate reaction energies using xTB")
    print("✓ Estimate activation energies")
    print("✓ Screen for thermodynamic and kinetic feasibility")
    print("✓ Filter and rank reactions")
    print("\nReady for reaction network generation!")


if __name__ == "__main__":
    main()
