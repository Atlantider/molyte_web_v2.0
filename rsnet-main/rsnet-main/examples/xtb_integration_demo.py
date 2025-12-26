#!/usr/bin/env python3
"""
xTB integration demonstration.

This example shows how to use xTB calculations for:
1. Single-point energy calculations
2. Geometry optimizations
3. Reaction energy screening
4. Thermodynamic and kinetic analysis
"""

import sys
import os
import time

# Add the parent directory to the path so we can import rsnet
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.core.reaction import Reaction
from rsnet.compute.xtb_calculator import XTBCalculator
from rsnet.compute.reaction_screener import ReactionScreener
from rsnet.operators.bond_breaking import BondBreakingOperator


def demo_xtb_calculator():
    """Demonstrate xTB calculator functionality."""
    print("=== xTB Calculator Demo ===")
    
    # Test molecules
    molecules = [
        Molecule.from_smiles("[H][H]", name="hydrogen"),
        Molecule.from_smiles("C", name="methane"),
        Molecule.from_smiles("CC", name="ethane"),
        Molecule.from_smiles("CCO", name="ethanol"),
        Molecule.from_smiles("O", name="water"),
    ]
    
    calc = XTBCalculator()
    
    print("\n--- Single-Point Energy Calculations ---")
    energies = {}
    
    for mol in molecules:
        print(f"\nCalculating {mol.name} ({mol.smiles})...")
        start_time = time.time()
        
        try:
            result = calc.single_point(mol)
            calc_time = time.time() - start_time
            
            if result['success']:
                energy = result['energy_kcal_mol']
                energies[mol.name] = energy
                print(f"  Energy: {energy:.2f} kcal/mol")
                print(f"  Time: {calc_time:.2f} seconds")
            else:
                print(f"  Calculation failed")
                
        except Exception as e:
            print(f"  Error: {e}")
    
    print(f"\n--- Energy Summary ---")
    for name, energy in sorted(energies.items(), key=lambda x: x[1]):
        print(f"  {name:10s}: {energy:8.2f} kcal/mol")
    
    # Test geometry optimization
    print(f"\n--- Geometry Optimization ---")
    mol = Molecule.from_smiles("CCO", name="ethanol")
    print(f"Optimizing {mol.name}...")
    
    start_time = time.time()
    result = calc.optimize(mol)
    opt_time = time.time() - start_time
    
    if result['success']:
        print(f"  Initial energy: {result['energy_kcal_mol']:.2f} kcal/mol")
        print(f"  Optimization steps: {len(result['optimization_steps'])}")
        print(f"  Converged: {result['converged']}")
        print(f"  Time: {opt_time:.2f} seconds")
        
        if result['optimization_steps']:
            final_step = result['optimization_steps'][-1]
            print(f"  Final gradient norm: {final_step['gradient_norm']:.6f}")
    else:
        print(f"  Optimization failed")
    
    # Test solvent effects
    print(f"\n--- Solvent Effects ---")
    mol = Molecule.from_smiles("CCO", name="ethanol")
    
    solvents = ["gas", "water", "acetone"]
    solvent_energies = {}
    
    for solvent in solvents:
        if solvent == "gas":
            calc_solv = XTBCalculator()
        else:
            calc_solv = XTBCalculator(solvent=solvent)
        
        try:
            result = calc_solv.single_point(mol)
            if result['success']:
                solvent_energies[solvent] = result['energy_kcal_mol']
                print(f"  {solvent:10s}: {result['energy_kcal_mol']:8.2f} kcal/mol")
        except Exception as e:
            print(f"  {solvent:10s}: Error - {e}")
    
    # Calculate solvation energies
    if "gas" in solvent_energies:
        print(f"\n--- Solvation Energies ---")
        gas_energy = solvent_energies["gas"]
        for solvent, energy in solvent_energies.items():
            if solvent != "gas":
                solvation_energy = energy - gas_energy
                print(f"  ΔG_solv({solvent}): {solvation_energy:6.2f} kcal/mol")


def demo_reaction_screening():
    """Demonstrate reaction screening functionality."""
    print(f"\n\n=== Reaction Screening Demo ===")
    
    # Generate some reactions using bond breaking operator
    print("\n--- Generating Reactions ---")
    
    molecules = [
        Molecule.from_smiles("CC", name="ethane"),
        Molecule.from_smiles("CCO", name="ethanol"),
        Molecule.from_smiles("CCC", name="propane"),
    ]
    
    env = Environment(temperature=500.0, pressure=1.0, solvent="gas")
    bond_break_op = BondBreakingOperator(config={'min_fragment_size': 1})
    
    all_reactions = []
    for mol in molecules:
        if bond_break_op.is_applicable([mol]):
            reactions = bond_break_op.apply([mol], env)
            all_reactions.extend(reactions)
            print(f"  {mol.name}: {len(reactions)} reactions")
    
    print(f"Total reactions generated: {len(all_reactions)}")
    
    # Screen reactions
    print(f"\n--- Screening Reactions ---")
    screener = ReactionScreener(
        max_workers=2,
        energy_threshold=100.0,  # Generous threshold for demo
        optimize_geometries=False  # Faster for demo
    )
    
    start_time = time.time()
    screening_results = screener.screen_reactions(all_reactions[:5], env)  # Limit for demo
    screening_time = time.time() - start_time
    
    print(f"Screened {len(screening_results)} reactions in {screening_time:.2f} seconds")
    
    # Analyze results
    print(f"\n--- Screening Results ---")
    feasible_count = 0
    
    for i, result in enumerate(screening_results, 1):
        if result['success']:
            reaction = result['reaction']
            reactant_smiles = " + ".join([r.smiles for r in reaction.reactants])
            product_smiles = " + ".join([p.smiles for p in reaction.products])
            
            print(f"\nReaction {i}: {reactant_smiles} → {product_smiles}")
            print(f"  ΔE: {result['reaction_energy']:6.2f} kcal/mol")
            print(f"  Ea: {result['activation_energy']:6.2f} kcal/mol")
            print(f"  Feasible: {result['feasible']}")
            print(f"  Rate constant: {result['rate_constant']:.2e} s⁻¹")
            
            if result['feasible']:
                feasible_count += 1
        else:
            print(f"\nReaction {i}: Failed - {result.get('error', 'Unknown error')}")
    
    print(f"\n--- Summary ---")
    print(f"Total reactions screened: {len(screening_results)}")
    print(f"Feasible reactions: {feasible_count}")
    print(f"Success rate: {feasible_count/len(screening_results)*100:.1f}%")
    
    # Energy statistics
    stats = screener.get_energy_statistics()
    if stats['cache_size'] > 0:
        print(f"\n--- Energy Cache Statistics ---")
        print(f"Molecules calculated: {stats['cache_size']}")
        print(f"Energy range: {stats['min_energy']:.2f} to {stats['max_energy']:.2f} kcal/mol")
        print(f"Mean energy: {stats['mean_energy']:.2f} kcal/mol")


def demo_thermodynamic_analysis():
    """Demonstrate thermodynamic analysis of reactions."""
    print(f"\n\n=== Thermodynamic Analysis Demo ===")
    
    # Create some test reactions
    reactions = [
        # H2 dissociation (highly endothermic)
        Reaction(
            [Molecule.from_smiles("[H][H]", name="H2")],
            [Molecule.from_smiles("[H]", name="H1"), Molecule.from_smiles("[H]", name="H2")],
            name="H2_dissociation"
        ),
        
        # Methane combustion (simplified - just C + O2 -> CO2)
        Reaction(
            [Molecule.from_smiles("C", name="C"), Molecule.from_smiles("O=O", name="O2")],
            [Molecule.from_smiles("O=C=O", name="CO2")],
            name="carbon_oxidation"
        ),
    ]
    
    screener = ReactionScreener(optimize_geometries=False)
    env = Environment(temperature=298.15, pressure=1.0, solvent="gas")
    
    print(f"\n--- Analyzing {len(reactions)} reactions ---")
    
    for reaction in reactions:
        print(f"\nReaction: {reaction.name}")
        reactant_smiles = " + ".join([r.smiles for r in reaction.reactants])
        product_smiles = " + ".join([p.smiles for p in reaction.products])
        print(f"  {reactant_smiles} → {product_smiles}")
        
        try:
            # Calculate reaction energy
            energetics = screener.calculate_reaction_energy(reaction, env)
            
            if energetics['success']:
                print(f"  ΔE: {energetics['reaction_energy']:6.2f} kcal/mol")
                print(f"  Type: {'Exothermic' if energetics['exothermic'] else 'Endothermic'}")
                
                # Estimate activation energy
                ea = screener.estimate_activation_energy(reaction, env)
                if ea is not None:
                    print(f"  Ea (estimated): {ea:6.2f} kcal/mol")
                    
                    # Calculate rate constant at different temperatures
                    print(f"  Rate constants:")
                    for temp in [298.15, 400.0, 500.0, 600.0]:
                        import math
                        R = 1.987e-3  # kcal/(mol·K)
                        k = 1e13 * math.exp(-ea / (R * temp))
                        print(f"    T={temp:5.1f}K: k={k:.2e} s⁻¹")
            else:
                print(f"  Error: {energetics.get('error', 'Unknown error')}")
                
        except Exception as e:
            print(f"  Error: {e}")


def main():
    """Run all demonstrations."""
    print("=== RSNet xTB Integration Comprehensive Demo ===")
    print("This demo showcases quantum chemical calculations using xTB")
    print("for reaction network analysis.\n")
    
    try:
        # Test xTB availability
        calc = XTBCalculator()
        print("✓ xTB calculator initialized successfully")
        
        # Run demonstrations
        demo_xtb_calculator()
        demo_reaction_screening()
        demo_thermodynamic_analysis()
        
        print(f"\n=== Demo Completed Successfully! ===")
        print("xTB integration is working properly and ready for")
        print("reaction network generation and analysis.")
        
    except Exception as e:
        print(f"❌ Demo failed: {e}")
        print("Please check xTB installation and try again.")


if __name__ == "__main__":
    main()
