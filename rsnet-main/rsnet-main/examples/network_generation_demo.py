#!/usr/bin/env python3
"""
Reaction Network Generation Demo

This example demonstrates the complete RSNet workflow:
1. Initialize seed molecules and environment
2. Generate reaction network using operators and xTB screening
3. Analyze network topology and pathways
4. Apply pruning strategies
5. Export results
"""

import sys
import os
import time
import json

# Add the parent directory to the path so we can import rsnet
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.operators.bond_breaking import BondBreakingOperator
from rsnet.operators.hydrogen_transfer import HydrogenTransferOperator
from rsnet.compute.reaction_screener import ReactionScreener
from rsnet.network.generator import NetworkGenerator
from rsnet.network.analyzer import NetworkAnalyzer
from rsnet.network.pruner import NetworkPruner


def main():
    """Demonstrate complete network generation workflow."""
    print("=== RSNet Complete Network Generation Demo ===")
    print("This demo shows the full workflow from molecules to analyzed networks.\n")
    
    # Step 1: Define seed molecules and environment
    print("--- Step 1: Setting up seed molecules and environment ---")
    
    seed_molecules = [
        Molecule.from_smiles("CC", name="ethane"),
        Molecule.from_smiles("CCO", name="ethanol"),
        Molecule.from_smiles("C=C", name="ethene"),
    ]
    
    environment = Environment(
        temperature=500.0,  # High temperature for more reactions
        pressure=1.0,
        solvent="gas"
    )
    
    print(f"Seed molecules: {[mol.name for mol in seed_molecules]}")
    print(f"Environment: T={environment.temperature}K, P={environment.pressure}atm")
    
    # Step 2: Configure network generator
    print(f"\n--- Step 2: Configuring network generator ---")
    
    # Set up operators
    operators = [
        BondBreakingOperator(config={'min_fragment_size': 1}),
        HydrogenTransferOperator()  # Will return empty for now
    ]
    
    # Set up screener (faster settings for demo)
    screener = ReactionScreener(
        max_workers=2,
        energy_threshold=80.0,  # More permissive for demo
        optimize_geometries=False  # Faster
    )
    
    # Configure generator
    config = {
        'max_generations': 3,  # Limit for demo
        'max_species': 50,     # Limit for demo
        'energy_cutoff': 80.0, # kcal/mol
        'min_new_reactions': 1
    }
    
    generator = NetworkGenerator(
        operators=operators,
        screener=screener,
        config=config
    )
    
    print(f"Operators: {[op.name for op in operators]}")
    print(f"Max generations: {config['max_generations']}")
    print(f"Max species: {config['max_species']}")
    
    # Step 3: Generate network
    print(f"\n--- Step 3: Generating reaction network ---")
    
    start_time = time.time()
    network = generator.generate_network(
        seed_molecules=seed_molecules,
        environment=environment,
        max_time=300  # 5 minutes max
    )
    generation_time = time.time() - start_time
    
    print(f"Network generation completed in {generation_time:.2f} seconds")
    
    # Step 4: Analyze network
    print(f"\n--- Step 4: Analyzing network ---")
    
    stats = network.get_statistics()
    print(f"Network statistics:")
    print(f"  Molecules: {stats['num_molecules']}")
    print(f"  Reactions: {stats['num_reactions']}")
    print(f"  Max generation: {stats['max_generation']}")
    print(f"  Molecules by generation: {stats['molecules_by_generation']}")
    
    # Detailed analysis
    analyzer = NetworkAnalyzer(network)
    
    # Topology analysis
    topology = analyzer.analyze_network_topology()
    print(f"\nTopology analysis:")
    print(f"  Connected: {topology['is_connected']}")
    print(f"  Components: {topology['num_components']}")
    print(f"  Average degree: {topology['degree_stats']['mean']:.2f}")
    print(f"  Clustering coefficient: {topology['clustering_coefficient']:.3f}")
    
    # Thermodynamic analysis
    thermo = analyzer.analyze_thermodynamics()
    print(f"\nThermodynamic analysis:")
    print(f"  Reactions with energies: {thermo['num_reactions_with_energies']}")
    if thermo['reaction_energy_stats']['mean'] is not None:
        print(f"  Mean reaction energy: {thermo['reaction_energy_stats']['mean']:.2f} kcal/mol")
        print(f"  Energy range: {thermo['reaction_energy_stats']['min']:.2f} to {thermo['reaction_energy_stats']['max']:.2f} kcal/mol")
    print(f"  Exothermic reactions: {thermo['exothermic_reactions']}")
    print(f"  Endothermic reactions: {thermo['endothermic_reactions']}")
    
    # Find key intermediates
    if stats['num_molecules'] > len(seed_molecules):
        intermediates = analyzer.identify_key_intermediates(top_n=5)
        print(f"\nKey intermediates:")
        for i, intermediate in enumerate(intermediates, 1):
            mol = intermediate['molecule']
            print(f"  {i}. {mol.name} ({mol.smiles}) - Score: {intermediate['centrality_score']:.3f}")
    
    # Find dominant pathways
    seed_smiles = [mol.smiles for mol in seed_molecules]
    pathways = analyzer.find_dominant_pathways(
        source_molecules=seed_smiles,
        max_pathways=3
    )
    
    if pathways:
        print(f"\nDominant pathways:")
        for i, pathway in enumerate(pathways, 1):
            print(f"  {i}. {pathway['source']} → {pathway['target']}")
            print(f"     Steps: {pathway['num_steps']}, Total energy: {pathway['total_energy']:.2f} kcal/mol")
            print(f"     Path: {' → '.join(pathway['path'][:3])}{'...' if len(pathway['path']) > 3 else ''}")
    
    # Step 5: Apply pruning
    print(f"\n--- Step 5: Applying network pruning ---")
    
    pruner = NetworkPruner(config={
        'energy_threshold': 60.0,  # More restrictive
        'barrier_threshold': 100.0,
        'min_degree': 1
    })
    
    pruned_network = pruner.prune_network(
        network,
        strategies=['energy_filter', 'dead_end_removal', 'isolated_node_removal'],
        preserve_seeds=True
    )
    
    pruned_stats = pruned_network.get_statistics()
    print(f"Pruned network statistics:")
    print(f"  Molecules: {stats['num_molecules']} → {pruned_stats['num_molecules']}")
    print(f"  Reactions: {stats['num_reactions']} → {pruned_stats['num_reactions']}")
    
    # Step 6: Export results
    print(f"\n--- Step 6: Exporting results ---")
    
    # Export network data
    analyzer_pruned = NetworkAnalyzer(pruned_network)
    network_data = analyzer_pruned.export_network_data(format='json')
    
    # Save to file
    output_file = 'network_generation_results.json'
    with open(output_file, 'w') as f:
        json.dump(network_data, f, indent=2, default=str)
    
    print(f"Network data exported to {output_file}")
    
    # Generation statistics
    gen_stats = generator.get_generation_statistics()
    if gen_stats:
        print(f"\nGeneration statistics:")
        for stat in gen_stats:
            print(f"  Gen {stat['generation']}: {stat['input_molecules']} molecules → "
                  f"{stat['new_reactions']} reactions ({stat['time']:.2f}s)")
    
    # Step 7: Summary
    print(f"\n--- Summary ---")
    print(f"✓ Successfully generated reaction network from {len(seed_molecules)} seed molecules")
    print(f"✓ Final network: {pruned_stats['num_molecules']} molecules, {pruned_stats['num_reactions']} reactions")
    print(f"✓ Network spans {pruned_stats['max_generation']} generations")
    print(f"✓ Total computation time: {generation_time:.2f} seconds")
    
    if pathways:
        print(f"✓ Identified {len(pathways)} dominant reaction pathways")
    
    if intermediates:
        print(f"✓ Found {len(intermediates)} key intermediate species")
    
    print(f"✓ Results exported to {output_file}")
    
    print(f"\n=== Network Generation Demo Completed Successfully! ===")
    print("RSNet is now capable of:")
    print("  • Automated reaction network generation")
    print("  • Quantum chemical screening with xTB")
    print("  • Network topology analysis")
    print("  • Dominant pathway identification")
    print("  • Intelligent network pruning")
    print("  • Comprehensive data export")


def quick_test():
    """Quick test with minimal settings."""
    print("=== Quick Network Generation Test ===")
    
    # Single molecule, minimal settings
    seed_molecules = [Molecule.from_smiles("CC", name="ethane")]
    environment = Environment(temperature=500.0)
    
    # Minimal configuration
    config = {
        'max_generations': 2,
        'max_species': 20,
        'energy_cutoff': 100.0
    }
    
    generator = NetworkGenerator(
        operators=[BondBreakingOperator(config={'min_fragment_size': 1})],
        screener=ReactionScreener(optimize_geometries=False, max_workers=1),
        config=config
    )
    
    print("Generating network...")
    network = generator.generate_network(seed_molecules, environment, max_time=60)
    
    stats = network.get_statistics()
    print(f"Generated network: {stats['num_molecules']} molecules, {stats['num_reactions']} reactions")
    
    return network


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "quick":
        quick_test()
    else:
        main()
