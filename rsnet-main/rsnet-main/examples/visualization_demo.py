#!/usr/bin/env python3
"""
RSNet Visualization Demo

This example demonstrates the visualization capabilities of RSNet,
including network graphs, reaction pathways, energy distributions,
and molecular structures.

Run with: python visualization_demo.py
"""

import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt

# Add rsnet to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.network.generator import NetworkGenerator
from rsnet.network.config import NetworkGenerationConfig
from rsnet.utils.visualization import (
    plot_network_graph,
    plot_reaction_pathway,
    plot_energy_distribution,
    plot_molecule_structure,
    plot_generation_statistics
)
from rsnet.utils.io import create_analysis_report


def demo_network_visualization():
    """Demonstrate network graph visualization."""
    print("=== Network Graph Visualization ===")
    
    # Create a small network for visualization
    ethanol = Molecule.from_smiles('CCO', name='ethanol')
    env = Environment(temperature=500.0, electrode_type='cathode', voltage=4.0)
    
    config = NetworkGenerationConfig(
        max_generations=2,
        max_species=15,
        energy_cutoff=60.0
    )
    
    generator = NetworkGenerator(config=config)
    network = generator.generate_network([ethanol], env, max_time=20.0)
    
    stats = network.get_statistics()
    print(f"Generated network: {stats['num_molecules']} molecules, {stats['num_reactions']} reactions")
    
    # Create different network visualizations
    layouts = ['spring', 'circular', 'kamada_kawai']
    
    fig, axes = plt.subplots(1, len(layouts), figsize=(15, 5))
    
    for i, layout in enumerate(layouts):
        ax = axes[i] if len(layouts) > 1 else axes
        
        # Plot network with different layouts
        plot_network_graph(
            network, 
            layout=layout,
            figsize=(5, 5),
            show_labels=True
        )
        
        plt.sca(ax)
        plt.title(f'{layout.title()} Layout')
    
    plt.tight_layout()
    plt.savefig('network_layouts.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("Network visualizations saved as 'network_layouts.png'")
    return network


def demo_energy_visualization(network):
    """Demonstrate energy distribution visualization."""
    print("\n=== Energy Distribution Visualization ===")
    
    # Plot energy distribution
    fig = plot_energy_distribution(network, figsize=(12, 6))
    plt.savefig('energy_distribution.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("Energy distribution saved as 'energy_distribution.png'")


def demo_generation_statistics(network):
    """Demonstrate generation statistics visualization."""
    print("\n=== Generation Statistics Visualization ===")
    
    # Plot generation statistics
    fig = plot_generation_statistics(network, figsize=(12, 5))
    plt.savefig('generation_stats.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("Generation statistics saved as 'generation_stats.png'")


def demo_molecule_visualization():
    """Demonstrate molecular structure visualization."""
    print("\n=== Molecular Structure Visualization ===")
    
    # Create some interesting molecules
    molecules = [
        Molecule.from_smiles('CCO', name='Ethanol'),
        Molecule.from_smiles('c1ccccc1', name='Benzene'),
        Molecule.from_smiles('CC(=O)O', name='Acetic Acid'),
        Molecule.from_smiles('C1CC1', name='Cyclopropane')
    ]
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    axes = axes.flatten()
    
    for i, mol in enumerate(molecules):
        plt.sca(axes[i])
        
        try:
            plot_molecule_structure(mol, figsize=(6, 6))
        except Exception as e:
            print(f"Could not plot {mol.name}: {e}")
            axes[i].text(0.5, 0.5, f'{mol.name}\n{mol.smiles}\nMW: {mol.molecular_weight:.1f}',
                        ha='center', va='center', transform=axes[i].transAxes)
            axes[i].axis('off')
    
    plt.tight_layout()
    plt.savefig('molecular_structures.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("Molecular structures saved as 'molecular_structures.png'")


def demo_pathway_visualization(network):
    """Demonstrate reaction pathway visualization."""
    print("\n=== Reaction Pathway Visualization ===")
    
    # Get molecules from network
    molecules = list(network.molecules.keys())
    
    if len(molecules) >= 2:
        start_mol = molecules[0]  # Starting molecule
        
        # Find a target molecule (preferably from a later generation)
        target_mol = None
        for mol_smiles in molecules[1:]:
            mol = network.molecules[mol_smiles]
            if getattr(mol, 'generation', 0) > 0:
                target_mol = mol_smiles
                break
        
        if target_mol:
            print(f"Finding pathways from {start_mol} to {target_mol}")
            
            try:
                fig = plot_reaction_pathway(
                    network, 
                    start_mol, 
                    target_mol,
                    max_pathways=3,
                    figsize=(14, 8)
                )
                plt.savefig('reaction_pathways.png', dpi=300, bbox_inches='tight')
                plt.show()
                
                print("Reaction pathways saved as 'reaction_pathways.png'")
            except Exception as e:
                print(f"Could not create pathway visualization: {e}")
        else:
            print("No suitable target molecule found for pathway analysis")
    else:
        print("Network too small for pathway analysis")


def demo_comprehensive_analysis():
    """Demonstrate comprehensive analysis report generation."""
    print("\n=== Comprehensive Analysis Report ===")
    
    # Generate a larger network for analysis
    molecules = [
        Molecule.from_smiles('CCO', name='ethanol'),
        Molecule.from_smiles('C=C', name='ethylene')
    ]
    
    env = Environment(
        temperature=450.0,
        electrode_type='cathode',
        voltage=4.2,
        li_activity=0.3
    )
    
    config = NetworkGenerationConfig(
        max_generations=3,
        max_species=25,
        energy_cutoff=55.0,
        driving_force_threshold=0.25
    )
    
    generator = NetworkGenerator(config=config)
    network = generator.generate_network(molecules, env, max_time=30.0)
    
    stats = network.get_statistics()
    print(f"Analysis network: {stats['num_molecules']} molecules, {stats['num_reactions']} reactions")
    
    # Create comprehensive analysis report
    output_dir = Path('analysis_report')
    create_analysis_report(network, output_dir, include_plots=True)
    
    print(f"Comprehensive analysis report created in '{output_dir}'")
    print("Report includes:")
    print("  - Network data (JSON, GraphML)")
    print("  - Reaction and molecule tables (CSV)")
    print("  - Network statistics and topology analysis")
    print("  - Visualization plots (PNG)")
    print("  - Summary report (TXT)")


def demo_custom_visualization():
    """Demonstrate custom visualization techniques."""
    print("\n=== Custom Visualization Techniques ===")
    
    # Create a simple network
    mol = Molecule.from_smiles('CC', name='ethane')
    env = Environment(temperature=600.0)
    
    config = NetworkGenerationConfig(max_generations=2, max_species=10)
    generator = NetworkGenerator(config=config)
    network = generator.generate_network([mol], env, max_time=15.0)
    
    # Custom network plot with enhanced styling
    fig, ax = plt.subplots(figsize=(10, 8))
    
    plot_network_graph(
        network,
        layout='spring',
        figsize=(10, 8),
        node_size_scale=2.0,
        edge_width_scale=1.5,
        show_labels=True
    )
    
    plt.title('Enhanced Network Visualization\nEthane Reaction Network at 600K', 
              fontsize=16, fontweight='bold')
    plt.savefig('custom_network.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("Custom network visualization saved as 'custom_network.png'")


def main():
    """Run all visualization demonstrations."""
    print("RSNet Visualization Demo")
    print("=" * 50)
    
    try:
        # Create output directory
        output_dir = Path('visualization_output')
        output_dir.mkdir(exist_ok=True)
        os.chdir(output_dir)
        
        # Run demonstrations
        network = demo_network_visualization()
        
        if len(network.reactions) > 0:
            demo_energy_visualization(network)
        else:
            print("Skipping energy visualization (no reactions with energy data)")
        
        demo_generation_statistics(network)
        demo_molecule_visualization()
        
        if len(network.molecules) > 1:
            demo_pathway_visualization(network)
        else:
            print("Skipping pathway visualization (network too small)")
        
        demo_comprehensive_analysis()
        demo_custom_visualization()
        
        print("\n" + "=" * 50)
        print("Visualization demo completed successfully!")
        print(f"All outputs saved in '{output_dir.absolute()}'")
        print("\nGenerated files:")
        for file in output_dir.glob('*.png'):
            print(f"  - {file.name}")
        
        if (output_dir / 'analysis_report').exists():
            print("  - analysis_report/ (comprehensive analysis)")
        
        print("\nVisualization capabilities demonstrated:")
        print("  ✓ Network graph layouts (spring, circular, kamada-kawai)")
        print("  ✓ Energy distribution analysis")
        print("  ✓ Generation statistics")
        print("  ✓ Molecular structure plots")
        print("  ✓ Reaction pathway visualization")
        print("  ✓ Comprehensive analysis reports")
        print("  ✓ Custom visualization styling")
        
    except Exception as e:
        print(f"\nError during visualization demo: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == '__main__':
    exit(main())
