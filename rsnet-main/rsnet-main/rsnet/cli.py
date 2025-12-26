#!/usr/bin/env python3
"""
Command-line interface for RSNet.

This module provides a simple CLI for running reaction network generation
and analysis tasks.
"""

import argparse
import sys
from pathlib import Path
from typing import List, Optional
import time

import sys
from pathlib import Path

# Add rsnet to path for standalone execution
if __name__ == '__main__':
    sys.path.insert(0, str(Path(__file__).parent.parent))

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.network.generator import NetworkGenerator
from rsnet.network.analyzer import NetworkAnalyzer
from rsnet.utils.config import load_config_with_overrides, get_preset_config, print_config_summary
from rsnet.utils.io import (load_molecules_from_smiles, load_molecules_from_smiles_list,
                           save_network, create_analysis_report,
                           export_reaction_table, export_molecule_table)
from rsnet.utils.visualization import plot_network_graph, plot_energy_distribution


def create_parser() -> argparse.ArgumentParser:
    """Create command-line argument parser."""
    parser = argparse.ArgumentParser(
        description='RSNet - Reaction Network Generation Engine',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate network from single SMILES
  rsnet generate -s "CCO" -o output/

  # Generate network from multiple SMILES
  rsnet generate -m "CCO" "CC" "C=C" -o output/

  # Use molecules from file
  rsnet generate -f molecules.csv -o output/

  # Use preset configuration
  rsnet generate -s "CCO" --preset battery_chemistry -o output/

  # Analyze existing network
  rsnet analyze -n network.json -o analysis/

  # Create configuration template
  rsnet config --template config.yaml --preset organic_synthesis
        """
    )
    
    parser.add_argument('--version', action='version', version='RSNet 0.1.0')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Generate command
    gen_parser = subparsers.add_parser('generate', help='Generate reaction network')
    gen_parser.add_argument('-s', '--smiles', type=str, help='Input SMILES string')
    gen_parser.add_argument('-m', '--multi-smiles', type=str, nargs='+', help='Multiple SMILES strings')
    gen_parser.add_argument('-f', '--file', type=str, help='Input file with molecules')
    gen_parser.add_argument('-o', '--output', type=str, required=True, help='Output directory')
    gen_parser.add_argument('-c', '--config', type=str, help='Configuration file')
    gen_parser.add_argument('--preset', type=str, help='Preset configuration')
    gen_parser.add_argument('--max-time', type=float, default=300.0, help='Maximum time (seconds)')
    gen_parser.add_argument('--temperature', type=float, help='Temperature (K)')
    gen_parser.add_argument('--no-plots', action='store_true', help='Skip plot generation')
    
    # Analyze command
    analyze_parser = subparsers.add_parser('analyze', help='Analyze reaction network')
    analyze_parser.add_argument('-n', '--network', type=str, required=True, help='Network file')
    analyze_parser.add_argument('-o', '--output', type=str, required=True, help='Output directory')
    analyze_parser.add_argument('--pathways', type=str, nargs=2, metavar=('START', 'END'),
                               help='Find pathways between two SMILES')
    analyze_parser.add_argument('--no-plots', action='store_true', help='Skip plot generation')
    
    # Config command
    config_parser = subparsers.add_parser('config', help='Configuration management')
    config_parser.add_argument('--template', type=str, help='Create configuration template')
    config_parser.add_argument('--preset', type=str, help='Use preset configuration')
    config_parser.add_argument('--validate', type=str, help='Validate configuration file')
    config_parser.add_argument('--show', type=str, help='Show configuration summary')
    
    return parser


def generate_network(args) -> int:
    """Generate reaction network."""
    try:
        # Load configuration
        if args.config:
            config = load_config_with_overrides(args.config)
        elif args.preset:
            config = get_preset_config(args.preset)
        else:
            config = load_config_with_overrides()
        
        if args.verbose:
            print_config_summary(config)
        
        # Load molecules
        if args.smiles:
            molecules = [Molecule.from_smiles(args.smiles, name='input_mol')]
        elif args.multi_smiles:
            molecules = load_molecules_from_smiles_list(args.multi_smiles)
            print(f"Loaded {len(molecules)} molecules from multi-SMILES input")
        elif args.file:
            molecules = load_molecules_from_smiles(args.file)
        else:
            print("Error: Must specify either --smiles, --multi-smiles, or --file")
            return 1
        
        if not molecules:
            print("Error: No valid molecules found")
            return 1
        
        print(f"Loaded {len(molecules)} molecules")
        
        # Create environment
        env_config = config.environment
        if args.temperature:
            env_config['temperature'] = args.temperature
        
        env = Environment(**env_config)
        
        # Create generator
        generator = NetworkGenerator(config=config.network)
        
        print(f"Starting network generation...")
        print(f"  Max generations: {config.network.max_generations}")
        print(f"  Max species: {config.network.max_species}")
        print(f"  Energy cutoff: {config.network.energy_cutoff} kcal/mol")
        
        # Generate network
        start_time = time.time()
        network = generator.generate_network(molecules, env, max_time=args.max_time)
        generation_time = time.time() - start_time
        
        # Show results
        stats = network.get_statistics()
        print(f"\nGeneration completed in {generation_time:.2f}s")
        print(f"  Final network: {stats['num_molecules']} molecules, {stats['num_reactions']} reactions")
        print(f"  Generations: {stats['max_generation'] + 1}")
        
        # Create output directory
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save network and analysis
        print(f"\nSaving results to {output_dir}")
        create_analysis_report(network, output_dir, include_plots=not args.no_plots)
        
        print("Generation completed successfully!")
        return 0
        
    except Exception as e:
        print(f"Error during generation: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


def analyze_network(args) -> int:
    """Analyze existing network."""
    try:
        from rsnet.utils.io import load_network
        
        # Load network
        network = load_network(args.network)
        print(f"Loaded network with {len(network.molecules)} molecules and {len(network.reactions)} reactions")
        
        # Create analyzer
        analyzer = NetworkAnalyzer(network)
        
        # Create output directory
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Basic analysis
        print("Performing network analysis...")
        create_analysis_report(network, output_dir, include_plots=not args.no_plots)
        
        # Pathway analysis if requested
        if args.pathways:
            start_smiles, end_smiles = args.pathways
            print(f"Finding pathways from {start_smiles} to {end_smiles}")
            
            pathways = analyzer.find_pathways(start_smiles, end_smiles, max_pathways=5)
            
            if pathways:
                print(f"Found {len(pathways)} pathways:")
                for i, pathway in enumerate(pathways):
                    print(f"  Pathway {i+1}: {len(pathway['molecules'])-1} steps")
                
                # Save pathway analysis
                import json
                with open(output_dir / 'pathways.json', 'w') as f:
                    json.dump(pathways, f, indent=2)
                
                # Plot pathways
                if not args.no_plots:
                    from rsnet.utils.visualization import plot_reaction_pathway
                    fig = plot_reaction_pathway(network, start_smiles, end_smiles)
                    fig.savefig(output_dir / 'pathways.png', dpi=300, bbox_inches='tight')
                    plt.close(fig)
            else:
                print("No pathways found")
        
        print("Analysis completed successfully!")
        return 0
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


def manage_config(args) -> int:
    """Manage configuration."""
    try:
        if args.template:
            from rsnet.utils.config import create_config_template
            create_config_template(args.template, preset=args.preset)
            print(f"Configuration template created: {args.template}")
            
        elif args.validate:
            from rsnet.utils.config import load_config, validate_config
            config = load_config(args.validate)
            issues = validate_config(config)
            
            if issues:
                print("Configuration validation failed:")
                for issue in issues:
                    print(f"  - {issue}")
                return 1
            else:
                print("Configuration is valid")
                
        elif args.show:
            from rsnet.utils.config import load_config
            config = load_config(args.show)
            print_config_summary(config)
            
        else:
            print("Error: Must specify --template, --validate, or --show")
            return 1
        
        return 0
        
    except Exception as e:
        print(f"Error managing configuration: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


def main() -> int:
    """Main CLI entry point."""
    parser = create_parser()
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 1
    
    # Set up logging if verbose
    if args.verbose:
        import logging
        logging.basicConfig(level=logging.INFO)
    
    # Dispatch to appropriate function
    if args.command == 'generate':
        return generate_network(args)
    elif args.command == 'analyze':
        return analyze_network(args)
    elif args.command == 'config':
        return manage_config(args)
    else:
        print(f"Unknown command: {args.command}")
        return 1


if __name__ == '__main__':
    sys.exit(main())
