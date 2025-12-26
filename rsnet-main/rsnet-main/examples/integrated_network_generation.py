#!/usr/bin/env python3
"""
Integrated Network Generation Example

This example demonstrates the complete integrated RSNet system with:
- Enhanced environment modeling
- Intelligent operator selection
- Driving force evaluation
- Structure-based filtering
- Advanced configuration options

Run with: python integrated_network_generation.py
"""

import sys
import os
import time
from pathlib import Path

# Add rsnet to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.network.generator import NetworkGenerator
from rsnet.network.config import (NetworkGenerationConfig, NetworkGenerationPresets,
                                   get_preset_config, OperatorSelectionStrategy)
from rsnet.operators.registry import OPERATOR_REGISTRY
from rsnet.features.driving_forces import get_driving_forces
from rsnet.features.structure_tags import get_structure_tags


def demonstrate_basic_integration():
    """Demonstrate basic integrated network generation."""
    print("=== Basic Integrated Network Generation ===")
    
    # Create molecules
    ethanol = Molecule.from_smiles('CCO', name='ethanol')
    print(f"Starting molecule: {ethanol.name} ({ethanol.smiles})")
    
    # Enhanced environment
    env = Environment(
        temperature=500.0,
        electrode_type='cathode',
        voltage=4.2,
        li_activity=0.5,
        interface_type='SEI'
    )
    print(f"Environment: T={env.temperature}K, electrode={env.electrode_type}, V={env.voltage}V")
    
    # Show driving forces
    forces = get_driving_forces([ethanol], env)
    strong_forces = {k: v for k, v in forces.items() if v > 0.3}
    print(f"Strong driving forces: {strong_forces}")
    
    # Show structure tags
    tags = get_structure_tags([ethanol])
    print(f"Molecular features: {len(tags.get('functional_groups', []))} functional groups, "
          f"MW={tags.get('molecular_weight', 0):.1f}")
    
    # Configure network generation
    config = NetworkGenerationConfig(
        max_generations=3,
        max_species=25,
        energy_cutoff=55.0,
        driving_force_threshold=0.25,
        max_operators_per_generation=4,
        use_structure_based_filtering=True
    )
    
    # Generate network
    generator = NetworkGenerator(config=config)
    print(f"Using intelligent selection: {generator.use_intelligent_selection}")
    
    start_time = time.time()
    network = generator.generate_network([ethanol], env, max_time=30.0)
    generation_time = time.time() - start_time
    
    # Show results
    stats = network.get_statistics()
    print(f"\nResults after {generation_time:.2f}s:")
    print(f"  Molecules: {stats['num_molecules']}")
    print(f"  Reactions: {stats['num_reactions']}")
    print(f"  Generations: {stats['max_generation'] + 1}")
    print(f"  Molecules by generation: {stats['molecules_by_generation']}")
    
    # Show generation statistics
    gen_stats = generator.get_generation_statistics()
    for stat in gen_stats:
        print(f"  Gen {stat['generation']}: {stat['new_reactions']} reactions, "
              f"{stat['total_molecules']} molecules, {stat['time']:.2f}s")
    
    return network


def demonstrate_preset_configurations():
    """Demonstrate different preset configurations."""
    print("\n=== Preset Configuration Comparison ===")
    
    # Test molecule
    acetaldehyde = Molecule.from_smiles('CC=O', name='acetaldehyde')
    env = Environment(temperature=450.0, solvent='water')
    
    presets = [
        ('battery_chemistry', 'Battery Chemistry'),
        ('organic_synthesis', 'Organic Synthesis'),
        ('fast_screening', 'Fast Screening')
    ]
    
    results = {}
    
    for preset_name, display_name in presets:
        print(f"\n--- {display_name} ---")
        
        config = get_preset_config(preset_name)
        generator = NetworkGenerator(config=config)
        
        start_time = time.time()
        network = generator.generate_network([acetaldehyde], env, max_time=15.0)
        elapsed_time = time.time() - start_time
        
        stats = network.get_statistics()
        results[preset_name] = {
            'molecules': stats['num_molecules'],
            'reactions': stats['num_reactions'],
            'time': elapsed_time,
            'generations': stats['max_generation'] + 1
        }
        
        print(f"  Config: max_gen={config.max_generations}, "
              f"max_species={config.max_species}, "
              f"strategy={config.generation_strategy.value}")
        print(f"  Results: {stats['num_molecules']} molecules, "
              f"{stats['num_reactions']} reactions in {elapsed_time:.2f}s")
    
    # Compare results
    print("\n--- Comparison ---")
    for preset_name, result in results.items():
        print(f"{preset_name:20}: {result['molecules']:3d} mol, "
              f"{result['reactions']:3d} rxn, {result['time']:5.2f}s")


def demonstrate_operator_intelligence():
    """Demonstrate intelligent operator selection."""
    print("\n=== Intelligent Operator Selection ===")
    
    # Different molecules with different characteristics
    molecules = [
        Molecule.from_smiles('CCO', name='ethanol'),
        Molecule.from_smiles('C=C', name='ethylene'),
        Molecule.from_smiles('c1ccccc1', name='benzene'),
        Molecule.from_smiles('CC(=O)O', name='acetic_acid')
    ]
    
    # Different environments
    environments = [
        ('High Temperature', Environment(temperature=700.0)),
        ('Electrochemical', Environment(temperature=350.0, electrode_type='cathode', voltage=4.0)),
        ('Mild Conditions', Environment(temperature=298.15, solvent='water'))
    ]
    
    for env_name, env in environments:
        print(f"\n--- {env_name} Environment ---")
        
        for mol in molecules:
            # Get recommended operators
            operators = OPERATOR_REGISTRY.get_recommended_operators([mol], env, max_operators=3)
            operator_names = [op.name for op in operators]
            
            # Get driving forces
            forces = get_driving_forces([mol], env)
            strong_forces = [k for k, v in forces.items() if v > 0.4]
            
            print(f"  {mol.name:12}: operators={operator_names}, "
                  f"strong_forces={strong_forces}")


def demonstrate_advanced_features():
    """Demonstrate advanced features and customization."""
    print("\n=== Advanced Features ===")
    
    # Complex molecule
    mol = Molecule.from_smiles('CC(C)(C)C(=O)O', name='tert-butyl_acetate')
    
    # Complex environment
    env = Environment(
        temperature=550.0,
        pressure=2.0,
        electrode_type='cathode',
        voltage=4.5,
        li_activity=0.8,
        interface_type='SEI'
    )
    
    # Advanced configuration
    config = NetworkGenerationConfig(
        max_generations=4,
        max_species=40,
        energy_cutoff=50.0,
        activation_energy_cutoff=75.0,
        
        # Intelligent selection
        operator_selection_strategy=OperatorSelectionStrategy.ADAPTIVE,
        max_operators_per_generation=3,
        driving_force_threshold=0.2,
        
        # Structure filtering
        use_structure_based_filtering=True,
        max_molecular_weight=300.0,
        filter_trivial_reactions=True,
        min_structural_change=0.15,
        
        # Performance
        parallel_screening=True,
        max_screening_workers=2,
        
        # Advanced options
        prioritize_novel_structures=True,
        bias_towards_stable_products=False,
        allow_unstable_intermediates=True,
        
        # Debugging
        verbose_logging=True
    )
    
    print(f"Molecule: {mol.name}")
    print(f"Environment: T={env.temperature}K, P={env.pressure}atm, "
          f"electrode={env.electrode_type}, V={env.voltage}V")
    
    # Show molecular analysis
    tags = get_structure_tags([mol])
    forces = get_driving_forces([mol], env)
    
    print(f"Molecular weight: {tags.get('molecular_weight', 0):.1f} g/mol")
    print(f"Functional groups: {tags.get('functional_groups', [])}")
    print(f"Strong driving forces: {[k for k, v in forces.items() if v > 0.4]}")
    
    # Generate network
    generator = NetworkGenerator(config=config)
    
    print(f"\nConfiguration:")
    print(f"  Strategy: {config.generation_strategy.value}")
    print(f"  Operator selection: {config.operator_selection_strategy.value}")
    print(f"  Max operators/gen: {config.max_operators_per_generation}")
    print(f"  Driving force threshold: {config.driving_force_threshold}")
    
    start_time = time.time()
    network = generator.generate_network([mol], env, max_time=45.0)
    generation_time = time.time() - start_time
    
    # Detailed results
    stats = network.get_statistics()
    gen_stats = generator.get_generation_statistics()
    
    print(f"\nDetailed Results:")
    print(f"  Total time: {generation_time:.2f}s")
    print(f"  Final network: {stats['num_molecules']} molecules, {stats['num_reactions']} reactions")
    print(f"  Generations completed: {len(gen_stats)}")
    
    for i, stat in enumerate(gen_stats):
        print(f"    Gen {i}: {stat['input_molecules']} input → "
              f"{stat['new_reactions']} reactions → "
              f"{stat['total_molecules']} total molecules ({stat['time']:.2f}s)")
    
    # Show some example reactions
    if network.reactions:
        print(f"\nExample reactions:")
        for i, (reaction_id, reaction) in enumerate(list(network.reactions.items())[:3]):
            reactants = [r.smiles for r in reaction.reactants]
            products = [p.smiles for p in reaction.products]
            energy = getattr(reaction, 'reaction_energy', 'N/A')
            driving_score = getattr(reaction, 'driving_force_score', 'N/A')
            
            print(f"  {i+1}. {' + '.join(reactants)} → {' + '.join(products)}")
            print(f"     Energy: {energy} kcal/mol, Driving force: {driving_score}")


def main():
    """Run all demonstrations."""
    print("RSNet Integrated Network Generation Demo")
    print("=" * 50)
    
    try:
        # Basic integration
        network = demonstrate_basic_integration()
        
        # Preset configurations
        demonstrate_preset_configurations()
        
        # Operator intelligence
        demonstrate_operator_intelligence()
        
        # Advanced features
        demonstrate_advanced_features()
        
        print("\n" + "=" * 50)
        print("Demo completed successfully!")
        print("The integrated RSNet system demonstrates:")
        print("- Intelligent operator selection based on molecular features")
        print("- Driving force evaluation for reaction feasibility")
        print("- Structure-based filtering and grouping")
        print("- Flexible configuration with presets")
        print("- Enhanced reaction screening and ranking")
        
    except Exception as e:
        print(f"\nError during demonstration: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == '__main__':
    exit(main())
