#!/usr/bin/env python3
"""
Enhanced Features Demo for RSNet

This script demonstrates the enhanced functionality of RSNet including:
1. Physical driving force system
2. Structure annotation system  
3. Expanded operator types
4. Operator registry and selection

Author: RSNet Development Team
Date: 2025-12-22
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.features.structure_tags import get_structure_tags
from rsnet.features.driving_forces import get_driving_forces
from rsnet.operators.registry import OPERATOR_REGISTRY


def demo_enhanced_environment():
    """Demonstrate enhanced environment with physical driving forces."""
    print("=" * 60)
    print("1. Enhanced Environment and Driving Forces")
    print("=" * 60)
    
    # Create different environments
    environments = [
        {
            'name': 'High Temperature Thermal',
            'env': Environment(temperature=800.0, pressure=1.0)
        },
        {
            'name': 'Electrochemical (Cathode)',
            'env': Environment(
                temperature=298.15,
                electrode_type='cathode',
                voltage=4.2,
                li_activity=0.5,
                interface_type='SEI'
            )
        },
        {
            'name': 'Electrochemical (Anode)',
            'env': Environment(
                temperature=298.15,
                electrode_type='anode',
                voltage=0.1,
                li_activity=1.0,
                interface_type='SEI'
            )
        }
    ]
    
    for env_info in environments:
        print(f"\n{env_info['name']} Environment:")
        env = env_info['env']
        drives = env.get_active_drives()
        
        print(f"  Temperature: {env.temperature} K")
        if hasattr(env, 'electrode_type') and env.electrode_type:
            print(f"  Electrode: {env.electrode_type}, Voltage: {env.voltage} V")
            print(f"  Li+ activity: {env.li_activity}")
            print(f"  Interface: {env.interface_type}")
        
        print("  Active driving forces:")
        for drive, active in drives.items():
            if active:
                print(f"    ✓ {drive}")


def demo_structure_annotation():
    """Demonstrate structure annotation system."""
    print("\n" + "=" * 60)
    print("2. Structure Annotation System")
    print("=" * 60)
    
    # Test molecules with different structural features
    test_molecules = [
        ('CCO', 'Ethanol'),
        ('c1ccccc1', 'Benzene'),
        ('CC(=O)C', 'Acetone'),
        ('C1CC1', 'Cyclopropane'),
        ('C=C', 'Ethylene'),
        ('CC(C)(C)C', 'tert-Butane'),
        ('[N+](=O)[O-]', 'Nitro group'),
        ('C#N', 'Acetonitrile')
    ]
    
    for smiles, name in test_molecules:
        print(f"\n{name} ({smiles}):")
        mol = Molecule.from_smiles(smiles)
        tags = get_structure_tags([mol])
        
        print("  Structural features:")
        for feature, present in tags.items():
            if present:
                print(f"    ✓ {feature}")


def demo_driving_force_evaluation():
    """Demonstrate driving force strength evaluation."""
    print("\n" + "=" * 60)
    print("3. Driving Force Strength Evaluation")
    print("=" * 60)
    
    # Test molecule
    mol = Molecule.from_smiles('CCO')  # Ethanol
    
    # Different environments
    environments = [
        ('Thermal (500K)', Environment(temperature=500.0)),
        ('High Voltage Cathode', Environment(
            temperature=298.15,
            electrode_type='cathode',
            voltage=4.5,
            li_activity=0.8
        )),
        ('Low Voltage Anode', Environment(
            temperature=298.15,
            electrode_type='anode',
            voltage=0.05,
            li_activity=1.0
        ))
    ]
    
    for env_name, env in environments:
        print(f"\n{env_name}:")
        drives = get_driving_forces([mol], env)
        
        print("  Driving force strengths (0.0-1.0):")
        for drive, strength in drives.items():
            if strength > 0.0:
                print(f"    {drive}: {strength:.3f}")


def demo_operator_registry():
    """Demonstrate operator registry and selection."""
    print("\n" + "=" * 60)
    print("4. Operator Registry and Selection")
    print("=" * 60)
    
    # Test molecules
    molecules = [
        ('CCO', 'Ethanol'),
        ('C=C', 'Ethylene'),
        ('C=O', 'Formaldehyde'),
        ('c1ccccc1', 'Benzene')
    ]
    
    # Test environment
    env = Environment(
        temperature=500.0,
        electrode_type='cathode',
        voltage=4.2,
        li_activity=0.5
    )
    
    print(f"Environment: {env.temperature}K, cathode at {env.voltage}V")
    print(f"Available operators: {list(OPERATOR_REGISTRY._operators.keys())}")
    
    for smiles, name in molecules:
        print(f"\n{name} ({smiles}):")
        mol = Molecule.from_smiles(smiles)
        
        # Get active operators
        active_ops = OPERATOR_REGISTRY.get_active_operators([mol], env)
        print(f"  Active operators: {[op.name for op in active_ops]}")
        
        # Get recommended operators
        recommended = OPERATOR_REGISTRY.get_recommended_operators([mol], env, max_operators=3)
        print(f"  Recommended (top 3): {[op.name for op in recommended]}")
        
        # Show operator priorities
        print("  Operator priorities:")
        for op in recommended:
            priority = OPERATOR_REGISTRY.get_operator_priority(op, [mol], env)
            print(f"    {op.name}: {priority:.3f}")


def demo_reaction_generation():
    """Demonstrate reaction generation with new operators."""
    print("\n" + "=" * 60)
    print("5. Reaction Generation with New Operators")
    print("=" * 60)
    
    # Test molecules
    test_cases = [
        ('CCO', 'Ethanol', 'oxidation'),
        ('C=O', 'Formaldehyde', 'reduction'),
        ('C=C', 'Ethylene', 'addition'),
        ('CCCC', 'Butane', 'bond_breaking')
    ]
    
    for smiles, name, expected_reaction in test_cases:
        print(f"\n{name} ({smiles}) - Expected: {expected_reaction}")
        mol = Molecule.from_smiles(smiles)
        
        # Create appropriate environment
        if expected_reaction == 'oxidation':
            env = Environment(electrode_type='cathode', voltage=4.5)
        elif expected_reaction == 'reduction':
            env = Environment(electrode_type='anode', voltage=0.1)
        else:
            env = Environment(temperature=500.0)
        
        # Get recommended operators
        recommended = OPERATOR_REGISTRY.get_recommended_operators([mol], env, max_operators=2)
        
        for op in recommended:
            print(f"  Testing {op.name} operator:")
            
            try:
                # Find reaction sites
                sites = op.find_reaction_sites([mol])
                print(f"    Reaction sites found: {len(sites)}")
                
                # Generate reactions
                reactions = op.apply([mol], env)
                print(f"    Reactions generated: {len(reactions)}")
                
                for i, reaction in enumerate(reactions[:2]):  # Show first 2 reactions
                    print(f"      Reaction {i+1}: {reaction.name}")
                    reactant_smiles = [r.smiles for r in reaction.reactants]
                    product_smiles = [p.smiles for p in reaction.products]
                    print(f"        {' + '.join(reactant_smiles)} → {' + '.join(product_smiles)}")
                    
            except Exception as e:
                print(f"    Error: {e}")


def main():
    """Run all demonstrations."""
    print("RSNet Enhanced Features Demonstration")
    print("====================================")
    
    try:
        demo_enhanced_environment()
        demo_structure_annotation()
        demo_driving_force_evaluation()
        demo_operator_registry()
        demo_reaction_generation()
        
        print("\n" + "=" * 60)
        print("Demo completed successfully!")
        print("=" * 60)
        
    except Exception as e:
        print(f"\nError during demonstration: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
