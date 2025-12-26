"""
Operator Activation Analysis
=============================
Understand why certain operators are/aren't active
"""

import sys
import os
sys.path.append(os.path.dirname(__file__) + '/..')

from rsnet.operators.registry import OperatorRegistry
from rsnet.core.environment import Environment
from rsnet.core.molecule import Molecule
from rdkit import Chem

def analyze_operator_activation():
    """Check which operators can be activated under SEI conditions."""
    
    print("="*70)
    print("OPERATOR ACTIVATION ANALYSIS")
    print("="*70)
    print()
    
    # Setup SEI environment (anode)
    env = Environment(
        temperature=298.0,
        voltage=0.1,  # vs Li/Li+
        electrode_type="anode",
        solvent="mixture"
    )
    
    print("Environment:")
    print(f"  Temperature: {env.temperature}K")
    print(f"  Voltage: {env.voltage}V (anode)")
    print(f"  Electrode: {env.electrode_type}")
    print()
    
    # Get active drives
    drives = env.get_active_drives()
    print("Active Driving Forces:")
    for drive, active in drives.items():
        if active:
            print(f"  ✓ {drive}")
    print()
    
    # Initialize registry
    registry = OperatorRegistry()
    
    # Test molecules
    test_molecules = [
        Molecule(Chem.MolFromSmiles("O=C1OCCO1"), name="EC"),
        Molecule(Chem.MolFromSmiles("COC(=O)OC"), name="DMC"),
        Molecule(Chem.MolFromSmiles("[Li+]"), name="Li+"),
        Molecule(Chem.MolFromSmiles("F[P-](F)(F)(F)(F)F"), name="PF6-"),
    ]
    
    print("Test Molecules: EC, DMC, Li+, PF6-")
    print()
    
    # Get all registered operators
    all_operators = registry.get_all_operators()
    
    print(f"Registered Operators: {len(all_operators)}")
    print("-"*70)
    
    activated = []
    not_activated = []
    
    for name, op in all_operators.items():
        # Test if can apply
        try:
            can_apply = op.can_apply(test_molecules, env)
            
            if can_apply:
                activated.append(name)
                print(f"✓ {name:<25} CAN APPLY")
            else:
                not_activated.append(name)
                print(f"✗ {name:<25} CANNOT APPLY")
                
        except Exception as e:
            not_activated.append(name)
            print(f"⚠ {name:<25} ERROR: {str(e)[:40]}")
    
    print()
    print("="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Activated: {len(activated)}/{len(all_operators)}")
    print(f"Not Activated: {len(not_activated)}/{len(all_operators)}")
    print()
    
    # Focus on RingOpening and Polymerization
    print("="*70)
    print("CRITICAL OPERATORS")
    print("="*70)
    
    critical = ['ring_opening', 'polymerization', 'bond_breaking', 'radical_reaction']
    
    for op_name in critical:
        if op_name in all_operators:
            op = all_operators[op_name]
            status = "✓ ACTIVE" if op_name in activated else "✗ NOT ACTIVE"
            print(f"{op_name:<25} {status}")
            
            # Get more details
            if hasattr(op, 'can_apply'):
                print(f"  Description: {op.description}")
                
                # Check with single EC molecule for ring opening
                if op_name == 'ring_opening':
                    ec = [Molecule(Chem.MolFromSmiles("O=C1OCCO1"), name="EC")]
                    try:
                        can_apply_ec = op.can_apply(ec, env)
                        print(f"  Can apply to EC? {can_apply_ec}")
                    except Exception as e:
                        print(f"  Error with EC: {e}")
        print()
    
    return {
        'activated': activated,
        'not_activated': not_activated,
        'environment': env,
        'drives': drives
    }

if __name__ == "__main__":
    results = analyze_operator_activation()
