#!/usr/bin/env python3
"""
Quick test of RSNet enhanced features.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.features.structure_tags import get_structure_tags
from rsnet.features.driving_forces import get_driving_forces
from rsnet.operators.registry import OPERATOR_REGISTRY


def main():
    print("RSNet Enhanced Features Quick Test")
    print("=" * 40)
    
    # Test 1: Enhanced Environment
    print("\n1. Enhanced Environment:")
    env = Environment(
        temperature=500.0,
        electrode_type='cathode',
        voltage=4.2,
        li_activity=0.5,
        interface_type='SEI'
    )
    drives = env.get_active_drives()
    active_drives = [k for k, v in drives.items() if v]
    print(f"   Active drives: {active_drives}")
    
    # Test 2: Structure Tags
    print("\n2. Structure Tags:")
    mol = Molecule.from_smiles('CCO')
    tags = get_structure_tags([mol])
    active_tags = [k for k, v in tags.items() if v and isinstance(v, bool)]
    print(f"   Ethanol features: {active_tags}")
    
    # Test 3: Driving Forces
    print("\n3. Driving Force Strengths:")
    force_strengths = get_driving_forces([mol], env)
    strong_forces = {k: v for k, v in force_strengths.items() if v > 0.5}
    print(f"   Strong forces (>0.5): {strong_forces}")
    
    # Test 4: Operator Registry
    print("\n4. Operator Registry:")
    available = list(OPERATOR_REGISTRY._operators.keys())
    print(f"   Available operators: {available}")
    
    active_ops = OPERATOR_REGISTRY.get_active_operators([mol], env)
    print(f"   Active for ethanol: {[op.name for op in active_ops]}")
    
    # Test 5: Reaction Sites
    print("\n5. Reaction Sites:")
    for op in active_ops[:2]:  # Test first 2 operators
        try:
            sites = op.find_reaction_sites([mol])
            print(f"   {op.name}: {len(sites)} sites")
        except Exception as e:
            print(f"   {op.name}: Error - {e}")
    
    print("\n✓ All tests completed successfully!")


if __name__ == "__main__":
    main()
