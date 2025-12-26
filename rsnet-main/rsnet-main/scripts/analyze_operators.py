"""
Analyze Which Operators Were Active in Simulation
===================================================
Check which operators were used and look for ring opening/polymerization products
"""

import json
from collections import Counter
from rdkit import Chem

def analyze_operators_and_reactions():
    """Analyze simulation to find what operators were active."""
    
    print("="*60)
    print("OPERATOR AND REACTION ANALYSIS")
    print("="*60)
    print()
    
    # Load simulation summary
    with open('simulation_summary.json', 'r') as f:
        data = json.load(f)
    
    print(f"Total species: {data['species_count']}")
    print(f"Total reactions: {data['reaction_count']}")
    print()
    
    # Analyze reaction types from names
    rxn_types = Counter()
    for rxn_name in data['reaction_list']:
        # Extract operator type from reaction name
        if ':' in rxn_name:
            op_type = rxn_name.split(':')[0]
            rxn_types[op_type] += 1
    
    print("Reaction Types Found:")
    print("-" * 60)
    for rxn_type, count in sorted(rxn_types.items(), key=lambda x: -x[1]):
        print(f"  {rxn_type:<30} {count:>5} reactions")
    print()
    
    # Check for ring-opening products
    print("Checking for ring-opening products...")
    print("-" * 60)
    
    ring_opened = []
    linear_products = []
    
    for smiles in data['species_list']:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            
            # Check if it's a linear chain (no rings)
            ring_info = mol.GetRingInfo()
            num_rings = ring_info.NumRings()
            
            # Check for specific ring-opening signatures
            # EC ring opening: O=C1OCCO1 -> OCCO (ethylene glycol derivatives)
            if 'OCCO' in smiles and 'O=C1OCCO1' not in smiles and num_rings == 0:
                ring_opened.append(smiles)
            # DMC doesn't have a ring, but can form linear chains
            elif num_rings == 0 and len(smiles) > 10 and '.' not in smiles:
                linear_products.append(smiles)
        except:
            continue
    
    print(f"Found {len(ring_opened)} potential ring-opened products")
    if ring_opened[:5]:
        print("Examples:")
        for ex in ring_opened[:5]:
            print(f"  {ex}")
    print()
    
    # Check for polymerization products
    print("Checking for polymerization products...")
    print("-" * 60)
    
    polymers = []
    oligomers = []
    
    for smiles in data['species_list']:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            
            num_heavy_atoms = mol.GetNumHeavyAtoms()
            
            # Oligomers: repeated units, >20 heavy atoms
            if num_heavy_atoms > 20 and '.' not in smiles:
                oligomers.append((smiles, num_heavy_atoms))
            
            # Check for polymer-like structures
            if num_heavy_atoms > 30:
                polymers.append((smiles, num_heavy_atoms))
        except:
            continue
    
    print(f"Found {len(oligomers)} oligomers (>20 atoms)")
    print(f"Found {len(polymers)} large polymers (>30 atoms)")
    
    if oligomers:
        print("\nTop 5 largest oligomers:")
        for smiles, atoms in sorted(oligomers, key=lambda x: -x[1])[:5]:
            print(f"  {atoms} atoms: {smiles[:80]}...")
    print()
    
    # Summary
    print("="*60)
    print("SUMMARY")
    print("="*60)
    
    if 'RingOpening' in rxn_types:
        print(f"✓ Ring opening IS active ({rxn_types['RingOpening']} reactions)")
    else:
        print("✗ Ring opening NOT found in reaction names")
    
    if 'Polymerization' in rxn_types:
        print(f"✓ Polymerization IS active ({rxn_types['Polymerization']} reactions)")
    else:
        print("✗ Polymerization NOT found in reaction names")
    
    print()
    print("Possible reasons if missing:")
    print("  1. Operators not enabled in config")
    print("  2. Conditions don't favor these reactions")
    print("  3. Energy cutoff filtered them out")
    print("  4. Products were generated but filtered from visualization (top 80)")
    print()
    
    # Check if we have the products but they're filtered
    if ring_opened or oligomers:
        print("⚠️  Ring-opened/oligomer products EXIST but may be filtered out!")
        print(f"   Ring-opened: {len(ring_opened)} species")
        print(f"   Oligomers: {len(oligomers)} species")
        print()
        print("   Recommendation: Create specialized view showing these products")
    
    return {
        'reaction_types': rxn_types,
        'ring_opened': ring_opened,
        'oligomers': oligomers,
        'polymers': polymers
    }

if __name__ == "__main__":
    results = analyze_operators_and_reactions()
