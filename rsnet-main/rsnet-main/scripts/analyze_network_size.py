"""
Network Analysis Script
========================
Analyzes the simulation output to check for:
1. Duplicate species (same SMILES)
2. Solvation complex explosion
3. Species distribution by generation
"""

import json
from collections import Counter, defaultdict
from rdkit import Chem

def analyze_simulation_output(json_path="simulation_summary.json"):
    """Analyze the simulation output for duplicates and patterns."""
    
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    species_list = data['species_list']
    total_species = len(species_list)
    
    print("="*60)
    print("NETWORK ANALYSIS REPORT")
    print("="*60)
    print(f"Total Species: {total_species}")
    print(f"Total Reactions: {data['reaction_count']}")
    print()
    
    # 1. Check for exact duplicates
    print("1. DUPLICATE CHECK")
    print("-"*60)
    species_counter = Counter(species_list)
    duplicates = {smiles: count for smiles, count in species_counter.items() if count > 1}
    
    if duplicates:
        print(f"⚠️  Found {len(duplicates)} duplicate SMILES:")
        for smiles, count in sorted(duplicates.items(), key=lambda x: -x[1])[:10]:
            print(f"   {smiles}: appears {count} times")
    else:
        print("✓ No exact duplicates found (SMILES are unique)")
    print()
    
    # 2. Analyze solvation complexes (dot-separated SMILES)
    print("2. SOLVATION COMPLEX ANALYSIS")
    print("-"*60)
    
    simple_species = []
    solvation_complexes = []
    
    for smiles in species_list:
        if '.' in smiles:
            solvation_complexes.append(smiles)
        else:
            simple_species.append(smiles)
    
    print(f"Simple molecules: {len(simple_species)} ({len(simple_species)/total_species*100:.1f}%)")
    print(f"Solvation complexes: {len(solvation_complexes)} ({len(solvation_complexes)/total_species*100:.1f}%)")
    print()
    
    # Analyze complexity of solvation complexes
    component_counts = Counter()
    for smiles in solvation_complexes:
        num_components = len(smiles.split('.'))
        component_counts[num_components] += 1
    
    print("Solvation complex distribution by number of components:")
    for num_comp in sorted(component_counts.keys()):
        count = component_counts[num_comp]
        print(f"   {num_comp} components: {count} species ({count/total_species*100:.1f}%)")
    print()
    
    # 3. Identify most common fragments in complexes
    print("3. MOST COMMON FRAGMENTS IN COMPLEXES")
    print("-"*60)
    fragment_counter = Counter()
    
    for smiles in solvation_complexes:
        fragments = smiles.split('.')
        for frag in fragments:
            fragment_counter[frag] += 1
    
    print("Top 10 most common fragments:")
    for frag, count in fragment_counter.most_common(10):
        print(f"   {frag}: appears in {count} complexes")
    print()
    
    # 4. Check for combinatorial explosion patterns
    print("4. COMBINATORIAL EXPLOSION CHECK")
    print("-"*60)
    
    # Group complexes by their sorted fragment composition
    composition_groups = defaultdict(list)
    for smiles in solvation_complexes:
        fragments = sorted(smiles.split('.'))
        composition_key = tuple(fragments)
        composition_groups[composition_key].append(smiles)
    
    # Find groups with multiple "duplicates" (same composition, different order)
    duplicate_compositions = {comp: smiles_list for comp, smiles_list in composition_groups.items() if len(smiles_list) > 1}
    
    if duplicate_compositions:
        print(f"⚠️  Found {len(duplicate_compositions)} compositions with multiple SMILES representations:")
        for comp, smiles_list in list(duplicate_compositions.items())[:5]:
            print(f"   Composition: {' + '.join(comp[:3])}{'...' if len(comp) > 3 else ''}")
            print(f"   Variants: {len(smiles_list)}")
            for s in smiles_list[:2]:
                print(f"      - {s[:80]}...")
    else:
        print("✓ No duplicate compositions found")
    print()
    
    # 5. Analyze species with Li+ and PF6-
    print("5. LITHIUM COMPLEX ANALYSIS")
    print("-"*60)
    
    li_complexes = [s for s in species_list if '[Li' in s or 'Li+' in s]
    pf6_complexes = [s for s in species_list if 'P' in s and 'F' in s]
    
    print(f"Species containing Li: {len(li_complexes)} ({len(li_complexes)/total_species*100:.1f}%)")
    print(f"Species containing PF6: {len(pf6_complexes)} ({len(pf6_complexes)/total_species*100:.1f}%)")
    print()
    
    # 6. Size distribution
    print("6. MOLECULAR SIZE DISTRIBUTION")
    print("-"*60)
    
    size_bins = defaultdict(int)
    for smiles in simple_species:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                num_heavy = mol.GetNumHeavyAtoms()
                bin_key = (num_heavy // 5) * 5  # Bin by 5s
                size_bins[bin_key] += 1
        except:
            pass
    
    print("Simple molecule size distribution (by heavy atoms):")
    for bin_start in sorted(size_bins.keys())[:10]:
        count = size_bins[bin_start]
        print(f"   {bin_start}-{bin_start+4} atoms: {count} species")
    print()
    
    # 7. Summary and recommendations
    print("="*60)
    print("SUMMARY AND RECOMMENDATIONS")
    print("="*60)
    
    if len(solvation_complexes) > total_species * 0.5:
        print("⚠️  ISSUE: >50% of species are solvation complexes")
        print("   This suggests the ClusteringOperator is generating too many combinations.")
        print()
        print("   RECOMMENDATIONS:")
        print("   1. Disable or limit the ClusteringOperator")
        print("   2. Set max_coordination to a lower value (e.g., 2)")
        print("   3. Filter complexes with >2 components")
    else:
        print("✓ Solvation complex ratio is reasonable")
    
    print()
    
    if duplicates:
        print("⚠️  ISSUE: Duplicate SMILES detected")
        print("   This indicates a bug in the deduplication logic.")
    else:
        print("✓ No duplicates - deduplication is working correctly")
    
    print()
    print(f"Network density: {data['reaction_count']/total_species:.1f} reactions per species")
    print()

if __name__ == "__main__":
    analyze_simulation_output()
