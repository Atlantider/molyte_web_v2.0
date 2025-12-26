"""
SEI Composition Analysis
========================
Classify species into inorganic/organic SEI components
Calculate reactivity metrics
"""

import json
from collections import Counter, defaultdict
from rdkit import Chem

def classify_sei_component(smiles):
    """
    Classify species as:
    - inorganic_sei: LiF, Li2O, Li2CO3
    - organic_sei: Polymers, alkyl carbonates
    - gas: CO2, CO, C2H4
    - intermediate: Everything else
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 'unknown', 0
        
        # Count atoms
        atoms = Counter([atom.GetSymbol() for atom in mol.GetAtoms()])
        num_heavy = mol.GetNumHeavyAtoms()
        
        # Inorganic SEI (no C-H bonds, contains Li)
        if atoms.get('Li', 0) > 0 and atoms.get('C', 0) == 0:
            # LiF, Li2O, LiOH, etc.
            return 'inorganic_sei', 10
        
        # Gases (volatiles)
        if smiles in ['O=C=O', '[C-]#[O+]', 'C=C', 'C']:
            return 'gas', 0
        
        # Li2CO3 (carbonate - inorganic)
        if 'Li' in smiles and 'C(=O)O' in smiles and 'C-C' not in smiles:
            # Simple carbonates
            if atoms.get('C', 0) == 1:
                return 'inorganic_sei', 9
        
        # Organic SEI
        if atoms.get('Li', 0) > 0 and atoms.get('C', 0) > 0:
            # Polymers/oligomers (>15 heavy atoms)
            if num_heavy > 15:
                return 'organic_sei_polymer', 8
            # Small Li-organic (alkyl carbonates, etc.)
            else:
                return 'organic_sei_small', 7
        
        # Pure organics (no Li - probably intermediates)
        if atoms.get('Li', 0) == 0 and atoms.get('C', 0) > 0:
            if '.' in smiles:
                return 'intermediate', 2
            return 'intermediate', 3
        
        return 'intermediate', 1
        
    except:
        return 'unknown', 0

def analyze_sei_composition(summary_file='simulation_summary.json'):
    """Analyze SEI composition and reactivity."""
    
    print("="*70)
    print("SEI COMPOSITION ANALYSIS")
    print("="*70)
    print()
    
    # Load data
    with open(summary_file, 'r') as f:
        data = json.load(f)
    
    # Build network graph to calculate degrees
    from collections import defaultdict
    out_degree = defaultdict(int)
    in_degree = defaultdict(int)
    
    # Parse reactions to get degrees
    for rxn in data['reaction_list']:
        if '->' in rxn:
            parts = rxn.split('->')
            if len(parts) == 2:
                reactants = parts[0].strip()
                products = parts[1].strip()
                
                # Simple parsing (assumes single reactant/product per reaction name)
                # This is approximate - would need better parsing
                for spec in data['species_list']:
                    if spec in reactants:
                        out_degree[spec] += 1
                    if spec in products:
                        in_degree[spec] += 1
    
    # Classify all species
    composition = {
        'inorganic_sei': [],
        'organic_sei_polymer': [],
        'organic_sei_small': [],
        'gas': [],
        'intermediate': []
    }
    
    species_info = {}
    
    for smiles in data['species_list']:
        category, priority = classify_sei_component(smiles)
        composition[category].append(smiles)
        
        # Calculate reactivity
        out_deg = out_degree.get(smiles, 0)
        in_deg = in_degree.get(smiles, 0)
        reactivity = out_deg / (in_deg + 1) if in_deg > 0 else out_deg
        
        species_info[smiles] = {
            'category': category,
            'priority': priority,
            'out_degree': out_deg,
            'in_degree': in_deg,
            'reactivity': reactivity,
            'is_stable': reactivity < 0.5,
            'is_active': reactivity > 2.0
        }
    
    # Print summary
    print("SEI COMPOSITION SUMMARY")
    print("-"*70)
    print()
    
    print("🟦 INORGANIC SEI (Bottom Layer):")
    print(f"   Total: {len(composition['inorganic_sei'])} species")
    if composition['inorganic_sei'][:10]:
        for smiles in composition['inorganic_sei'][:10]:
            info = species_info[smiles]
            stability = "稳定" if info['is_stable'] else "活性"
            print(f"   • {smiles[:40]} ({stability})")
    print()
    
    print("🟩 ORGANIC SEI (Outer Layer):")
    print(f"   Polymers: {len(composition['organic_sei_polymer'])}")
    print(f"   Small organics: {len(composition['organic_sei_small'])}")
    print(f"   Total: {len(composition['organic_sei_polymer']) + len(composition['organic_sei_small'])}")
    if composition['organic_sei_polymer'][:5]:
        print("   Top polymers:")
        for smiles in composition['organic_sei_polymer'][:5]:
            mol = Chem.MolFromSmiles(smiles)
            atoms = mol.GetNumHeavyAtoms() if mol else 0
            print(f"   • {smiles[:40]} ({atoms} atoms)")
    print()
    
    print("⚪ GAS PRODUCTS (Not in SEI):")
    print(f"   Total: {len(composition['gas'])} species")
    for smiles in composition['gas']:
        print(f"   • {smiles}")
    print()
    
    print("🟨 ACTIVE INTERMEDIATES (High Reactivity):")
    active = [(s, info) for s, info in species_info.items() if info['is_active']]
    print(f"   Total: {len(active)} species")
    for smiles, info in sorted(active, key=lambda x: -x[1]['reactivity'])[:10]:
        print(f"   • {smiles[:40]} (出度: {info['out_degree']}, 反应性: {info['reactivity']:.2f})")
    print()
    
    print("="*70)
    print("KEY FINDINGS")
    print("="*70)
    
    inorg_count = len(composition['inorganic_sei'])
    org_count = len(composition['organic_sei_polymer']) + len(composition['organic_sei_small'])
    total_sei = inorg_count + org_count
    
    print(f"SEI Total: {total_sei} species")
    print(f"  ├─ Inorganic: {inorg_count} ({inorg_count/total_sei*100:.1f}%)")
    print(f"  └─ Organic: {org_count} ({org_count/total_sei*100:.1f}%)")
    print()
    
    # Save results
    results = {
        'composition': {k: len(v) for k, v in composition.items()},
        'species_details': species_info,
        'inorganic_list': composition['inorganic_sei'],
        'organic_polymer_list': composition['organic_sei_polymer'],
        'organic_small_list': composition['organic_sei_small'],
        'gas_list': composition['gas']
    }
    
    output_file = 'sei_composition_analysis.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"Detailed results saved to: {output_file}")
    print()
    
    return results

if __name__ == "__main__":
    analyze_sei_composition()
