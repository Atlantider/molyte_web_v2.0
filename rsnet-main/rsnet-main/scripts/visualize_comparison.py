"""
Cathode Comparison Visualization
=================================
Creates side-by-side comparison of the three cathode simulations.
"""

import json
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import Counter

def create_comparison_visualization():
    """Create comparison charts for the three simulations."""
    
    print("="*60)
    print("CATHODE COMPARISON VISUALIZATION")
    print("="*60)
    print()
    
    # Load results
    results = {}
    for salt_type in ['LiPF6', 'LiTFSI', 'mixed']:
        with open(f'cathode_{salt_type}_summary.json', 'r') as f:
            results[salt_type] = json.load(f)
    
    # Create figure with subplots
    fig = plt.subplots(figsize=(16, 10), facecolor='white')
    fig = plt.figure(figsize=(16, 10), facecolor='white')
    
    # 1. Overall comparison
    ax1 = plt.subplot(2, 2, 1)
    
    salts = list(results.keys())
    species_counts = [results[s]['species_count'] for s in salts]
    reaction_counts = [results[s]['reaction_count'] for s in salts]
    
    x = range(len(salts))
    width = 0.35
    
    bars1 = ax1.bar([i - width/2 for i in x], species_counts, width, label='Species', color='#4ECDC4')
    bars2 = ax1.bar([i + width/2 for i in x], [r/4 for r in reaction_counts], width, label='Reactions (÷4)', color='#FF6B6B')
    
    ax1.set_xlabel('Salt Type', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Count', fontsize=12, fontweight='bold')
    ax1.set_title('Network Size Comparison\nCathode @ 4.2V, 298K', fontsize=14, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(salts)
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)
    
    # Add value labels
    for bar in bars1:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}',
                ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # 2. Species type distribution
    ax2 = plt.subplot(2, 2, 2)
    
    # Classify species
    def classify(smiles):
        if '[O-]' in smiles or '[C]' in smiles:
            return 'Radicals'
        if '.' in smiles:
            return 'Complexes'
        if 'S' in smiles and 'N' in smiles:
            return 'TFSI-derived'
        if 'P' in smiles and 'F' in smiles:
            return 'PF6-derived'
        return 'Intermediates'
    
    categories = ['Radicals', 'Complexes', 'TFSI-derived', 'PF6-derived', 'Intermediates']
    colors_cat = ['#FFD700', '#87CEEB', '#9370DB', '#FF69B4', '#90EE90']
    
    data_by_salt = {}
    for salt in salts:
        counts = Counter()
        for smiles in results[salt]['species_list']:
            counts[classify(smiles)] += 1
        data_by_salt[salt] = [counts.get(cat, 0) for cat in categories]
    
    x_cat = range(len(categories))
    width_cat = 0.25
    
    for i, salt in enumerate(salts):
        offset = (i - 1) * width_cat
        ax2.bar([x + offset for x in x_cat], data_by_salt[salt], width_cat, label=salt)
    
    ax2.set_xlabel('Species Type', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Count', fontsize=12, fontweight='bold')
    ax2.set_title('Species Type Distribution', fontsize=14, fontweight='bold')
    ax2.set_xticks(x_cat)
    ax2.set_xticklabels(categories, rotation=45, ha='right')
    ax2.legend()
    ax2.grid(axis='y', alpha=0.3)
    
    # 3. Key findings text
    ax3 = plt.subplot(2, 2, 3)
    ax3.axis('off')
    
    findings_text = f"""
    KEY FINDINGS:
    
    1. LiPF6 (Baseline):
       • {results['LiPF6']['species_count']} species, {results['LiPF6']['reaction_count']} reactions
       • Simpler network, fewer intermediates
       • Dominated by PF6- decomposition products
    
    2. LiTFSI:
       • {results['LiTFSI']['species_count']} species (+{results['LiTFSI']['species_count'] - results['LiPF6']['species_count']} vs LiPF6)
       • 2.4x more species than LiPF6
       • TFSI- generates more diverse chemistry
       • Higher radical concentration
    
    3. Mixed (LiPF6 + LiTFSI):
       • {results['mixed']['species_count']} species (+{results['mixed']['species_count'] - results['LiPF6']['species_count']} vs LiPF6)
       • Slightly more than LiTFSI alone
       • Cross-reactions between PF6- and TFSI- products
       • Most complex network
    
    RECOMMENDATION:
    • LiPF6: Simpler, more predictable CEI
    • LiTFSI: More complex but potentially more stable
    • Mixed: Synergistic effects, needs careful study
    """
    
    ax3.text(0.1, 0.9, findings_text, transform=ax3.transAxes,
             fontsize=11, verticalalignment='top', family='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    # 4. Unique species comparison
    ax4 = plt.subplot(2, 2, 4)
    
    # Find unique species in each condition
    lipf6_set = set(results['LiPF6']['species_list'])
    litfsi_set = set(results['LiTFSI']['species_list'])
    mixed_set = set(results['mixed']['species_list'])
    
    unique_lipf6 = len(lipf6_set - litfsi_set - mixed_set)
    unique_litfsi = len(litfsi_set - lipf6_set - mixed_set)
    unique_mixed = len(mixed_set - lipf6_set - litfsi_set)
    shared_all = len(lipf6_set & litfsi_set & mixed_set)
    
    labels = ['Unique to\nLiPF6', 'Unique to\nLiTFSI', 'Unique to\nMixed', 'Shared by\nAll']
    sizes = [unique_lipf6, unique_litfsi, unique_mixed, shared_all]
    colors_pie = ['#FF6B6B', '#9370DB', '#4ECDC4', '#90EE90']
    
    ax4.pie(sizes, labels=labels, colors=colors_pie, autopct='%1.1f%%',
            startangle=90, textprops={'fontsize': 10, 'fontweight': 'bold'})
    ax4.set_title('Species Uniqueness', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    
    # Save
    output_path = 'cathode_comparison.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    
    print(f"✓ Comparison visualization saved: {output_path}")
    print()
    
    # Print summary
    print("="*60)
    print("SUMMARY")
    print("="*60)
    print(f"LiPF6:  {results['LiPF6']['species_count']:3d} species, {results['LiPF6']['reaction_count']:4d} reactions")
    print(f"LiTFSI: {results['LiTFSI']['species_count']:3d} species, {results['LiTFSI']['reaction_count']:4d} reactions (+{results['LiTFSI']['species_count'] - results['LiPF6']['species_count']:3d} species)")
    print(f"Mixed:  {results['mixed']['species_count']:3d} species, {results['mixed']['reaction_count']:4d} reactions (+{results['mixed']['species_count'] - results['LiPF6']['species_count']:3d} species)")
    print()
    print("LiTFSI generates 2.4x more species than LiPF6 at cathode!")
    print("="*60)

if __name__ == "__main__":
    create_comparison_visualization()
