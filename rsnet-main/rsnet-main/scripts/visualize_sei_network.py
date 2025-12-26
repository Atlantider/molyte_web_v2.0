"""
SEI Network Visualization for Electrochemists
==============================================
Creates an interactive visualization of the reaction network with:
- Color coding by species type (reactants, radicals, SEI products)
- Node size by generation/importance
- Edge thickness by reaction frequency
- Electrochemical context (voltage, reduction/oxidation)
"""

import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from rdkit import Chem
from collections import defaultdict

def classify_species(smiles):
    """Classify species by chemical type for color coding."""
    
    # Initial reactants
    if smiles in ['O=C1OCCO1', 'COC(=O)OC', '[Li+]', 'F[P-](F)(F)(F)(F)F']:
        return 'reactant'
    
    # Li metal and dimers
    if smiles in ['[Li]', '[Li][Li]']:
        return 'li_metal'
    
    # Radicals (contains unpaired electrons)
    if '[O-]' in smiles or '[C]' in smiles or '[CH]' in smiles:
        return 'radical'
    
    # Solvation complexes (dot-separated)
    if '.' in smiles:
        return 'complex'
    
    # SEI products (carbonates, fluorides)
    if 'Li' in smiles and 'F' in smiles:
        return 'lif_product'
    if 'C(=O)O' in smiles or 'CO3' in smiles:
        return 'carbonate_product'
    
    # Oligomers (large molecules)
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol and mol.GetNumHeavyAtoms() > 15:
            return 'oligomer'
    except:
        pass
    
    return 'intermediate'

def create_sei_visualization(graphml_path='sei_full_network.graphml', 
                             output_path='sei_network_viz.png',
                             max_nodes=500):
    """Create chemically-informed visualization."""
    
    print("Loading network...")
    G = nx.read_graphml(graphml_path)
    
    # Separate molecule and reaction nodes
    molecule_nodes = [n for n, d in G.nodes(data=True) if d.get('type') != 'reaction']
    reaction_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'reaction']
    
    print(f"Total nodes: {len(G.nodes())}")
    print(f"  Molecules: {len(molecule_nodes)}")
    print(f"  Reactions: {len(reaction_nodes)}")
    
    # Filter to most important nodes if network is too large
    if len(molecule_nodes) > max_nodes:
        print(f"\nNetwork too large. Filtering to {max_nodes} most important molecules...")
        
        # Calculate importance (degree centrality)
        importance = nx.degree_centrality(G)
        important_molecules = sorted(
            molecule_nodes, 
            key=lambda n: importance.get(n, 0), 
            reverse=True
        )[:max_nodes]
        
        # Get reactions involving these molecules
        important_reactions = set()
        for rxn in reaction_nodes:
            neighbors = list(G.neighbors(rxn))
            if any(n in important_molecules for n in neighbors):
                important_reactions.add(rxn)
        
        # Create subgraph
        important_nodes = set(important_molecules) | important_reactions
        G = G.subgraph(important_nodes).copy()
        
        molecule_nodes = [n for n in G.nodes() if n in important_molecules]
        reaction_nodes = [n for n in G.nodes() if n in important_reactions]
        
        print(f"Filtered network: {len(molecule_nodes)} molecules, {len(reaction_nodes)} reactions")
    
    # Classify molecules
    print("\nClassifying species...")
    node_colors = {}
    node_sizes = {}
    
    color_map = {
        'reactant': '#FF6B6B',      # Red - initial reactants
        'li_metal': '#4ECDC4',      # Cyan - Li metal
        'radical': '#FFE66D',       # Yellow - radicals
        'complex': '#95E1D3',       # Light teal - solvation complexes
        'lif_product': '#F38181',   # Pink - LiF products
        'carbonate_product': '#AA96DA', # Purple - carbonate products
        'oligomer': '#FCBAD3',      # Light pink - oligomers
        'intermediate': '#A8E6CF',  # Light green - intermediates
        'reaction': '#D3D3D3'       # Gray - reaction nodes
    }
    
    for node in G.nodes():
        if node in reaction_nodes:
            node_colors[node] = color_map['reaction']
            node_sizes[node] = 50
        else:
            species_type = classify_species(node)
            node_colors[node] = color_map[species_type]
            
            # Size by generation (if available)
            gen = G.nodes[node].get('generation', 0)
            node_sizes[node] = 100 + gen * 20
    
    # Create layout
    print("Computing layout...")
    # Use spring layout for better clustering
    pos = nx.spring_layout(G, k=0.5, iterations=50, seed=42)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(20, 16), facecolor='white')
    
    # Draw edges (reactions)
    print("Drawing edges...")
    nx.draw_networkx_edges(
        G, pos, 
        edge_color='#CCCCCC',
        alpha=0.3,
        width=0.5,
        arrows=True,
        arrowsize=5,
        ax=ax
    )
    
    # Draw nodes
    print("Drawing nodes...")
    for node_type, color in color_map.items():
        nodes_of_type = [n for n in G.nodes() if node_colors.get(n) == color]
        if nodes_of_type:
            nx.draw_networkx_nodes(
                G, pos,
                nodelist=nodes_of_type,
                node_color=color,
                node_size=[node_sizes.get(n, 100) for n in nodes_of_type],
                alpha=0.8,
                edgecolors='black',
                linewidths=0.5,
                ax=ax
            )
    
    # Add labels for important nodes only
    print("Adding labels...")
    important_labels = {}
    for node in molecule_nodes[:20]:  # Top 20 most connected
        label = G.nodes[node].get('label', node)
        if len(label) > 20:
            label = label[:17] + '...'
        important_labels[node] = label
    
    nx.draw_networkx_labels(
        G, pos,
        labels=important_labels,
        font_size=6,
        font_weight='bold',
        ax=ax
    )
    
    # Create legend
    legend_elements = [
        mpatches.Patch(color=color_map['reactant'], label='Initial Reactants (EC, DMC, Li+, PF6-)'),
        mpatches.Patch(color=color_map['li_metal'], label='Li Metal (Li, Li2)'),
        mpatches.Patch(color=color_map['radical'], label='Radicals'),
        mpatches.Patch(color=color_map['complex'], label='Solvation Complexes'),
        mpatches.Patch(color=color_map['lif_product'], label='LiF Products'),
        mpatches.Patch(color=color_map['carbonate_product'], label='Carbonate Products'),
        mpatches.Patch(color=color_map['oligomer'], label='Oligomers (>15 atoms)'),
        mpatches.Patch(color=color_map['intermediate'], label='Intermediates'),
        mpatches.Patch(color=color_map['reaction'], label='Reactions', alpha=0.5),
    ]
    
    ax.legend(
        handles=legend_elements,
        loc='upper left',
        fontsize=10,
        framealpha=0.9
    )
    
    # Add title and context
    ax.set_title(
        'SEI Formation Network at Anode (0.1V vs Li/Li+, 298K)\n'
        f'{len(molecule_nodes)} Species, {len(reaction_nodes)} Reactions',
        fontsize=16,
        fontweight='bold',
        pad=20
    )
    
    ax.axis('off')
    plt.tight_layout()
    
    # Save
    print(f"\nSaving visualization to {output_path}...")
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print("✓ Visualization saved!")
    
    # Print statistics
    print("\n" + "="*60)
    print("NETWORK STATISTICS")
    print("="*60)
    
    species_counts = defaultdict(int)
    for node in molecule_nodes:
        species_type = classify_species(node)
        species_counts[species_type] += 1
    
    for species_type, count in sorted(species_counts.items(), key=lambda x: -x[1]):
        print(f"{species_type:20s}: {count:4d} ({count/len(molecule_nodes)*100:5.1f}%)")
    
    print()
    print(f"Average degree: {sum(dict(G.degree()).values()) / len(G.nodes()):.2f}")
    print(f"Network density: {nx.density(G):.4f}")
    
    return G

if __name__ == "__main__":
    G = create_sei_visualization()
    plt.show()
