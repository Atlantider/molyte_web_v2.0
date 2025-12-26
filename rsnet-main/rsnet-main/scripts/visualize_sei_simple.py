"""
SEI Network Visualization (Simplified)
=======================================
Creates a clear, publication-quality visualization using only matplotlib and networkx.
"""

import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict

def classify_species(smiles):
    """Classify species by chemical type."""
    # Initial reactants
    if smiles in ['O=C1OCCO1', 'COC(=O)OC', '[Li+]', 'F[P-](F)(F)(F)(F)F']:
        return 'reactant', '#FF4444', 'Initial\nReactants'
    
    # Li metal
    if smiles in ['[Li]', '[Li][Li]']:
        return 'li_metal', '#00CED1', 'Li Metal'
    
    # Radicals
    if '[O-]' in smiles or ('[C]' in smiles and '(' not in smiles.split('[C]')[1][:3]):
        return 'radical', '#FFD700', 'Radicals'
    
    # SEI products
    if 'Li' in smiles and 'F' in smiles and '.' not in smiles:
        return 'lif', '#FF69B4', 'LiF Products'
    
    if ('C(=O)O' in smiles or 'CO3' in smiles) and '.' not in smiles:
        return 'carbonate', '#9370DB', 'Carbonates'
    
    # Solvation complexes
    if '.' in smiles:
        return 'complex', '#87CEEB', 'Solvation\nComplexes'
    
    # Intermediates
    return 'intermediate', '#90EE90', 'Intermediates'

def create_simplified_visualization(
    graphml_path='sei_full_network.graphml',
    output_path='sei_network_simple.png',
    max_nodes=200
):
    """Create simplified but informative visualization."""
    
    print("="*60)
    print("SEI NETWORK VISUALIZATION")
    print("="*60)
    print()
    print("Loading network...")
    
    G = nx.read_graphml(graphml_path)
    
    # Separate nodes
    molecule_nodes = [n for n, d in G.nodes(data=True) if d.get('type') != 'reaction']
    reaction_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'reaction']
    
    print(f"  Total molecules: {len(molecule_nodes)}")
    print(f"  Total reactions: {len(reaction_nodes)}")
    print()
    
    # Create molecule-only network (remove reaction nodes for clarity)
    print("Creating molecule-only network...")
    mol_graph = nx.DiGraph()
    
    # Add molecules
    for mol in molecule_nodes:
        mol_graph.add_node(mol, **G.nodes[mol])
    
    # Add direct molecule-to-molecule edges (through reactions)
    for rxn in reaction_nodes:
        predecessors = [n for n in G.predecessors(rxn) if n in molecule_nodes]
        successors = [n for n in G.successors(rxn) if n in molecule_nodes]
        
        for pred in predecessors:
            for succ in successors:
                if pred != succ:
                    mol_graph.add_edge(pred, succ)
    
    # Filter to most important nodes if too large
    if len(mol_graph.nodes()) > max_nodes:
        print(f"Network too large. Filtering to top {max_nodes} molecules...")
        
        # Calculate importance
        pagerank = nx.pagerank(mol_graph)
        important_nodes = sorted(
            mol_graph.nodes(),
            key=lambda n: pagerank.get(n, 0),
            reverse=True
        )[:max_nodes]
        
        mol_graph = mol_graph.subgraph(important_nodes).copy()
        print(f"  Filtered to {len(mol_graph.nodes())} molecules")
    
    print()
    print("Classifying species...")
    
    # Classify all nodes
    node_data = {}
    species_counts = defaultdict(int)
    
    for node in mol_graph.nodes():
        species_type, color, label = classify_species(node)
        node_data[node] = {
            'type': species_type,
            'color': color,
            'label': label
        }
        species_counts[species_type] += 1
    
    # Print distribution
    print()
    print("Species distribution:")
    for stype, count in sorted(species_counts.items(), key=lambda x: -x[1]):
        print(f"  {stype:15s}: {count:4d}")
    
    # Create visualization
    print()
    print("Creating visualization...")
    
    fig, ax = plt.subplots(figsize=(16, 12), facecolor='white')
    
    # Compute layout
    print("  Computing layout...")
    pos = nx.spring_layout(mol_graph, k=1.5, iterations=50, seed=42)
    
    # Draw edges
    print("  Drawing edges...")
    nx.draw_networkx_edges(
        mol_graph, pos,
        edge_color='#DDDDDD',
        alpha=0.4,
        width=0.8,
        arrows=True,
        arrowsize=8,
        arrowstyle='->',
        ax=ax,
        connectionstyle='arc3,rad=0.1'
    )
    
    # Draw nodes by type
    print("  Drawing nodes...")
    for species_type in set(d['type'] for d in node_data.values()):
        nodes_of_type = [n for n, d in node_data.items() if d['type'] == species_type]
        if nodes_of_type:
            color = node_data[nodes_of_type[0]]['color']
            
            # Size by generation
            sizes = []
            for n in nodes_of_type:
                gen = mol_graph.nodes[n].get('generation', 0)
                sizes.append(200 + int(gen) * 40)
            
            nx.draw_networkx_nodes(
                mol_graph, pos,
                nodelist=nodes_of_type,
                node_color=color,
                node_size=sizes,
                alpha=0.85,
                edgecolors='black',
                linewidths=1.5,
                ax=ax
            )
    
    # Add labels for important nodes
    print("  Adding labels...")
    important_nodes = sorted(
        mol_graph.nodes(),
        key=lambda n: mol_graph.degree(n),
        reverse=True
    )[:15]
    
    labels = {}
    for node in important_nodes:
        label = mol_graph.nodes[node].get('label', node)
        if len(label) > 15:
            label = label[:12] + '...'
        labels[node] = label
    
    nx.draw_networkx_labels(
        mol_graph, pos,
        labels=labels,
        font_size=7,
        font_weight='bold',
        font_color='black',
        ax=ax
    )
    
    # Create legend
    print("  Creating legend...")
    legend_elements = []
    seen_types = set()
    
    for node, data in node_data.items():
        if data['type'] not in seen_types:
            legend_elements.append(
                mpatches.Patch(
                    color=data['color'],
                    label=data['label']
                )
            )
            seen_types.add(data['type'])
    
    ax.legend(
        handles=legend_elements,
        loc='upper right',
        fontsize=11,
        framealpha=0.95,
        edgecolor='black',
        title='Species Types',
        title_fontsize=12
    )
    
    # Add title
    ax.set_title(
        'SEI Formation Network\n'
        f'Anode: 0.1V vs Li/Li+, 298K | {len(mol_graph.nodes())} Species',
        fontsize=16,
        fontweight='bold',
        pad=20
    )
    
    # Add info box
    info_text = (
        f"Network Statistics:\n"
        f"• Molecules: {len(mol_graph.nodes())}\n"
        f"• Connections: {len(mol_graph.edges())}\n"
        f"• Avg. Degree: {sum(dict(mol_graph.degree()).values())/len(mol_graph.nodes()):.1f}"
    )
    
    ax.text(
        0.02, 0.02, info_text,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment='bottom',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='black')
    )
    
    ax.axis('off')
    plt.tight_layout()
    
    # Save
    print()
    print(f"Saving to {output_path}...")
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print()
    print("="*60)
    print(f"✓ Visualization saved to: {output_path}")
    print("="*60)
    
    return mol_graph

if __name__ == "__main__":
    G = create_simplified_visualization()
