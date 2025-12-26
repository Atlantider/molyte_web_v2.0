"""
Generation-Based Network Views
================================
Creates separate visualizations for each generation to show network evolution.
"""

import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict

def classify_species(smiles):
    """Classify species by chemical type."""
    if smiles in ['O=C1OCCO1', 'COC(=O)OC', '[Li+]', 'F[P-](F)(F)(F)(F)F']:
        return 'reactant', '#FF4444'
    if smiles in ['[Li]', '[Li][Li]']:
        return 'li_metal', '#00CED1'
    if '[O-]' in smiles or ('[C]' in smiles and '(' not in smiles.split('[C]')[1][:3] if '[C]' in smiles else False):
        return 'radical', '#FFD700'
    if 'Li' in smiles and 'F' in smiles and '.' not in smiles:
        return 'lif', '#FF69B4'
    if ('C(=O)O' in smiles or 'CO3' in smiles) and '.' not in smiles:
        return 'carbonate', '#9370DB'
    if '.' in smiles:
        return 'complex', '#87CEEB'
    return 'intermediate', '#90EE90'

def create_generation_views(graphml_path='sei_full_network.graphml'):
    """Create separate views for each generation."""
    
    print("="*60)
    print("GENERATION-BASED NETWORK VIEWS")
    print("="*60)
    print()
    
    # Load network
    print("Loading network...")
    G = nx.read_graphml(graphml_path)
    
    molecule_nodes = [n for n, d in G.nodes(data=True) if d.get('type') != 'reaction']
    
    # Group by generation
    by_generation = defaultdict(list)
    for node in molecule_nodes:
        gen = int(G.nodes[node].get('generation', 0))
        by_generation[gen].append(node)
    
    print(f"Generations found: {sorted(by_generation.keys())}")
    print()
    
    # Create subgraph for each generation
    for gen in sorted(by_generation.keys())[:5]:  # First 5 generations
        print(f"Creating view for Generation {gen}...")
        
        # Get nodes up to this generation
        nodes_up_to_gen = []
        for g in range(gen + 1):
            nodes_up_to_gen.extend(by_generation[g])
        
        # Create subgraph
        subG = G.subgraph(nodes_up_to_gen).copy()
        
        # Remove reaction nodes
        mol_subG = nx.DiGraph()
        for node in subG.nodes():
            if node in molecule_nodes:
                mol_subG.add_node(node, **subG.nodes[node])
        
        # Add edges
        reaction_nodes = [n for n in subG.nodes() if subG.nodes[n].get('type') == 'reaction']
        for rxn in reaction_nodes:
            preds = [n for n in subG.predecessors(rxn) if n in molecule_nodes]
            succs = [n for n in subG.successors(rxn) if n in molecule_nodes]
            for p in preds:
                for s in succs:
                    if p != s:
                        mol_subG.add_edge(p, s)
        
        # Limit to most important nodes if too large
        if len(mol_subG.nodes()) > 100:
            pagerank = nx.pagerank(mol_subG)
            important = sorted(mol_subG.nodes(), key=lambda n: pagerank.get(n, 0), reverse=True)[:100]
            mol_subG = mol_subG.subgraph(important).copy()
        
        # Create visualization
        fig, ax = plt.subplots(figsize=(14, 10), facecolor='white')
        
        # Layout
        pos = nx.spring_layout(mol_subG, k=1.2, iterations=50, seed=42)
        
        # Draw edges
        nx.draw_networkx_edges(
            mol_subG, pos,
            edge_color='#DDDDDD',
            alpha=0.5,
            width=1,
            arrows=True,
            arrowsize=10,
            ax=ax
        )
        
        # Draw nodes by type
        node_data = {}
        for node in mol_subG.nodes():
            species_type, color = classify_species(node)
            node_data[node] = {'type': species_type, 'color': color}
        
        for species_type in set(d['type'] for d in node_data.values()):
            nodes_of_type = [n for n, d in node_data.items() if d['type'] == species_type]
            if nodes_of_type:
                color = node_data[nodes_of_type[0]]['color']
                
                # Highlight new species in this generation
                sizes = []
                alphas = []
                for n in nodes_of_type:
                    node_gen = int(mol_subG.nodes[n].get('generation', 0))
                    if node_gen == gen:
                        sizes.append(400)  # Larger for new species
                        alphas.append(1.0)
                    else:
                        sizes.append(200)
                        alphas.append(0.6)
                
                nx.draw_networkx_nodes(
                    mol_subG, pos,
                    nodelist=nodes_of_type,
                    node_color=color,
                    node_size=sizes,
                    alpha=0.85,
                    edgecolors='black',
                    linewidths=2,
                    ax=ax
                )
        
        # Add labels for new species
        new_species = [n for n in mol_subG.nodes() if int(mol_subG.nodes[n].get('generation', 0)) == gen]
        labels = {}
        for node in new_species[:10]:  # Top 10 new species
            label = mol_subG.nodes[node].get('label', node)
            if len(label) > 12:
                label = label[:9] + '...'
            labels[node] = label
        
        nx.draw_networkx_labels(
            mol_subG, pos,
            labels=labels,
            font_size=7,
            font_weight='bold',
            ax=ax
        )
        
        # Title
        ax.set_title(
            f'Generation {gen} Network Evolution\n'
            f'{len(new_species)} New Species (highlighted), {len(mol_subG.nodes())} Total',
            fontsize=14,
            fontweight='bold',
            pad=15
        )
        
        # Stats box
        stats_text = (
            f"Generation {gen} Statistics:\n"
            f"• New species: {len(new_species)}\n"
            f"• Total species: {len(mol_subG.nodes())}\n"
            f"• Connections: {len(mol_subG.edges())}"
        )
        
        ax.text(
            0.02, 0.98, stats_text,
            transform=ax.transAxes,
            fontsize=9,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='black')
        )
        
        ax.axis('off')
        plt.tight_layout()
        
        # Save
        output_path = f'sei_network_gen{gen}.png'
        plt.savefig(output_path, dpi=200, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print(f"  ✓ Saved: {output_path}")
    
    print()
    print("="*60)
    print("✓ Generation views created!")
    print("  Files: sei_network_gen0.png to sei_network_gen4.png")
    print("="*60)

if __name__ == "__main__":
    create_generation_views()
