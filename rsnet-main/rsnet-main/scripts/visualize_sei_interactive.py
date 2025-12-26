"""
Interactive SEI Network Visualization
======================================
Creates an interactive HTML visualization using Pyvis for better exploration.
"""

from pyvis.network import Network
import networkx as nx
from collections import defaultdict
from rdkit import Chem

def classify_species(smiles):
    """Classify species by chemical type."""
    if smiles in ['O=C1OCCO1', 'COC(=O)OC', '[Li+]', 'F[P-](F)(F)(F)(F)F']:
        return 'reactant', '#FF6B6B'
    if smiles in ['[Li]', '[Li][Li]']:
        return 'li_metal', '#4ECDC4'
    if '[O-]' in smiles or '[C]' in smiles or '[CH]' in smiles:
        return 'radical', '#FFE66D'
    if '.' in smiles:
        return 'complex', '#95E1D3'
    if 'Li' in smiles and 'F' in smiles:
        return 'lif_product', '#F38181'
    if 'C(=O)O' in smiles:
        return 'carbonate_product', '#AA96DA'
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol and mol.GetNumHeavyAtoms() > 15:
            return 'oligomer', '#FCBAD3'
    except:
        pass
    
    return 'intermediate', '#A8E6CF'

def create_interactive_visualization(
    graphml_path='sei_full_network.graphml',
    output_path='sei_network_interactive.html',
    max_nodes=300
):
    """Create interactive HTML visualization."""
    
    print("Loading network...")
    G = nx.read_graphml(graphml_path)
    
    # Separate nodes
    molecule_nodes = [n for n, d in G.nodes(data=True) if d.get('type') != 'reaction']
    reaction_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'reaction']
    
    print(f"Total: {len(molecule_nodes)} molecules, {len(reaction_nodes)} reactions")
    
    # Filter to most important nodes
    if len(molecule_nodes) > max_nodes:
        print(f"Filtering to top {max_nodes} molecules by connectivity...")
        importance = nx.degree_centrality(G)
        important_molecules = sorted(
            molecule_nodes,
            key=lambda n: importance.get(n, 0),
            reverse=True
        )[:max_nodes]
        
        important_reactions = set()
        for rxn in reaction_nodes:
            neighbors = list(G.neighbors(rxn))
            if any(n in important_molecules for n in neighbors):
                important_reactions.add(rxn)
        
        important_nodes = set(important_molecules) | important_reactions
        G = G.subgraph(important_nodes).copy()
        
        molecule_nodes = [n for n in G.nodes() if n in important_molecules]
        reaction_nodes = [n for n in G.nodes() if n in important_reactions]
    
    # Create Pyvis network
    print("Creating interactive visualization...")
    net = Network(
        height='900px',
        width='100%',
        bgcolor='#ffffff',
        font_color='black',
        directed=True
    )
    
    # Configure physics
    net.set_options("""
    {
      "physics": {
        "forceAtlas2Based": {
          "gravitationalConstant": -50,
          "centralGravity": 0.01,
          "springLength": 100,
          "springConstant": 0.08
        },
        "maxVelocity": 50,
        "solver": "forceAtlas2Based",
        "timestep": 0.35,
        "stabilization": {"iterations": 150}
      },
      "interaction": {
        "hover": true,
        "tooltipDelay": 100
      }
    }
    """)
    
    # Add nodes
    print("Adding nodes...")
    for node in G.nodes():
        if node in reaction_nodes:
            # Reaction node
            operator = G.nodes[node].get('operator', 'reaction')
            net.add_node(
                node,
                label='',
                title=f"Reaction: {operator}",
                color='#D3D3D3',
                size=5,
                shape='diamond'
            )
        else:
            # Molecule node
            species_type, color = classify_species(node)
            label = G.nodes[node].get('label', node)
            
            # Truncate long labels
            if len(label) > 30:
                display_label = label[:27] + '...'
            else:
                display_label = label
            
            # Size by generation
            gen = int(G.nodes[node].get('generation', 0))
            size = 15 + gen * 3
            
            # Create tooltip
            tooltip = f"<b>{label}</b><br>Type: {species_type}<br>Generation: {gen}<br>SMILES: {node}"
            
            net.add_node(
                node,
                label=display_label if gen == 0 or species_type in ['reactant', 'li_metal'] else '',
                title=tooltip,
                color=color,
                size=size,
                shape='dot'
            )
    
    # Add edges
    print("Adding edges...")
    for edge in G.edges():
        net.add_edge(edge[0], edge[1], color='#CCCCCC', width=0.5)
    
    # Add legend as HTML
    legend_html = """
    <div style="position: fixed; top: 10px; right: 10px; background: white; 
                padding: 15px; border: 2px solid #333; border-radius: 5px; 
                font-family: Arial; font-size: 12px; z-index: 1000;">
        <h3 style="margin-top: 0;">SEI Network Legend</h3>
        <div><span style="color: #FF6B6B;">●</span> Initial Reactants</div>
        <div><span style="color: #4ECDC4;">●</span> Li Metal</div>
        <div><span style="color: #FFE66D;">●</span> Radicals</div>
        <div><span style="color: #95E1D3;">●</span> Solvation Complexes</div>
        <div><span style="color: #F38181;">●</span> LiF Products</div>
        <div><span style="color: #AA96DA;">●</span> Carbonate Products</div>
        <div><span style="color: #FCBAD3;">●</span> Oligomers</div>
        <div><span style="color: #A8E6CF;">●</span> Intermediates</div>
        <div><span style="color: #D3D3D3;">◆</span> Reactions</div>
        <hr>
        <div><b>Anode: 0.1V vs Li/Li+, 298K</b></div>
        <div>{} molecules, {} reactions</div>
    </div>
    """.format(len(molecule_nodes), len(reaction_nodes))
    
    # Save
    print(f"Saving to {output_path}...")
    net.save_graph(output_path)
    
    # Add legend to HTML
    with open(output_path, 'r', encoding='utf-8') as f:
        html = f.read()
    
    html = html.replace('</body>', legend_html + '</body>')
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)
    
    print(f"✓ Interactive visualization saved!")
    print(f"  Open {output_path} in a web browser to explore")
    
    # Statistics
    print("\n" + "="*60)
    print("SPECIES DISTRIBUTION")
    print("="*60)
    
    species_counts = defaultdict(int)
    for node in molecule_nodes:
        species_type, _ = classify_species(node)
        species_counts[species_type] += 1
    
    for species_type, count in sorted(species_counts.items(), key=lambda x: -x[1]):
        print(f"{species_type:20s}: {count:4d}")

if __name__ == "__main__":
    create_interactive_visualization()
