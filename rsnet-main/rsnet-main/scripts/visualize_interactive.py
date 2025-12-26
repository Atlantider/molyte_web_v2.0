"""
Interactive SEI Network Visualization with Pyvis
=================================================
Creates multiple interactive views for better understanding:
1. Hierarchical view by generation
2. Clustered view by species type
3. Filterable interface
"""

from pyvis.network import Network
import networkx as nx
from collections import defaultdict
from rdkit import Chem

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

def create_interactive_network(
    graphml_path='sei_full_network.graphml',
    output_path='sei_network_interactive.html',
    max_nodes=300
):
    """Create interactive visualization with hierarchical layout."""
    
    print("="*60)
    print("INTERACTIVE SEI NETWORK VISUALIZATION")
    print("="*60)
    print()
    print("Loading network...")
    
    G = nx.read_graphml(graphml_path)
    
    # Separate nodes
    molecule_nodes = [n for n, d in G.nodes(data=True) if d.get('type') != 'reaction']
    reaction_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'reaction']
    
    print(f"  Molecules: {len(molecule_nodes)}")
    print(f"  Reactions: {len(reaction_nodes)}")
    print()
    
    # Create molecule-only network
    print("Creating molecule network...")
    mol_graph = nx.DiGraph()
    
    for mol in molecule_nodes:
        mol_graph.add_node(mol, **G.nodes[mol])
    
    # Add edges through reactions
    for rxn in reaction_nodes:
        predecessors = [n for n in G.predecessors(rxn) if n in molecule_nodes]
        successors = [n for n in G.successors(rxn) if n in molecule_nodes]
        for pred in predecessors:
            for succ in successors:
                if pred != succ:
                    mol_graph.add_edge(pred, succ)
    
    # Filter to important nodes
    if len(mol_graph.nodes()) > max_nodes:
        print(f"Filtering to top {max_nodes} molecules...")
        pagerank = nx.pagerank(mol_graph)
        important_nodes = sorted(
            mol_graph.nodes(),
            key=lambda n: pagerank.get(n, 0),
            reverse=True
        )[:max_nodes]
        mol_graph = mol_graph.subgraph(important_nodes).copy()
    
    print(f"  Final network: {len(mol_graph.nodes())} molecules")
    print()
    
    # Create Pyvis network
    print("Creating interactive visualization...")
    net = Network(
        height='900px',
        width='100%',
        bgcolor='#ffffff',
        font_color='black',
        directed=True,
        notebook=False,
        cdn_resources='remote'
    )
    
    # Configure physics for hierarchical layout
    net.set_options("""
    {
      "physics": {
        "hierarchicalRepulsion": {
          "centralGravity": 0.0,
          "springLength": 150,
          "springConstant": 0.01,
          "nodeDistance": 200,
          "damping": 0.09
        },
        "maxVelocity": 50,
        "solver": "hierarchicalRepulsion",
        "timestep": 0.5,
        "stabilization": {"iterations": 200}
      },
      "layout": {
        "hierarchical": {
          "enabled": true,
          "levelSeparation": 200,
          "nodeSpacing": 150,
          "treeSpacing": 200,
          "blockShifting": true,
          "edgeMinimization": true,
          "parentCentralization": true,
          "direction": "UD",
          "sortMethod": "directed"
        }
      },
      "interaction": {
        "hover": true,
        "tooltipDelay": 100,
        "navigationButtons": true,
        "keyboard": true
      },
      "nodes": {
        "borderWidth": 2,
        "borderWidthSelected": 4,
        "font": {
          "size": 12,
          "face": "Arial"
        }
      },
      "edges": {
        "smooth": {
          "type": "cubicBezier",
          "forceDirection": "vertical"
        },
        "arrows": {
          "to": {
            "enabled": true,
            "scaleFactor": 0.5
          }
        }
      }
    }
    """)
    
    # Add nodes with hierarchical levels
    print("Adding nodes...")
    species_stats = defaultdict(int)
    
    for node in mol_graph.nodes():
        species_type, color = classify_species(node)
        species_stats[species_type] += 1
        
        label = mol_graph.nodes[node].get('label', node)
        gen = int(mol_graph.nodes[node].get('generation', 0))
        
        # Truncate long labels
        if len(label) > 25:
            display_label = label[:22] + '...'
        else:
            display_label = label
        
        # Size by importance
        degree = mol_graph.degree(node)
        size = 15 + degree * 2
        
        # Tooltip with full info
        tooltip = (
            f"<b>{label}</b><br>"
            f"Type: {species_type}<br>"
            f"Generation: {gen}<br>"
            f"Connections: {degree}<br>"
            f"SMILES: {node[:50]}{'...' if len(node) > 50 else ''}"
        )
        
        # Add node with hierarchical level
        net.add_node(
            node,
            label=display_label if gen == 0 or degree > 5 else '',
            title=tooltip,
            color=color,
            size=size,
            level=gen,  # Hierarchical level
            group=species_type,  # For filtering
            borderWidth=2,
            borderWidthSelected=4
        )
    
    # Add edges
    print("Adding edges...")
    for edge in mol_graph.edges():
        net.add_edge(
            edge[0], edge[1],
            color={'color': '#CCCCCC', 'highlight': '#FF0000'},
            width=1,
            smooth={'type': 'cubicBezier'}
        )
    
    # Create custom HTML with controls
    print("Creating interactive controls...")
    
    controls_html = """
    <div style="position: fixed; top: 10px; left: 10px; background: white; 
                padding: 15px; border: 2px solid #333; border-radius: 8px; 
                font-family: Arial; font-size: 13px; z-index: 1000; max-width: 300px;">
        <h3 style="margin-top: 0; color: #333;">SEI Network Explorer</h3>
        
        <div style="margin-bottom: 10px;">
            <b>Environment:</b> Anode, 0.1V vs Li/Li+, 298K<br>
            <b>Species:</b> """ + str(len(mol_graph.nodes())) + """<br>
            <b>Connections:</b> """ + str(len(mol_graph.edges())) + """
        </div>
        
        <hr style="border: 1px solid #ddd;">
        
        <div style="margin-bottom: 10px;">
            <b>Species Types:</b><br>
            <div style="margin-left: 10px; font-size: 11px;">
    """
    
    for stype, count in sorted(species_stats.items(), key=lambda x: -x[1]):
        _, color = classify_species('dummy' if stype != 'reactant' else 'O=C1OCCO1')
        if stype == 'li_metal':
            _, color = classify_species('[Li]')
        elif stype == 'radical':
            _, color = classify_species('[O-]C')
        elif stype == 'complex':
            _, color = classify_species('A.B')
        elif stype == 'lif':
            _, color = classify_species('LiF')
        elif stype == 'carbonate':
            _, color = classify_species('C(=O)O')
        
        controls_html += f'<span style="color: {color};">●</span> {stype.title()}: {count}<br>'
    
    controls_html += """
            </div>
        </div>
        
        <hr style="border: 1px solid #ddd;">
        
        <div style="font-size: 11px; color: #666;">
            <b>Controls:</b><br>
            • Scroll to zoom<br>
            • Drag to pan<br>
            • Click node for details<br>
            • Hover for tooltip
        </div>
    </div>
    """
    
    # Save network
    print(f"Saving to {output_path}...")
    net.save_graph(output_path)
    
    # Add controls to HTML
    with open(output_path, 'r', encoding='utf-8') as f:
        html = f.read()
    
    html = html.replace('</body>', controls_html + '</body>')
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)
    
    print()
    print("="*60)
    print(f"✓ Interactive visualization saved!")
    print(f"  File: {output_path}")
    print(f"  Open in web browser to explore")
    print("="*60)
    print()
    print("Features:")
    print("  • Hierarchical layout by generation (top to bottom)")
    print("  • Color-coded by species type")
    print("  • Interactive zoom and pan")
    print("  • Click nodes for details")
    print("  • Hover for tooltips")
    print()
    
    return net

if __name__ == "__main__":
    create_interactive_network()
