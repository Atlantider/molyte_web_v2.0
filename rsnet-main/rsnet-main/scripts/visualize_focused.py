"""
Focused SEI Network Visualization
==================================
Creates a simplified view highlighting only the most important reaction pathways.
"""

from pyvis.network import Network
import networkx as nx
from collections import defaultdict

def classify_species(smiles):
    """Classify species by chemical type."""
    if smiles in ['O=C1OCCO1', 'COC(=O)OC', '[Li+]', 'F[P-](F)(F)(F)(F)F']:
        return 'reactant', '#FF4444', 3
    if smiles in ['[Li]', '[Li][Li]']:
        return 'li_metal', '#00CED1', 2
    if 'Li' in smiles and 'F' in smiles and '.' not in smiles:
        return 'lif', '#FF69B4', 1  # High priority SEI product
    if ('C(=O)O' in smiles or 'CO3' in smiles) and '.' not in smiles:
        return 'carbonate', '#9370DB', 1  # High priority SEI product
    if '[O-]' in smiles or ('[C]' in smiles and '(' not in smiles.split('[C]')[1][:3] if '[C]' in smiles else False):
        return 'radical', '#FFD700', 2
    if '.' in smiles:
        return 'complex', '#87CEEB', 0  # Low priority
    return 'intermediate', '#90EE90', 1

def create_focused_visualization(
    graphml_path='sei_full_network.graphml',
    output_path='sei_network_focused.html',
    max_nodes=100
):
    """Create focused visualization highlighting key pathways."""
    
    print("="*60)
    print("FOCUSED SEI NETWORK VISUALIZATION")
    print("="*60)
    print()
    
    # Load network
    print("Loading network...")
    G = nx.read_graphml(graphml_path)
    
    molecule_nodes = [n for n, d in G.nodes(data=True) if d.get('type') != 'reaction']
    reaction_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'reaction']
    
    # Create molecule network
    print("Creating molecule network...")
    mol_graph = nx.DiGraph()
    
    for mol in molecule_nodes:
        mol_graph.add_node(mol, **G.nodes[mol])
    
    for rxn in reaction_nodes:
        preds = [n for n in G.predecessors(rxn) if n in molecule_nodes]
        succs = [n for n in G.successors(rxn) if n in molecule_nodes]
        for p in preds:
            for s in succs:
                if p != s:
                    mol_graph.add_edge(p, s)
    
    # Filter by priority
    print("Filtering to key species...")
    
    # Score each node
    node_scores = {}
    for node in mol_graph.nodes():
        _, _, priority = classify_species(node)
        degree = mol_graph.degree(node)
        gen = int(mol_graph.nodes[node].get('generation', 0))
        
        # Score = priority * degree - generation_penalty
        score = priority * degree - gen * 0.5
        node_scores[node] = score
    
    # Select top nodes
    important_nodes = sorted(node_scores.keys(), key=lambda n: node_scores[n], reverse=True)[:max_nodes]
    mol_graph = mol_graph.subgraph(important_nodes).copy()
    
    print(f"  Selected {len(mol_graph.nodes())} key species")
    print()
    
    # Create Pyvis network
    print("Creating interactive visualization...")
    net = Network(
        height='900px',
        width='100%',
        bgcolor='#ffffff',
        font_color='black',
        directed=True,
        notebook=False
    )
    
    # Hierarchical layout
    net.set_options("""
    {
      "layout": {
        "hierarchical": {
          "enabled": true,
          "levelSeparation": 250,
          "nodeSpacing": 200,
          "treeSpacing": 250,
          "direction": "UD",
          "sortMethod": "directed"
        }
      },
      "physics": {
        "enabled": false
      },
      "interaction": {
        "hover": true,
        "navigationButtons": true,
        "keyboard": true
      },
      "nodes": {
        "font": {"size": 14, "face": "Arial", "bold": true}
      },
      "edges": {
        "smooth": {"type": "cubicBezier"},
        "arrows": {"to": {"enabled": true, "scaleFactor": 0.8}},
        "width": 2
      }
    }
    """)
    
    # Add nodes
    print("Adding nodes...")
    species_stats = defaultdict(int)
    
    for node in mol_graph.nodes():
        species_type, color, priority = classify_species(node)
        species_stats[species_type] += 1
        
        label = mol_graph.nodes[node].get('label', node)
        gen = int(mol_graph.nodes[node].get('generation', 0))
        degree = mol_graph.degree(node)
        
        # Truncate label
        if len(label) > 20:
            display_label = label[:17] + '...'
        else:
            display_label = label
        
        # Size by importance
        size = 20 + priority * 10 + degree * 3
        
        # Tooltip
        tooltip = (
            f"<b>{label}</b><br>"
            f"Type: {species_type}<br>"
            f"Generation: {gen}<br>"
            f"Connections: {degree}<br>"
            f"Priority: {'★' * priority}<br>"
            f"SMILES: {node[:60]}{'...' if len(node) > 60 else ''}"
        )
        
        net.add_node(
            node,
            label=display_label,
            title=tooltip,
            color=color,
            size=size,
            level=gen,
            borderWidth=3
        )
    
    # Add edges
    print("Adding edges...")
    for edge in mol_graph.edges():
        net.add_edge(
            edge[0], edge[1],
            color={'color': '#999999', 'highlight': '#FF0000'},
            width=2
        )
    
    # Custom HTML
    controls_html = f"""
    <div style="position: fixed; top: 10px; left: 10px; background: white; 
                padding: 20px; border: 3px solid #333; border-radius: 10px; 
                font-family: Arial; font-size: 14px; z-index: 1000; max-width: 350px;
                box-shadow: 0 4px 6px rgba(0,0,0,0.1);">
        <h2 style="margin-top: 0; color: #333; border-bottom: 2px solid #333; padding-bottom: 10px;">
            🔬 SEI Network (Focused View)
        </h2>
        
        <div style="margin-bottom: 15px; background: #f0f0f0; padding: 10px; border-radius: 5px;">
            <b>Environment:</b> Anode, 0.1V, 298K<br>
            <b>Species:</b> {len(mol_graph.nodes())} (key pathways only)<br>
            <b>Connections:</b> {len(mol_graph.edges())}
        </div>
        
        <div style="margin-bottom: 15px;">
            <b>📊 Species Distribution:</b><br>
            <div style="margin-left: 10px; font-size: 12px; margin-top: 5px;">
    """
    
    for stype, count in sorted(species_stats.items(), key=lambda x: -x[1]):
        _, color, _ = classify_species('dummy' if stype != 'reactant' else 'O=C1OCCO1')
        if stype == 'li_metal': _, color, _ = classify_species('[Li]')
        elif stype == 'radical': _, color, _ = classify_species('[O-]C')
        elif stype == 'lif': _, color, _ = classify_species('LiF')
        elif stype == 'carbonate': _, color, _ = classify_species('C(=O)O')
        
        controls_html += f'<span style="color: {color}; font-size: 16px;">●</span> <b>{stype.title()}</b>: {count}<br>'
    
    controls_html += """
            </div>
        </div>
        
        <div style="background: #fff3cd; padding: 10px; border-radius: 5px; border-left: 4px solid #ffc107;">
            <b>💡 层次化布局说明:</b><br>
            <div style="font-size: 12px; margin-top: 5px;">
            • <b>从上到下</b> = 从初始反应物到产物<br>
            • <b>同一层</b> = 同一代数生成<br>
            • <b>箭头方向</b> = 反应路径<br>
            • <b>节点大小</b> = 重要性(连接数×优先级)
            </div>
        </div>
        
        <div style="margin-top: 15px; font-size: 12px; color: #666; border-top: 1px solid #ddd; padding-top: 10px;">
            <b>🖱️ 操作:</b> 滚轮缩放 | 拖动平移 | 点击查看详情
        </div>
    </div>
    """
    
    # Save
    print(f"Saving to {output_path}...")
    net.save_graph(output_path)
    
    with open(output_path, 'r', encoding='utf-8') as f:
        html = f.read()
    
    html = html.replace('</body>', controls_html + '</body>')
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)
    
    print()
    print("="*60)
    print("✓ Focused visualization created!")
    print(f"  File: {output_path}")
    print("="*60)
    print()
    print("Key improvements:")
    print("  • Only 100 most important species (vs 300 before)")
    print("  • Larger labels for better readability")
    print("  • Clearer hierarchical layout (generation-based)")
    print("  • Priority-based filtering (SEI products prioritized)")
    print()

if __name__ == "__main__":
    create_focused_visualization()
