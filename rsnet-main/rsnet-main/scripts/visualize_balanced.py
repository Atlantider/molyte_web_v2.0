"""
SEI Functional Visualization - Balanced Layout
================================================
均衡的布局设计，使用力导向布局而非层次化，实现方形显示
"""

from pyvis.network import Network
import networkx as nx
import json
from rdkit import Chem
from collections import Counter

def get_sei_layer(smiles):
    """Determine which SEI layer species belongs to."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 'other', '#CCCCCC', 20
        
        atoms = Counter([atom.GetSymbol() for atom in mol.GetAtoms()])
        num_heavy = mol.GetNumHeavyAtoms()
        
        # Inorganic bottom layer
        if atoms.get('Li', 0) > 0 and atoms.get('C', 0) == 0:
            return 'inorganic', '#3366FF', 45  # Blue
        
        # Organic outer layer - polymers
        if atoms.get('Li', 0) > 0 and atoms.get('C', 0) > 0 and num_heavy > 15:
            return 'organic_polymer', '#33CC33', 40  # Green
        
        # Organic outer layer - small
        if atoms.get('Li', 0) > 0 and atoms.get('C', 0) > 0:
            return 'organic_small', '#66FF66', 35  # Light green
        
        # Gases
        if smiles in ['O=C=O', '[C-]#[O+]', 'C=C']:
            return 'gas', '#AAAAAA', 30  # Gray
        
        return 'intermediate', '#FFD700', 25  # Gold
    except:
        return 'other', '#CCCCCC', 20

def get_chemical_name_sei(smiles):
    """Get chemical name for SEI components."""
    names = {
        'O=C1OCCO1': 'EC',
        'COC(=O)OC': 'DMC',
        '[Li]': 'Li',
        '[Li][Li]': 'Li₂',
        '[Li+]': 'Li⁺',
        'F[Li]': 'LiF',
        '[F-]': 'F⁻',
        'F[P-](F)(F)(F)(F)F': 'PF₆⁻',
        'O=C=O': 'CO₂',
        '[C-]#[O+]': 'CO',
        'C=C': 'C₂H₄',
    }
    
    if smiles in names:
        return names[smiles]
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            layer, _, _ = get_sei_layer(smiles)
            atoms = mol.GetNumHeavyAtoms()
            
            if layer == 'inorganic':
                return f'无机-{atoms}A'
            elif 'organic' in layer:
                return f'有机-{atoms}A'
            else:
                return f'中间体'
    except:
        pass
    
    return smiles[:8] + '...'

def create_balanced_sei_view(
    graphml_path='sei_full_network.graphml',
    composition_file='sei_composition_analysis.json',
    output_path='sei_network_balanced.html'
):
    """Create balanced layout SEI visualization."""
    
    print("="*70)
    print("SEI BALANCED LAYOUT VISUALIZATION")
    print("="*70)
    print()
    
    # Load composition
    with open(composition_file, 'r') as f:
        composition = json.load(f)
    
    # Load network
    G = nx.read_graphml(graphml_path)
    
    molecule_nodes = [n for n, d in G.nodes(data=True) if d.get('type') != 'reaction']
    reaction_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'reaction']
    
    # Filter to key SEI species (reduced for clarity)
    sei_species = set()
    
    # Initial reactants
    initial = ['O=C1OCCO1', 'COC(=O)OC', '[Li+]', 'F[P-](F)(F)(F)(F)F']
    for s in initial:
        if s in molecule_nodes:
            sei_species.add(s)
    
    # Top inorganic
    sei_species.update(composition['inorganic_list'][:20])
    
    # Top organic
    sei_species.update(composition['organic_polymer_list'][:40])
    sei_species.update(composition['organic_small_list'][:30])
    
    # Gases
    sei_species.update(composition['gas_list'])
    
    # Key intermediates (by connectivity)
    species_details = composition['species_details']
    active = [
        s for s, info in species_details.items()
        if info['out_degree'] > 5
    ]
    sei_species.update(active[:20])
    
    print(f"Selected {len(sei_species)} species for balanced view")
    
    # Create subgraph
    mol_graph = nx.DiGraph()
    for mol in sei_species:
        if mol in molecule_nodes:
            mol_graph.add_node(mol, **G.nodes[mol])
    
    # Add edges
    for rxn in reaction_nodes:
        preds = [n for n in G.predecessors(rxn) if n in sei_species]
        succs = [n for n in G.successors(rxn) if n in sei_species]
        
        for p in preds:
            for s in succs:
                if p != s:
                    mol_graph.add_edge(p, s)
    
    print(f"Network: {len(mol_graph.nodes())} species, {len(mol_graph.edges())} edges")
    print()
    
    # Calculate metrics
    out_deg = dict(mol_graph.out_degree())
    in_deg = dict(mol_graph.in_degree())
    
    # Create Pyvis network with BALANCED force-directed layout
    net = Network(
        height='900px',
        width='1400px',  # Balanced aspect ratio
        bgcolor='#FAFAFA',
        font_color='black',
        directed=True,
        notebook=False
    )
    
    # Use physics-based layout for balance
    net.set_options("""
    {
      "physics": {
        "enabled": true,
        "barnesHut": {
          "gravitationalConstant": -15000,
          "centralGravity": 0.3,
          "springLength": 180,
          "springConstant": 0.04,
          "damping": 0.09,
          "avoidOverlap": 0.2
        },
        "maxVelocity": 50,
        "minVelocity": 0.1,
        "solver": "barnesHut",
        "stabilization": {
          "enabled": true,
          "iterations": 300,
          "updateInterval": 25
        }
      },
      "interaction": {
        "hover": true,
        "navigationButtons": true,
        "keyboard": true,
        "tooltipDelay": 100
      },
      "nodes": {
        "font": {"size": 14, "face": "Arial", "bold": true}
      },
      "edges": {
        "smooth": {
          "type": "continuous",
          "roundness": 0.5
        },
        "arrows": {
          "to": {"enabled": true, "scaleFactor": 0.8}
        },
        "width": 2
      }
    }
    """)
    
    # Add nodes with grouping for layout
    layer_stats = Counter()
    
    for node in mol_graph.nodes():
        layer, color, base_size = get_sei_layer(node)
        layer_stats[layer] += 1
        
        label = get_chemical_name_sei(node)
        
        # Size by out-degree
        out_d = out_deg.get(node, 0)
        size = base_size + out_d * 3
        
        # Border by stability
        reactivity = out_d / (in_deg.get(node, 0) + 1)
        if reactivity > 2.0:
            border_width = 5
            border_color = '#FF6B6B'  # Red border for high reactivity
        elif reactivity < 0.5:
            border_width = 2
            border_color = '#4CAF50'  # Green border for stable
        else:
            border_width = 3
            border_color = color
        
        # Tooltip
        tooltip = (
            f"<b>{label}</b><br>"
            f"SEI层: {layer}<br>"
            f"出度: {out_d} →<br>"
            f"入度: {in_deg.get(node, 0)} ←<br>"
            f"{'🔴 高活性' if reactivity > 2.0 else '🟢 稳定' if reactivity < 0.5 else '⚪ 参与反应'}<br>"
            f"<br>SMILES: {node[:70]}"
        )
        
        # Group for physics
        group = {
            'inorganic': 0,
            'organic_polymer': 1,
            'organic_small': 2,
            'gas': 3,
            'intermediate': 4,
            'other': 5
        }.get(layer, 5)
        
        net.add_node(
            node,
            label=label,
            title=tooltip,
            color={'background': color, 'border': border_color},
            size=size,
            borderWidth=border_width,
            shape='dot',
            group=group
        )
    
    # Add edges with reactivity coloring
    for edge in mol_graph.edges():
        source_out = out_deg.get(edge[0], 0)
        
        # Edge color based on source reactivity
        if source_out > 10:
            edge_color = '#FF6B6B'  # Red for highly reactive source
            width = 3
        elif source_out > 5:
            edge_color = '#FFA726'  # Orange
            width = 2
        else:
            edge_color = '#999999'  # Gray
            width = 1.5
        
        net.add_edge(
            edge[0], edge[1],
            width=width,
            color={'color': edge_color}
        )
    
    # Enhanced legend
    legend_html = f"""
    <div style="position: fixed; top: 10px; right: 10px; background: white; 
                padding: 18px; border: 3px solid #333; border-radius: 12px; 
                font-family: Arial; font-size: 13px; z-index: 1000; max-width: 350px;
                box-shadow: 0 6px 12px rgba(0,0,0,0.2);">
        
        <h2 style="margin: 0 0 12px 0; color: #333; border-bottom: 3px solid #333; padding-bottom: 8px; font-size: 16px;">
            🔋 SEI 成分分析
        </h2>
        
        <div style="background: #E3F2FD; padding: 10px; border-radius: 6px; margin-bottom: 12px;">
            <b style="color: #1565C0;">📊 SEI 层次:</b><br>
            <div style="margin-top: 6px; font-size: 12px; line-height: 1.8;">
                <span style="color: #3366FF; font-size: 16px;">●</span> 
                <b>无机底层</b>: {layer_stats['inorganic']}<br>
                
                <span style="color: #33CC33; font-size: 16px;">●</span> 
                <b>有机聚合物</b>: {layer_stats['organic_polymer']}<br>
                
                <span style="color: #66FF66; font-size: 16px;">●</span> 
                <b>有机小分子</b>: {layer_stats['organic_small']}<br>
                
                <span style="color: #AAAAAA; font-size: 16px;">●</span> 
                <b>气体</b>: {layer_stats['gas']}<br>
            </div>
        </div>
        
        <div style="background: #FFF3E0; padding: 10px; border-radius: 6px; margin-bottom: 12px;">
            <b style="color: #E65100;">⚡ 反应性编码:</b><br>
            <div style="font-size: 11px; margin-top: 6px; line-height: 1.5;">
            
            <b>节点大小</b> = 出度 (→箭头数量)<br>
            <span style="font-size: 10px;">大 = 引入更多反应</span><br><br>
            
            <b>节点边框</b>:<br>
            <span style="color: #FF6B6B;">● 粗红框</span> = 高活性<br>
            <span style="color: #4CAF50;">● 细绿框</span> = 稳定产物<br>
            <span>● 普通框</span> = 参与反应<br><br>
            
            <b>箭头颜色/粗细</b>:<br>
            <span style="color: #FF6B6B;">━━</span> 红粗 = 高活性源<br>
            <span style="color: #FFA726;">━━</span> 橙中 = 中活性<br>
            <span style="color: #999999;">━━</span> 灰细 = 低活性
            </div>
        </div>
        
        <div style="background: #E8F5E9; padding: 10px; border-radius: 6px; margin-bottom: 12px;">
            <b style="color: #2E7D32;">💡 关键指标:</b><br>
            <div style="font-size: 11px; margin-top: 6px; line-height: 1.4;">
            • <b>无机:</b> 底层，机械强度<br>
            • <b>有机聚合物:</b> 外层，柔韧<br>
            • <b>高出度:</b> 引入更多反应<br>
            • <b>低出度:</b> 稳定 SEI 组分
            </div>
        </div>
        
        <div style="background: #F5F5F5; padding: 8px; border-radius: 6px; font-size: 11px;">
            <b>布局:</b> 力导向 (均衡)<br>
            <b>物种:</b> {len(mol_graph.nodes())}<br>
            <b>连接:</b> {len(mol_graph.edges())}
        </div>
        
    </div>
    """
    
    # Save
    print(f"Saving to {output_path}...")
    net.save_graph(output_path)
    
    with open(output_path, 'r', encoding='utf-8') as f:
        html = f.read()
    
    html = html.replace('</body>', legend_html + '</body>')
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)
    
    print()
    print("="*70)
    print("✓ Balanced SEI visualization created!")
    print(f"  File: {output_path}")
    print("="*70)
    print()
    print("Layout: BALANCED force-directed (900×1400px)")
    print(f"Inorganic: {layer_stats['inorganic']} | Organic: {layer_stats['organic_polymer'] + layer_stats['organic_small']}")
    print()

if __name__ == "__main__":
    create_balanced_sei_view()
