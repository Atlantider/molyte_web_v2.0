"""
SEI Functional Visualization - Vertical Layout
================================================
Focus on final SEI composition with vertical (tall) layout
Shows inorganic vs organic layers, reactivity, and stability
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
        '[Li]': 'Li',
        '[Li][Li]': 'Li₂',
        '[Li+]': 'Li⁺',
        'F[Li]': 'LiF',
        '[F-]': 'F⁻',
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
                return f'Inorg-{atoms}A'
            elif 'organic' in layer:
                return f'Org-{atoms}A'
            else:
                return f'Int-{atoms}A'
    except:
        pass
    
    return smiles[:12] + '...'

def create_sei_functional_view(
    graphml_path='sei_full_network.graphml',
    composition_file='sei_composition_analysis.json',
    output_path='sei_network_functional.html'
):
    """Create SEI functional visualization with vertical layout."""
    
    print("="*70)
    print("SEI FUNCTIONAL VISUALIZATION (Vertical Layout)")
    print("="*70)
    print()
    
    # Load composition analysis
    print("Loading composition analysis...")
    with open(composition_file, 'r') as f:
        composition = json.load(f)
    
    # Load network
    print("Loading network...")
    G = nx.read_graphml(graphml_path)
    
    molecule_nodes = [n for n, d in G.nodes(data=True) if d.get('type') != 'reaction']
    reaction_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'reaction']
    
    # Filter to SEI-relevant species
    print("Filtering to SEI components...")
    
    sei_species = set()
    
    # Include inorganic and organic SEI
    sei_species.update(composition['inorganic_list'])
    sei_species.update(composition['organic_polymer_list'][:100])  # Top 100 polymers
    sei_species.update(composition['organic_small_list'][:50])     # Top 50 small organics
    
    # Include initial reactants
    initial = ['O=C1OCCO1', 'COC(=O)OC', '[Li+]', 'F[P-](F)(F)(F)(F)F']
    for s in initial:
        if s in molecule_nodes:
            sei_species.add(s)
    
    # Include gases
    sei_species.update(composition['gas_list'])
    
    # Add key intermediates (high out-degree)
    species_details = composition['species_details']
    active_intermediates = [
        s for s, info in species_details.items()
        if info['out_degree'] > 10 and info['category'] == 'intermediate'
    ]
    sei_species.update(active_intermediates[:30])
    
    print(f"  Selected {len(sei_species)} SEI-relevant species")
    
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
    
    print(f"  Network: {len(mol_graph.nodes())} species, {len(mol_graph.edges())} edges")
    print()
    
    # Calculate reactivity
    out_deg = dict(mol_graph.out_degree())
    in_deg = dict(mol_graph.in_degree())
    
    reactivity = {}
    for node in mol_graph.nodes():
        out_d = out_deg.get(node, 0)
        in_d = in_deg.get(node, 0)
        reactivity[node] = out_d / (in_d + 1) if in_d > 0 else out_d
    
    # Create Pyvis network - VERTICAL LAYOUT
    print("Creating visualization...")
    net = Network(
        height='1200px',  # Tall!
        width='100%',
        bgcolor='#FAFAFA',
        font_color='black',
        directed=True,
        notebook=False
    )
    
    # CRITICAL: Use LEFT-RIGHT direction for vertical appearance
    net.set_options("""
    {
      "layout": {
        "hierarchical": {
          "enabled": true,
          "levelSeparation": 200,
          "nodeSpacing": 150,
          "treeSpacing": 180,
          "direction": "LR",
          "sortMethod": "directed"
        }
      },
      "physics": {"enabled": false},
      "interaction": {
        "hover": true,
        "navigationButtons": true,
        "keyboard": true
      },
      "nodes": {"font": {"size": 12, "face": "Arial", "bold": true}},
      "edges": {
        "smooth": {"type": "cubicBezier"},
        "arrows": {"to": {"enabled": true, "scaleFactor": 0.7}},
        "width": 1.5
      }
    }
    """)
    
    # Add nodes
    layer_stats = Counter()
    
    for node in mol_graph.nodes():
        layer, color, base_size = get_sei_layer(node)
        layer_stats[layer] += 1
        
        label = get_chemical_name_sei(node)
        gen = int(mol_graph.nodes[node].get('generation', 0))
        
        # Node size by out-degree (reactivity)
        out_d = out_deg.get(node, 0)
        size = base_size + out_d * 2
        
        # Border width by reactivity
        react = reactivity.get(node, 0)
        if react > 2.0:
            border_width = 6  # Very reactive
        elif react < 0.5:
            border_width = 1  # Very stable
        else:
            border_width = 3
        
        # Tooltip
        tooltip = (
            f"<b>{label}</b><br>"
            f"Layer: {layer}<br>"
            f"Generation: {gen}<br>"
            f"出度: {out_d} (引入反应)<br>"
            f"入度: {in_deg.get(node, 0)}<br>"
            f"反应性: {react:.2f}<br>"
            f"状态: {'活性中间体' if react > 2.0 else '稳定产物' if react < 0.5 else '参与反应'}<br>"
            f"<br>SMILES: {node[:60]}"
        )
        
        net.add_node(
            node,
            label=label,
            title=tooltip,
            color=color,
            size=size,
            level=gen,
            borderWidth=border_width,
            shape='dot'
        )
    
    # Add edges
    for edge in mol_graph.edges():
        net.add_edge(
            edge[0], edge[1],
            width=1.5,
            color={'color': '#999999'}
        )
    
    # Comprehensive legend
    legend_html = f"""
    <div style="position: fixed; top: 10px; right: 10px; background: white; 
                padding: 20px; border: 3px solid #333; border-radius: 10px; 
                font-family: Arial; font-size: 13px; z-index: 1000; max-width: 380px;
                box-shadow: 0 6px 12px rgba(0,0,0,0.15); max-height: 95vh; overflow-y: auto;">
        
        <h2 style="margin: 0 0 15px 0; color: #333; border-bottom: 3px solid #333; padding-bottom: 10px;">
            🔋 SEI 功能分析
        </h2>
        
        <!-- SEI层次 -->
        <div style="background: #E3F2FD; padding: 12px; border-radius: 6px; margin-bottom: 15px; border-left: 4px solid #2196F3;">
            <b style="color: #1565C0;">📊 SEI 层次组成:</b><br>
            <div style="margin-top: 8px; font-size: 12px; line-height: 2.0;">
                <span style="display: inline-block; width: 14px; height: 14px; 
                             background: #3366FF; border-radius: 50%; margin-right: 6px;"></span>
                <b>无机底层</b>: {layer_stats['inorganic']} 种<br>
                <span style="font-size: 10px; margin-left: 20px; color: #666;">
                LiF, Li₂O, Li₂CO₃ 等</span><br>
                
                <span style="display: inline-block; width: 14px; height: 14px; 
                             background: #33CC33; border-radius: 50%; margin-right: 6px;"></span>
                <b>有机层-聚合物</b>: {layer_stats['organic_polymer']} 种<br>
                <span style="font-size: 10px; margin-left: 20px; color: #666;">
                \u003e15 原子的聚合体</span><br>
                
                <span style="display: inline-block; width: 14px; height: 14px; 
                             background: #66FF66; border-radius: 50%; margin-right: 6px;"></span>
                <b>有机层-小分子</b>: {layer_stats['organic_small']} 种<br>
                <span style="font-size: 10px; margin-left: 20px; color: #666;">
                烷基碳酸锂等</span><br>
                
                <span style="display: inline-block; width: 14px; height: 14px; 
                             background: #AAAAAA; border-radius: 50%; margin-right: 6px;"></span>
                <b>气体产物</b>: {layer_stats['gas']} 种<br>
                <span style="font-size: 10px; margin-left: 20px; color: #666;">
                CO₂, CO, C₂H₄</span><br>
            </div>
        </div>
        
        <!-- 反应性 -->
        <div style="background: #FFF3E0; padding: 12px; border-radius: 6px; margin-bottom: 15px; border-left: 4px solid #FF9800;">
            <b style="color: #E65100;">⚡ 反应性指标:</b><br>
            <div style="font-size: 12px; margin-top: 8px; line-height: 1.6;">
            
            <b>节点大小 = 出度</b><br>
            <span style="font-size: 11px;">• 大节点 = 引入更多反应<br>
            • 小节点 = 反应少或稳定</span><br><br>
            
            <b>边框粗细 = 反应性</b><br>
            <span style="font-size: 11px;">
            • 粗边框 (6px) = 非常活跃 (反应性\u003e2)<br>
            • 细边框 (1px) = 很稳定 (反应性\u003c0.5)<br>
            • 中等 (3px) = 参与反应
            </span>
            </div>
        </div>
        
        <!-- 布局 -->
        <div style="background: #E8F5E9; padding: 12px; border-radius: 6px; margin-bottom: 15px; border-left: 4px solid #4CAF50;">
            <b style="color: #2E7D32;">📐 纵向布局:</b><br>
            <div style="font-size: 12px; margin-top: 8px;">
            • <b>左侧</b>: 初始反应物<br>
            • <b>中间</b>: 中间体/活性物种<br>
            • <b>右侧</b>: 最终 SEI 产物<br>
            • <b>从左到右</b>: 反应进程
            </div>
        </div>
        
        <!-- 关键指标 -->
        <div style="background: #FCE4EC; padding: 12px; border-radius: 6px; margin-bottom: 15px;">
            <b style="color: #C2185B;">💡 如何解读:</b><br>
            <div style="font-size: 11px; margin-top: 8px; line-height: 1.5;">
            1. <b>无机物</b> → 底层，提供机械强度<br>
            2. <b>有机聚合物</b> → 外层，柔韧性<br>
            3. <b>高出度</b> → 关键中间体<br>
            4. <b>低出度+高入度</b> → 最终产物<br>
            5. <b>气体</b> → 不在 SEI 中
            </div>
        </div>
        
        <!-- 统计 -->
        <div style="background: #F5F5F5; padding: 10px; border-radius: 6px; font-size: 11px;">
            <b>网络统计:</b><br>
            • 物种: {len(mol_graph.nodes())}<br>
            • 连接: {len(mol_graph.edges())}<br>
            • 布局: 纵向 (左→右)
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
    print("✓ SEI functional visualization created!")
    print(f"  File: {output_path}")
    print("="*70)
    print()
    print("Layout: VERTICAL (tall, left→right = reaction progress)")
    print(f"Inorganic SEI: {layer_stats['inorganic']} species")
    print(f"Organic SEI: {layer_stats['organic_polymer'] + layer_stats['organic_small']} species")
    print()

if __name__ == "__main__":
    create_sei_functional_view()
