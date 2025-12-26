"""
Enhanced SEI Visualization - Reaction-Focused Tooltips
========================================================
Tooltip显示具体反应列表
"""

from pyvis.network import Network
import networkx as nx
import json
from rdkit import Chem
from collections import Counter, defaultdict

def get_sei_layer_enhanced(smiles):
    """分类并返回层次级别."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 'other', '#CCCCCC', 20, 4
        
        atoms = Counter([atom.GetSymbol() for atom in mol.GetAtoms()])
        num_heavy = mol.GetNumHeavyAtoms()
        
        if smiles in ['O=C1OCCO1', 'COC(=O)OC', '[Li+]', 'F[P-](F)(F)(F)(F)F']:
            return 'initial', '#FF6B6B', 50, 0
        
        if atoms.get('Li', 0) > 0 and atoms.get('C', 0) == 0:
            return 'inorganic_sei', '#3366FF', 45, 3
        
        if atoms.get('Li', 0) > 0 and atoms.get('C', 0) > 0 and num_heavy > 15:
            return 'organic_polymer', '#33CC33', 40, 3
        
        if atoms.get('Li', 0) > 0 and atoms.get('C', 0) > 0:
            return 'organic_small', '#66FF66', 35, 3
        
        if smiles in ['O=C=O', '[C-]#[O+]', 'C=C']:
            return 'gas', '#AAAAAA', 30, 2
        
        return 'intermediate', '#FFD700', 25, 1
    except:
        return 'other', '#CCCCCC', 20, 4

def get_chemical_name_sei(smiles):
    """化学名称."""
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
            layer, _, _, _ = get_sei_layer_enhanced(smiles)
            atoms = mol.GetNumHeavyAtoms()
            
            if 'sei' in layer:
                return f'SEI-{atoms}A'
            elif layer == 'gas':
                return f'Gas'
            else:
                return f'{atoms}A'
    except:
        pass
    
    return smiles[:8] + '...'

def create_enhanced_sei_view(
    graphml_path='sei_full_network.graphml',
    composition_file='sei_composition_analysis.json',
    output_path='sei_network_enhanced.html'
):
    """创建增强版SEI可视化."""
    
    print("="*70)
    print("ENHANCED SEI VISUALIZATION")
    print("="*70)
    print()
    
    # Load data
    with open(composition_file, 'r') as f:
        composition = json.load(f)
    
    G = nx.read_graphml(graphml_path)
    
    molecule_nodes = [n for n, d in G.nodes(data=True) if d.get('type') != 'reaction']
    reaction_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'reaction']
    
    # 筛选关键物种
    sei_species = set()
    
    # 初始
    initial = ['O=C1OCCO1', 'COC(=O)OC', '[Li+]', 'F[P-](F)(F)(F)(F)F']
    sei_species.update(s for s in initial if s in molecule_nodes)
    
    # SEI组分
    sei_species.update(composition['inorganic_list'][:25])
    sei_species.update(composition['organic_polymer_list'][:45])
    sei_species.update(composition['organic_small_list'][:35])
    sei_species.update(composition['gas_list'])
    
    # 关键中间体
    species_details = composition['species_details']
    active = [s for s, info in species_details.items() if info['out_degree'] > 5]
    sei_species.update(active[:25])
    
    # 构建子图
    mol_graph = nx.DiGraph()
    for mol in sei_species:
        if mol in molecule_nodes:
            mol_graph.add_node(mol, **G.nodes[mol])
    
    # 构建反应映射 (molecule -> reactions)
    print("Building reaction mappings...")
    incoming_reactions = defaultdict(list)  # 生成该节点的反应
    outgoing_reactions = defaultdict(list)  # 消耗该节点的反应
    
    # 初始反应物列表 (这些不应该有incoming reactions)
    initial_reactants = {'O=C1OCCO1', 'COC(=O)OC', '[Li+]', 'F[P-](F)(F)(F)(F)F'}
    
    for rxn in reaction_nodes:
        rxn_label = G.nodes[rxn].get('label', rxn)
        
        # 过滤掉coordination/clustering反应 (solvation equilibria)
        # 这些反应只是物理结合/分离,不是真正的化学反应
        if any(keyword in rxn_label.lower() for keyword in ['clustering', 'coordination']):
            continue  # 跳过溶剂化反应
        
        preds = [n for n in G.predecessors(rxn) if n in sei_species]
        succs = [n for n in G.successors(rxn) if n in sei_species]
        
        # 对于产物,这是incoming reaction
        for succ in succs:
            # 初始反应物不应该有incoming reactions
            if succ not in initial_reactants:
                incoming_reactions[succ].append(rxn_label)
        
        # 对于反应物,这是outgoing reaction  
        for pred in preds:
            outgoing_reactions[pred].append(rxn_label)
        
        # 添加边
        for p in preds:
            for s in succs:
                if p != s:
                    mol_graph.add_edge(p, s)
    
    print(f"Network: {len(mol_graph.nodes())} species, {len(mol_graph.edges())} edges")
    
    # 计算度数
    out_deg = dict(mol_graph.out_degree())
    in_deg = dict(mol_graph.in_degree())
    
    # 创建Pyvis网络
    net = Network(
        height='900px',
        width='1400px',
        bgcolor='#FAFAFA',
        font_color='black',
        directed=True,
        notebook=False
    )
    
    # 使用均衡的力导向布局
    net.set_options("""
    {
      "physics": {
        "enabled": true,
        "barnesHut": {
          "gravitationalConstant": -20000,
          "centralGravity": 0.3,
          "springLength": 200,
          "springConstant": 0.04,
          "damping": 0.09,
          "avoidOverlap": 0.3
        },
        "maxVelocity": 50,
        "minVelocity": 0.1,
        "solver": "barnesHut",
        "stabilization": {
          "enabled": true,
          "iterations": 400,
          "updateInterval": 25
        }
      },
      "interaction": {
        "hover": true,
        "navigationButtons": true,
        "keyboard": true,
        "tooltipDelay": 100,
        "zoomView": true,
        "dragView": true
      },
      "nodes": {
        "font": {"size": 13, "face": "Arial", "bold": true}
      },
      "edges": {
        "smooth": {"type": "continuous", "roundness": 0.5},
        "arrows": {"to": {"enabled": true, "scaleFactor": 0.8}},
        "width": 2
      }
    }
    """)
    
    # 添加节点 (with reaction lists in tooltips)
    layer_stats = Counter()
    
    for node in mol_graph.nodes():
        layer, color, base_size, level = get_sei_layer_enhanced(node)
        layer_stats[layer] += 1
        
        label = get_chemical_name_sei(node)
        
        out_d = out_deg.get(node, 0)
        in_d = in_deg.get(node, 0)
        size = base_size + out_d * 3
        
        reactivity = out_d / (in_d + 1) if in_d > 0 else out_d
        
        if reactivity > 2.0:
            border_width = 5
            border_color = '#FF6B6B'
            stability = '🔴 高活性'
        elif reactivity < 0.5:
            border_width = 2
            border_color = '#4CAF50'
            stability = '🟢 稳定产物'
        else:
            border_width = 3
            border_color = color
            stability = '⚪ 参与反应'
        
        # 计算连接节点的类型统计
        predecessors = list(mol_graph.predecessors(node))
        successors = list(mol_graph.successors(node))
        
        # Layer name mapping
        layer_name_cn = {
            'initial': '初始反应物',
            'inorganic_sei': '无机SEI(底层)',
            'organic_polymer': '有机SEI(聚合物)',
            'organic_small': '有机SEI(小分子)',
            'gas': '气体产物',
            'intermediate': '中间体'
        }.get(layer, '其他')
        
        # 统计前驱节点的类型
        pred_stats = Counter()
        for pred in predecessors:
            pred_layer, _, _, _ = get_sei_layer_enhanced(pred)
            pred_stats[pred_layer] += 1
        
        # 统计后继节点的类型
        succ_stats = Counter()
        for succ in successors:
            succ_layer, _, _, _ = get_sei_layer_enhanced(succ)
            succ_stats[succ_layer] += 1
        
        # 按优先级组织tooltip
        tooltip_lines = [
            f"【{label}】",
            "",
            f"类型: {layer_name_cn}",
            f"状态: {stability}",
            "",
        ]
        
        # 前驱统计
        if predecessors:
            tooltip_lines.append(f"⬅️ 前驱物种 (来自{len(predecessors)}个):")
            if pred_stats.get('inorganic_sei', 0):
                tooltip_lines.append(f"  无机SEI: {pred_stats['inorganic_sei']}个")
            if pred_stats.get('organic_polymer', 0):
                tooltip_lines.append(f"  有机聚合物: {pred_stats['organic_polymer']}个")
            if pred_stats.get('organic_small', 0):
                tooltip_lines.append(f"  有机小分子: {pred_stats['organic_small']}个")
            if pred_stats.get('intermediate', 0):
                tooltip_lines.append(f"  中间体: {pred_stats['intermediate']}个")
            if pred_stats.get('initial', 0):
                tooltip_lines.append(f"  初始反应物: {pred_stats['initial']}个")
            if pred_stats.get('gas', 0):
                tooltip_lines.append(f"  气体: {pred_stats['gas']}个")
        else:
            tooltip_lines.append("⬅️ 无前驱 (初始物种)")
        
        tooltip_lines.append("")
        
        # 后继统计
        if successors:
            tooltip_lines.append(f"➡️ 后继物种 (到{len(successors)}个):")
            if succ_stats.get('inorganic_sei', 0):
                tooltip_lines.append(f"  无机SEI: {succ_stats['inorganic_sei']}个")
            if succ_stats.get('organic_polymer', 0):
                tooltip_lines.append(f"  有机聚合物: {succ_stats['organic_polymer']}个")
            if succ_stats.get('organic_small', 0):
                tooltip_lines.append(f"  有机小分子: {succ_stats['organic_small']}个")
            if succ_stats.get('intermediate', 0):
                tooltip_lines.append(f"  中间体: {succ_stats['intermediate']}个")
            if succ_stats.get('gas', 0):
                tooltip_lines.append(f"  气体: {succ_stats['gas']}个")
        else:
            tooltip_lines.append("➡️ 无后继 (最终产物)")
        
        tooltip_lines.append("")
        tooltip_lines.append(f"SMILES: {node[:70]}")
        
        tooltip = "\n".join(tooltip_lines)
        
        net.add_node(
            node,
            label=label,
            title=tooltip,
            color={'background': color, 'border': border_color},
            size=size,
            borderWidth=border_width,
            shape='dot',
            group=layer
        )
    
    # 添加边
    for edge in mol_graph.edges():
        source_out = out_deg.get(edge[0], 0)
        
        if source_out > 10:
            edge_color = '#FF6B6B'
            width = 3
        elif source_out > 5:
            edge_color = '#FFA726'
            width = 2
        else:
            edge_color = '#999999'
            width = 1.5
        
        net.add_edge(edge[0], edge[1], width=width, color={'color': edge_color})
    
    # 控制面板
    controls_html = f"""
    <div style="position: fixed; top: 10px; right: 10px; background: white; 
                padding: 18px; border: 3px solid #333; border-radius: 12px; 
                font-family: Arial; font-size: 13px; z-index: 1000; max-width: 360px;
                box-shadow: 0 8px 16px rgba(0,0,0,0.2);">
        
        <h2 style="margin: 0 0 12px 0; color: #333; border-bottom: 3px solid #333; padding-bottom: 8px;">
            🔋 SEI 网络分析
        </h2>
        
        <div style="background: #E3F2FD; padding: 12px; border-radius: 6px; margin-bottom: 12px;">
            <b style="color: #1565C0;">📊 SEI 层次 (颜色分组):</b><br>
            <div style="margin-top: 8px; font-size: 12px; line-height: 2.0;">
                <span style="color: #FF6B6B; font-size: 16px;">●</span> 
                <b>初始反应物</b>: {layer_stats.get('initial', 0)}<br>
                
                <span style="color: #FFD700; font-size: 16px;">●</span> 
                <b>中间体</b>: {layer_stats.get('intermediate', 0)}<br>
                
                <span style="color: #AAAAAA; font-size: 16px;">●</span> 
                <b>气体产物</b>: {layer_stats.get('gas', 0)}<br>
                
                <span style="color: #3366FF; font-size: 16px;">●</span> 
                <b>无机SEI</b>: {layer_stats.get('inorganic_sei', 0)}<br>
                
                <span style="color: #33CC33; font-size: 16px;">●</span> 
                <b>有机聚合物</b>: {layer_stats.get('organic_polymer', 0)}<br>
                
                <span style="color: #66FF66; font-size: 16px;">●</span> 
                <b>有机小分子</b>: {layer_stats.get('organic_small', 0)}<br>
            </div>
        </div>
        
        <div style="background: #FFF3E0; padding: 12px; border-radius: 6px; margin-bottom: 12px;">
            <b style="color: #E65100;">⚡ 编码说明:</b><br>
            <div style="font-size: 11px; margin-top: 6px; line-height: 1.5;">
            <b>节点大小</b> = 参与反应数<br>
            <b>边框</b>: 🔴粗红=高活性 | 🟢细绿=稳定<br>
            <b>箭头</b>: 🔴红粗=活性源 | ⚪灰=普通
            </div>
        </div>
        
        <div style="background: #E8F5E9; padding: 12px; border-radius: 6px;">
            <b style="color: #2E7D32;">🖱️ 悬停查看:</b><br>
            <div style="font-size: 11px; margin-top: 6px; line-height: 1.4;">
            • <b>⬅️ 前驱统计</b>: 有多少无机/有机<br>
            • <b>➡️ 后继统计</b>: 连到哪些类型<br>
            • 最终SEI组成分析
            </div>
        </div>
        
    </div>
    """
    
    # 保存
    print(f"Saving to {output_path}...")
    net.save_graph(output_path)
    
    with open(output_path, 'r', encoding='utf-8') as f:
        html = f.read()
    
    html = html.replace('</body>', controls_html + '</body>')
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)
    
    print()
    print("="*70)
    print("✓ Enhanced SEI visualization created!")
    print("="*70)
    print()
    print("Tooltip shows:")
    print("  ⬅️ Connected predecessors by SEI type")
    print("  ➡️ Connected successors by SEI type")
    print()

if __name__ == "__main__":
    create_enhanced_sei_view()
