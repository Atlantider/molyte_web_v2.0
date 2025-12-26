"""
Improved SEI Network Visualization
====================================
Clear chemical names, reaction type labels, comprehensive legend
"""

from pyvis.network import Network
import networkx as nx
from collections import defaultdict, Counter
from rdkit import Chem

# Comprehensive SMILES to chemical name mapping
CHEM_NAMES = {
    # Initial reactants
    'O=C1OCCO1': 'EC',
    'COC(=O)OC': 'DMC', 
    '[Li+]': 'Li⁺',
    'F[P-](F)(F)(F)(F)F': 'PF₆⁻',
    
    # Li metal
    '[Li]': 'Li',
    '[Li][Li]': 'Li₂',
    
    # Key SEI products
    'F[Li]': 'LiF',
    '[Li]OC([Li])=O': 'Li₂CO₃',
    'C=C': 'C₂H₄',
    'O=C=O': 'CO₂',
    '[C-]#[O+]': 'CO',
    
    # Radicals
    '[O-]C(=O)OC': 'DMC-O•',
    '[CH2]OC(=O)OC': 'DMC-C•',
    'O=C1OCC[O]1': 'EC-O•',
    'O=C1[CH]CCO1': 'EC-C•',
    
    # Decomposition products
    'F[P](F)(F)(F)F': 'PF₅',
    '[F-]': 'F⁻',
    'FP(F)(F)F': 'PF₄',
    
    # Add more as needed
}

def get_chemical_name(smiles):
    """Convert SMILES to readable chemical name."""
    if smiles in CHEM_NAMES:
        return CHEM_NAMES[smiles]
    
    # Pattern-based naming
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles[:15] + '...' if len(smiles) > 15 else smiles
        
        # Count atoms
        atoms = Counter([atom.GetSymbol() for atom in mol.GetAtoms()])
        
        # Simple formulas
        if '.' in smiles:
            return 'Complex'
        elif atoms.get('Li', 0) > 0 and atoms.get('F', 0) > 0:
            return f'Li{atoms["Li"]}F{atoms["F"]}'
        elif atoms.get('Li', 0) > 0 and (atoms.get('O', 0) > 0 or atoms.get('C', 0) > 0):
            return 'Li-organic'
        elif '[' in smiles or ']' in smiles:
            return 'Radical' if '•' not in smiles else smiles[:10]
        else:
            # Use molecular formula
            formula = ''
            for elem in ['C', 'H', 'O', 'N', 'F', 'P', 'S', 'Li']:
                count = atoms.get(elem, 0)
                if count > 0:
                    formula += elem + (str(count) if count > 1 else '')
            return formula if formula else smiles[:10]
    except:
        return smiles[:15] + '...' if len(smiles) > 15 else smiles

def classify_species_detailed(smiles):
    """Classify species with detailed categories."""
    if smiles in ['O=C1OCCO1', 'COC(=O)OC', '[Li+]', 'F[P-](F)(F)(F)(F)F']:
        return 'reactant', '#FF4444', 50, 'Initial Reactant'
    if smiles in ['[Li]', '[Li][Li]']:
        return 'li_metal', '#000000', 45, 'Li Metal'
    if 'Li' in smiles and 'F' in smiles and '.' not in smiles:
        return 'lif', '#0066CC', 45, 'LiF Product'
    if smiles == 'O=C=O':
        return 'gas', '#00CC66', 40, 'Gas Product (CO₂)'
    if smiles in ['C=C', '[C-]#[O+]']:
        return 'gas', '#00CC66', 40, 'Gas Product'
    if ('C(=O)O' in smiles or 'CO3' in smiles) and 'Li' in smiles and '.' not in smiles:
        return 'carbonate', '#6600CC', 45, 'Li₂CO₃ Product'
    if '[O-]' in smiles or '[C]' in smiles or '[O]' in smiles:
        return 'radical', '#FF9900', 35, 'Radical'
    if '.' in smiles:
        return 'complex', '#9966FF', 25, 'Solvation Complex'
    return 'intermediate', '#99CC00', 30, 'Intermediate'

def get_reaction_type_from_nodes(source_smiles, target_smiles):
    """Infer reaction type from source and target."""
    source_type, _, _, _ = classify_species_detailed(source_smiles)
    target_type, _, _, _ = classify_species_detailed(target_smiles)
    
    # Pattern matching
    if source_type == 'reactant' and target_type == 'radical':
        if '[Li+]' in source_smiles:
            return 'e⁻ Reduction', '#0066FF'
        return 'Decomposition', '#CC0000'
    elif source_type == 'reactant' and '[Li]' in target_smiles:
        return 'e⁻ Reduction', '#0066FF'
    elif source_type == 'radical' and target_type == 'radical':
        return 'Radical Coupling', '#00AA00'
    elif source_type == 'radical' and '.' not in source_smiles and '.' in target_smiles:
        return 'Coordination', '#888888'
    elif 'PF6' in source_smiles or 'F[P-]' in source_smiles:
        return 'Dissociation', '#CC00CC'
    elif source_type == 'li_metal' and target_type == 'complex':
        return 'Coordination', '#888888'
    elif target_type == 'gas':
        return 'Decomposition', '#CC0000'
    elif target_type in ['lif', 'carbonate']:
        return 'SEI Formation', '#0000CC'
    else:
        return 'Reaction', '#666666'

def create_improved_visualization(
    graphml_path='sei_full_network.graphml',
    output_path='sei_network_clear.html',
    max_nodes=80
):
    """Create improved visualization with clear labels."""
    
    print("="*60)
    print("IMPROVED SEI NETWORK VISUALIZATION")
    print("="*60)
    print()
    
    # Load network
    print("Loading network...")
    G = nx.read_graphml(graphml_path)
    
    molecule_nodes = [n for n, d in G.nodes(data=True) if d.get('type') != 'reaction']
    reaction_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'reaction']
    
    # Create molecule network with reaction info
    print("Creating molecule network...")
    mol_graph = nx.DiGraph()
    
    for mol in molecule_nodes:
        mol_graph.add_node(mol, **G.nodes[mol])
    
    # Store reaction types on edges
    edge_reactions = defaultdict(list)
    for rxn in reaction_nodes:
        preds = [n for n in G.predecessors(rxn) if n in molecule_nodes]
        succs = [n for n in G.successors(rxn) if n in molecule_nodes]
        rxn_name = G.nodes[rxn].get('label', 'reaction')
        
        for p in preds:
            for s in succs:
                if p != s:
                    mol_graph.add_edge(p, s)
                    edge_reactions[(p, s)].append(rxn_name)
    
    # Filter by priority
    print("Filtering to key species...")
    node_scores = {}
    for node in mol_graph.nodes():
        _, _, _, category = classify_species_detailed(node)
        degree = mol_graph.degree(node)
        gen = int(mol_graph.nodes[node].get('generation', 0))
        
        # Priority: reactants > products > radicals > intermediates > complexes
        priority_map = {
            'Initial Reactant': 5,
            'Li Metal': 4,
            'LiF Product': 4,
            'Li₂CO₃ Product': 4,
            'Gas Product': 3,
            'Radical': 3,
            'Intermediate': 2,
            'Solvation Complex': 1
        }
        priority = priority_map.get(category, 1)
        
        score = priority * degree - gen * 0.3
        node_scores[node] = score
    
    important_nodes = sorted(node_scores.keys(), key=lambda n: node_scores[n], reverse=True)[:max_nodes]
    mol_graph = mol_graph.subgraph(important_nodes).copy()
    
    print(f"  Selected {len(mol_graph.nodes())} key species")
    print()
    
    # Create Pyvis network
    print("Creating interactive visualization...")
    net = Network(
        height='950px',
        width='100%',
        bgcolor='#FAFAFA',
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
          "levelSeparation": 200,
          "nodeSpacing": 180,
          "treeSpacing": 220,
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
        "keyboard": true,
        "tooltipDelay": 100
      },
      "nodes": {
        "font": {"size": 16, "face": "Arial", "bold": true},
        "borderWidth": 3,
        "borderWidthSelected": 5
      },
      "edges": {
        "smooth": {"type": "cubicBezier"},
        "arrows": {"to": {"enabled": true, "scaleFactor": 1.0}},
        "width": 2.5,
        "font": {"size": 11, "align": "top", "background": "white"}
      }
    }
    """)
    
    # Add nodes
    print("Adding nodes...")
    species_stats = defaultdict(int)
    
    for node in mol_graph.nodes():
        species_type, color, size, category = classify_species_detailed(node)
        species_stats[category] += 1
        
        chem_name = get_chemical_name(node)
        gen = int(mol_graph.nodes[node].get('generation', 0))
        degree = mol_graph.degree(node)
        
        # Tooltip with full info
        tooltip = (
            f"<b>{chem_name}</b><br>"
            f"Category: {category}<br>"
            f"Generation: {gen}<br>"
            f"Connections: {degree}<br>"
            f"<br>"
            f"SMILES: {node[:70]}{'...' if len(node) > 70 else ''}"
        )
        
        net.add_node(
            node,
            label=chem_name,
            title=tooltip,
            color=color,
            size=size,
            level=gen,
            borderWidth=3,
            shape='dot'
        )
    
    # Add edges with reaction type labels
    print("Adding edges with reaction types...")
    edge_stats = defaultdict(int)
    
    for edge in mol_graph.edges():
        rxn_type, rxn_color = get_reaction_type_from_nodes(edge[0], edge[1])
        edge_stats[rxn_type] += 1
        
        net.add_edge(
            edge[0], edge[1],
            label=rxn_type,
            color={'color': rxn_color, 'highlight': '#FF0000'},
            width=2.5,
            font={'size': 10, 'align': 'top'}
        )
    
    # Create comprehensive legend
    legend_html = f"""
    <div style="position: fixed; top: 10px; left: 10px; background: white; 
                padding: 20px; border: 3px solid #333; border-radius: 12px; 
                font-family: Arial; font-size: 13px; z-index: 1000; max-width: 380px;
                box-shadow: 0 6px 12px rgba(0,0,0,0.15); max-height: 90vh; overflow-y: auto;">
        
        <h2 style="margin: 0 0 15px 0; color: #333; border-bottom: 3px solid #333; padding-bottom: 10px; font-size: 18px;">
            🔬 SEI Network Viewer
        </h2>
        
        <!-- Environment Info -->
        <div style="margin-bottom: 15px; background: #E8F5E9; padding: 12px; border-radius: 6px; border-left: 4px solid #4CAF50;">
            <b style="color: #2E7D32;">⚡ Conditions:</b><br>
            <span style="font-size: 12px;">
            • Electrode: <b>Anode</b><br>
            • Voltage: <b>0.1V vs Li/Li⁺</b><br>
            • Temperature: <b>298K</b><br>
            • Species shown: <b>{len(mol_graph.nodes())}</b> (most important)
            </span>
        </div>
        
        <!-- Node Legend -->
        <div style="margin-bottom: 15px; background: #FFF9C4; padding: 12px; border-radius: 6px; border-left: 4px solid #FFC107;">
            <b style="color: #F57C00;">🔵 Node (Species) Legend:</b><br>
            <div style="margin-top: 8px; font-size: 12px; line-height: 1.8;">
    """
    
    # Add species categories
    category_example = {
        'Initial Reactant': ('#FF4444', 'EC, DMC, Li⁺, PF₆⁻'),
        'Li Metal': ('#000000', 'Li, Li₂'),
        'LiF Product': ('#0066CC', 'LiF (key SEI)'),
        'Li₂CO₃ Product': ('#6600CC', 'Li₂CO₃ (key SEI)'),
        'Gas Product': ('#00CC66', 'CO₂, CO, C₂H₄'),
        'Radical': ('#FF9900', 'EC-O•, DMC-C•'),
        'Intermediate': ('#99CC00', 'Reaction intermediates'),
        'Solvation Complex': ('#9966FF', 'Li⁺(EC)₂'),
    }
    
    for cat, (color, example) in category_example.items():
        count = species_stats.get(cat, 0)
        if count > 0 or cat in ['Initial Reactant', 'Li Metal', 'Radical']:
            legend_html += f"""
                <span style="display: inline-block; width: 14px; height: 14px; 
                             background: {color}; border: 2px solid black; 
                             border-radius: 50%; margin-right: 6px;"></span>
                <b>{cat}</b> ({count})<br>
                <span style="margin-left: 24px; color: #666; font-size: 10px;">{example}</span><br>
            """
    
    legend_html += """
            </div>
        </div>
        
        <!-- Arrow Legend -->
        <div style="margin-bottom: 15px; background: #E3F2FD; padding: 12px; border-radius: 6px; border-left: 4px solid #2196F3;">
            <b style="color: #1565C0;">➡️ Arrow (Reaction) Legend:</b><br>
            <div style="margin-top: 8px; font-size: 12px; line-height: 1.8;">
                <span style="color: #0066FF; font-weight: bold;">━━➤</span> e⁻ Reduction<br>
                <span style="font-size: 10px; margin-left: 20px; color: #666;">Li⁺ → Li, EC → EC-O•</span><br>
                
                <span style="color: #00AA00; font-weight: bold;">━━➤</span> Radical Coupling<br>
                <span style="font-size: 10px; margin-left: 20px; color: #666;">R• + R• → R-R</span><br>
                
                <span style="color: #CC0000; font-weight: bold;">━━➤</span> Decomposition<br>
                <span style="font-size: 10px; margin-left: 20px; color: #666;">EC → CO₂ + C₂H₄</span><br>
                
                <span style="color: #CC00CC; font-weight: bold;">━━➤</span> Dissociation<br>
                <span style="font-size: 10px; margin-left: 20px; color: #666;">PF₆⁻ → PF₅ + F⁻</span><br>
                
                <span style="color: #888888; font-weight: bold;">━━➤</span> Coordination<br>
                <span style="font-size: 10px; margin-left: 20px; color: #666;">Li + EC → Li(EC)</span><br>
                
                <span style="color: #0000CC; font-weight: bold;">━━➤</span> SEI Formation<br>
                <span style="font-size: 10px; margin-left: 20px; color: #666;">→ LiF, Li₂CO₃</span><br>
            </div>
        </div>
        
        <!-- Node Size -->
        <div style="margin-bottom: 15px; background: #FCE4EC; padding: 12px; border-radius: 6px; border-left: 4px solid #E91E63;">
            <b style="color: #C2185B;">📏 Node Size:</b><br>
            <div style="font-size: 12px; margin-top: 6px;">
                Large ● = High importance (reactants, key products)<br>
                Medium ● = Intermediates, radicals<br>
                Small ● = Complexes
            </div>
        </div>
        
        <!-- Instructions -->
        <div style="background: #F5F5F5; padding: 12px; border-radius: 6px; border: 1px solid #DDD; font-size: 12px;">
            <b style="color: #424242;">🖱️ How to Use:</b><br>
            <div style="margin-top: 6px; line-height: 1.6;">
            • <b>Scroll</b> = Zoom in/out<br>
            • <b>Drag</b> = Pan the view<br>
            • <b>Click node</b> = View details<br>
            • <b>Hover</b> = Quick info
            </div>
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
    print("="*60)
    print("✓ Improved visualization created!")
    print(f"  File: {output_path}")
    print("="*60)
    print()
    print("Key Improvements:")
    print("  ✓ Chemical names instead of SMILES (EC, DMC, Li⁺)")
    print("  ✓ Reaction types labeled on arrows")
    print("  ✓ Comprehensive legend with all symbols explained")
    print("  ✓ Clear node sizes (importance-based)")
    print("  ✓ Color-coded by chemical function")
    print()
    print("Legend includes:")
    print("  • Node colors (8 categories)")
    print("  • Arrow colors (6 reaction types)")
    print("  • Node size meaning")
    print("  • Interactive instructions")
    print()

if __name__ == "__main__":
    create_improved_visualization()
