"""
Specialized View: Ring Opening and Polymerization Products
============================================================
Focus on showing EC ring opening and oligomerization pathways
"""

from pyvis.network import Network
import networkx as nx
from collections import defaultdict, Counter
from rdkit import Chem

def is_ring_opened(smiles):
    """Check if species is a ring-opened product."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        
        ring_info = mol.GetRingInfo()
        # Ring-opened EC: contains OCCO but no ring
        if 'OCCO' in smiles and ring_info.NumRings() == 0 and '.' not in smiles:
            return True
        return False
    except:
        return False

def is_oligomer(smiles):
    """Check if species is an oligomer (>15 heavy atoms)."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or '.' in smiles:
            return False
        return mol.GetNumHeavyAtoms() > 15
    except:
        return False

def get_chemical_name_simple(smiles):
    """Get simple chemical name."""
    names = {
        'O=C1OCCO1': 'EC',
        'COC(=O)OC': 'DMC',
        'OCCO': 'Ethylene Glycol',
        'CO[C](OC)OCCO': 'Ring-Opened DMC-EC',
    }
    
    if smiles in names:
        return names[smiles]
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            atoms = mol.GetNumHeavyAtoms()
            if is_oligomer(smiles):
                return f'Oligomer ({atoms}A)'
            elif is_ring_opened(smiles):
                return f'Ring-Opened ({atoms}A)'
        return smiles[:15] + '...'
    except:
        return smiles[:15] + '...'

def create_ring_opening_view(
    graphml_path='sei_full_network.graphml',
    output_path='sei_network_ring_opening.html'
):
    """Create specialized view for ring opening and polymerization."""
    
    print("="*60)
    print("RING OPENING & POLYMERIZATION VIEW")
    print("="*60)
    print()
    
    # Load network
    print("Loading network...")
    G = nx.read_graphml(graphml_path)
    
    molecule_nodes = [n for n, d in G.nodes(data=True) if d.get('type') != 'reaction']
    reaction_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'reaction']
    
    # Filter to relevant species
    print("Filtering to ring-opening and polymerization species...")
    
    relevant_species = set()
    
    # Add initial reactants
    initial = ['O=C1OCCO1', 'COC(=O)OC', '[Li+]', 'F[P-](F)(F)(F)(F)F']
    for s in initial:
        if s in molecule_nodes:
            relevant_species.add(s)
    
    # Add ring-opened and oligomer species
    for node in molecule_nodes:
        if is_ring_opened(node) or is_oligomer(node):
            relevant_species.add(node)
    
    print(f"  Found {len(relevant_species)} relevant species")
    
    # Build subgraph with connections
    mol_graph = nx.DiGraph()
    for mol in relevant_species:
        if mol in molecule_nodes:
            mol_graph.add_node(mol, **G.nodes[mol])
    
    # Add edges through reactions
    for rxn in reaction_nodes:
        preds = [n for n in G.predecessors(rxn) if n in relevant_species]
        succs = [n for n in G.successors(rxn) if n in relevant_species]
        
        for p in preds:
            for s in succs:
                if p != s:
                    mol_graph.add_edge(p, s)
    
    # Add intermediate connections (pathway nodes)
    # Find species that connect ring-intact to ring-opened
    for mol in molecule_nodes:
        if mol in relevant_species:
            continue
        
        # Check if this connects relevant species
        has_relevant_pred = any(p in relevant_species for p in G.predecessors(mol) if p in molecule_nodes)
        has_relevant_succ = any(s in relevant_species for s in G.successors(mol) if s in molecule_nodes)
        
        if has_relevant_pred and has_relevant_succ:
            # This is a bridge species
            relevant_species.add(mol)
    
    # Rebuild edges with bridge species
    mol_graph = nx.DiGraph()
    for mol in relevant_species:
        if mol in molecule_nodes:
            mol_graph.add_node(mol, **G.nodes[mol])
    
    for rxn in reaction_nodes:
        preds = [n for n in G.predecessors(rxn) if n in relevant_species]
        succs = [n for n in G.successors(rxn) if n in relevant_species]
        
        for p in preds:
            for s in succs:
                if p != s:
                    mol_graph.add_edge(p, s)
    
    print(f"  Final network: {len(mol_graph.nodes())} species, {len(mol_graph.edges())} connections")
    print()
    
    # Create Pyvis network
    print("Creating visualization...")
    net = Network(
        height='950px',
        width='100%',
        bgcolor='#FAFAFA',
        font_color='black',
        directed=True,
        notebook=False
    )
    
    net.set_options("""
    {
      "layout": {
        "hierarchical": {
          "enabled": true,
          "levelSeparation": 220,
          "nodeSpacing": 180,
          "direction": "UD"
        }
      },
      "physics": {"enabled": false},
      "interaction": {"hover": true, "navigationButtons": true},
      "nodes": {"font": {"size": 14, "face": "Arial", "bold": true}},
      "edges": {
        "smooth": {"type": "cubicBezier"},
        "arrows": {"to": {"enabled": true, "scaleFactor": 0.9}},
        "width": 2.5
      }
    }
    """)
    
    # Add nodes
    ring_opened_count = 0
    oligomer_count = 0
    
    for node in mol_graph.nodes():
        gen = int(mol_graph.nodes[node].get('generation', 0))
        
        # Classify
        if node in initial:
            color = '#FF4444'
            size = 50
            category = 'Initial Reactant'
        elif is_oligomer(node):
            color = '#9933FF'
            size = 45
            category = 'Oligomer'
            oligomer_count += 1
        elif is_ring_opened(node):
            color = '#33CCFF'
            size = 40
            category = 'Ring-Opened'
            ring_opened_count += 1
        else:
            color = '#FFAA00'
            size = 30
            category = 'Intermediate'
        
        label = get_chemical_name_simple(node)
        
        tooltip = (
            f"<b>{label}</b><br>"
            f"Category: {category}<br>"
            f"Generation: {gen}<br>"
            f"SMILES: {node[:80]}{'...' if len(node) > 80 else ''}"
        )
        
        net.add_node(
            node,
            label=label,
            title=tooltip,
            color=color,
            size=size,
            level=gen
        )
    
    # Add edges
    for edge in mol_graph.edges():
        net.add_edge(
            edge[0], edge[1],
            color={'color': '#666666', 'highlight': '#FF0000'},
            width=2.5
        )
    
    # Legend
    legend_html = f"""
    <div style="position: fixed; top: 10px; left: 10px; background: white; 
                padding: 20px; border: 3px solid #333; border-radius: 12px; 
                font-family: Arial; font-size: 14px; z-index: 1000; max-width: 400px;
                box-shadow: 0 6px 12px rgba(0,0,0,0.15);">
        
        <h2 style="margin: 0 0 15px 0; color: #333; border-bottom: 3px solid #333; padding-bottom: 10px;">
            🔓 Ring Opening & Polymerization
        </h2>
        
        <div style="background: #FFF3E0; padding: 12px; border-radius: 6px; margin-bottom: 15px; border-left: 4px solid #FF9800;">
            <b style="color: #E65100;">📊 Statistics:</b><br>
            <div style="margin-top: 8px; font-size: 13px;">
            • Ring-opened: <b>{ring_opened_count}</b> species<br>
            • Oligomers: <b>{oligomer_count}</b> species<br>
            • Total shown: <b>{len(mol_graph.nodes())}</b>
            </div>
        </div>
        
        <div style="background: #E8F5E9; padding: 12px; border-radius: 6px; margin-bottom: 15px;">
            <b style="color: #2E7D32;">🔵 Node Categories:</b><br>
            <div style="margin-top: 8px; font-size: 13px; line-height: 2.0;">
                <span style="display: inline-block; width: 16px; height: 16px; 
                             background: #FF4444; border: 2px solid black; 
                             border-radius: 50%; margin-right: 8px;"></span>
                <b>Initial Reactant</b> (EC, DMC)<br>
                
                <span style="display: inline-block; width: 16px; height: 16px; 
                             background: #33CCFF; border: 2px solid black; 
                             border-radius: 50%; margin-right: 8px;"></span>
                <b>Ring-Opened</b> (OCCO derivatives)<br>
                
                <span style="display: inline-block; width: 16px; height: 16px; 
                             background: #9933FF; border: 2px solid black; 
                             border-radius: 50%; margin-right: 8px;"></span>
                <b>Oligomer</b> (\u003e15 atoms)<br>
                
                <span style="display: inline-block; width: 16px; height: 16px; 
                             background: #FFAA00; border: 2px solid black; 
                             border-radius: 50%; margin-right: 8px;"></span>
                <b>Intermediate</b>
            </div>
        </div>
        
        <div style="background: #FCE4EC; padding: 12px; border-radius: 6px; margin-bottom: 15px;">
            <b style="color: #C2185B;">⚙️ Key Pathways:</b><br>
            <div style="font-size: 12px; margin-top: 8px; line-height: 1.6;">
            1. <b>EC Ring Opening:</b><br>
               <span style="margin-left: 15px;">EC → EC-O• → Ring-Opened</span><br>
            2. <b>Oligomerization:</b><br>
               <span style="margin-left: 15px;">Fragment + Fragment → Oligomer</span>
            </div>
        </div>
        
        <div style="background: #FFF9C4; padding: 12px; border-radius: 6px; border-left: 4px solid #FFC107;">
            <b style="color: #F57C00;">⚠️ Note:</b><br>
            <div style="font-size: 12px; margin-top: 6px;">
            开环和聚合反应 <b>存在于网络中</b>,<br>
            但通过 <b>Bond Breaking</b> operator<br>
            产生,而非单独的 RingOpening/<br>
            Polymerization operators。
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
    print("✓ Ring opening & polymerization view created!")
    print(f"  File: {output_path}")
    print("="*60)
    print()
    print(f"Showing {ring_opened_count} ring-opened products")
    print(f"Showing {oligomer_count} oligomers")
    print()

if __name__ == "__main__":
    create_ring_opening_view()
