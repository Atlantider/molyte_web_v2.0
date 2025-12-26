"""
Analyze what types of reactions produce EC
"""
import networkx as nx

G = nx.read_graphml('sei_full_network.graphml')

rxns = [n for n, d in G.nodes(data=True) if d.get('type') == 'reaction']
ec_prod = []

for r in rxns:
    succs = list(G.successors(r))
    if 'O=C1OCCO1' in succs:
        ec_prod.append(r)

print(f"Total reactions producing EC: {len(ec_prod)}")
print()
print("样例 (前10个):")
for i, r in enumerate(ec_prod[:10], 1):
    label = G.nodes[r].get('label', r)
    print(f"{i}. {label[:100]}")
