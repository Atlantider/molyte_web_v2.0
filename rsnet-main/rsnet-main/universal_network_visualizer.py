#!/usr/bin/env python3
"""
通用反应网络可视化器
自动识别反应模板类型，适用于任何化学体系
"""

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import networkx as nx
import numpy as np
from pathlib import Path
from collections import defaultdict, Counter
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

sys.path.insert(0, '.')

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.core.reaction import Reaction

# 导入生成器
from iterative_network_demo import IterativeNetworkGenerator


class UniversalNetworkVisualizer:
    """通用反应网络可视化器"""
    
    def __init__(self, generator, output_dir):
        self.generator = generator
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # 自动分析反应模板类型
        self.template_types = self._analyze_reaction_templates()
        self.template_colors = self._generate_template_colors()
        
        # 自动识别重要物种
        self.key_species = self._identify_key_species()
        
        # 计算重要性分数
        self._calculate_importance_scores()
        self._identify_dominant_pathways()
    
    def _analyze_reaction_templates(self):
        """自动分析反应模板类型"""
        
        template_analysis = {}
        
        for rxn in self.generator.all_reactions:
            template_name = rxn.name
            
            if template_name not in template_analysis:
                template_analysis[template_name] = {
                    'count': 0,
                    'reactions': [],
                    'avg_reactants': 0,
                    'avg_products': 0,
                    'avg_energy': 0,
                    'avg_barrier': 0,
                    'involves_metal': False,
                    'involves_organic': False,
                    'bond_changes': []
                }
            
            analysis = template_analysis[template_name]
            analysis['count'] += 1
            analysis['reactions'].append(rxn)
            analysis['avg_reactants'] += len(rxn.reactants)
            analysis['avg_products'] += len(rxn.products)
            analysis['avg_energy'] += rxn.reaction_energy
            analysis['avg_barrier'] += rxn.activation_energy
            
            # 分析涉及的元素类型
            all_molecules = rxn.reactants + rxn.products
            for mol in all_molecules:
                if self._contains_metal(mol):
                    analysis['involves_metal'] = True
                if self._contains_organic_carbon(mol):
                    analysis['involves_organic'] = True
        
        # 计算平均值
        for template_name, analysis in template_analysis.items():
            count = analysis['count']
            analysis['avg_reactants'] /= count
            analysis['avg_products'] /= count
            analysis['avg_energy'] /= count
            analysis['avg_barrier'] /= count
        
        return template_analysis
    
    def _contains_metal(self, mol):
        """检查分子是否包含金属"""
        try:
            if hasattr(mol, 'rdkit_mol') and mol.rdkit_mol:
                for atom in mol.rdkit_mol.GetAtoms():
                    atomic_num = atom.GetAtomicNum()
                    # 常见金属原子序数
                    metals = [3, 11, 19, 37, 55, 87,  # 碱金属
                             4, 12, 20, 38, 56, 88,   # 碱土金属
                             21, 22, 23, 24, 25, 26, 27, 28, 29, 30,  # 过渡金属
                             13, 31, 49, 50, 81, 82, 83]  # 其他金属
                    if atomic_num in metals:
                        return True
            # 简单SMILES检查
            metal_symbols = ['Li', 'Na', 'K', 'Mg', 'Ca', 'Al', 'Zn', 'Fe', 'Cu', 'Ni']
            for symbol in metal_symbols:
                if symbol in mol.smiles:
                    return True
        except:
            pass
        return False
    
    def _contains_organic_carbon(self, mol):
        """检查是否包含有机碳"""
        try:
            if hasattr(mol, 'rdkit_mol') and mol.rdkit_mol:
                carbon_count = sum(1 for atom in mol.rdkit_mol.GetAtoms() 
                                 if atom.GetAtomicNum() == 6)
                return carbon_count > 0
            # 简单检查
            return 'C' in mol.smiles and mol.smiles not in ['CO', 'CO2', '[CO]', 'O=C=O']
        except:
            return 'C' in mol.smiles
    
    def _generate_template_colors(self):
        """为反应模板生成颜色"""
        
        templates = list(self.template_types.keys())
        n_templates = len(templates)
        
        # 使用色彩丰富的调色板
        if n_templates <= 10:
            colors = sns.color_palette("Set3", n_templates)
        else:
            colors = sns.color_palette("husl", n_templates)
        
        template_colors = {}
        for i, template in enumerate(templates):
            template_colors[template] = colors[i]
        
        return template_colors
    
    def _identify_key_species(self):
        """自动识别关键物种"""
        
        key_species = {
            'initial': [],      # 初始物种
            'final': [],        # 最终产物
            'intermediate': [], # 重要中间体
            'metal_containing': [],  # 含金属物种
            'organic': []       # 有机物种
        }
        
        # 按代数分类
        by_generation = defaultdict(list)
        for smiles, (mol, gen) in self.generator.all_molecules.items():
            by_generation[gen].append(mol)
        
        # 初始物种
        key_species['initial'] = by_generation[0]
        
        # 最终产物（最高代数）
        max_gen = max(by_generation.keys())
        key_species['final'] = by_generation[max_gen]
        
        # 按化学性质分类
        for smiles, (mol, gen) in self.generator.all_molecules.items():
            if self._contains_metal(mol):
                key_species['metal_containing'].append(mol)
            if self._contains_organic_carbon(mol):
                key_species['organic'].append(mol)
            
            # 中间体：参与多个反应的物种
            reaction_count = 0
            for rxn in self.generator.all_reactions:
                if mol in rxn.reactants or mol in rxn.products:
                    reaction_count += 1
            
            if reaction_count >= 3 and gen > 0:  # 参与3个以上反应的非初始物种
                key_species['intermediate'].append(mol)
        
        return key_species
    
    def _calculate_importance_scores(self):
        """计算重要性分数（通用算法）"""
        
        # 物种重要性 = 度数 + 代数权重 + 化学性质权重
        self.species_importance = {}
        max_gen = max(gen for smiles, (mol, gen) in self.generator.all_molecules.items())
        
        for smiles, (mol, gen) in self.generator.all_molecules.items():
            # 计算度数
            degree = sum(1 for rxn in self.generator.all_reactions 
                        if mol in rxn.reactants or mol in rxn.products)
            
            # 代数权重（初始和最终产物更重要）
            gen_weight = 3 if gen == 0 else (2 if gen == max_gen else 1)
            
            # 化学性质权重
            chem_weight = 0
            if self._contains_metal(mol):
                chem_weight += 2  # 金属物种通常重要
            if mol in self.key_species['intermediate']:
                chem_weight += 2  # 中间体重要
            
            importance = degree * 2 + gen_weight + chem_weight
            self.species_importance[mol.name] = importance
        
        # 反应重要性 = 能量权重 + 模板频次权重 + 产物重要性
        self.reaction_importance = {}
        template_counts = Counter(rxn.name for rxn in self.generator.all_reactions)
        
        for i, rxn in enumerate(self.generator.all_reactions):
            # 能量权重（放热反应更重要）
            energy_weight = max(0, -rxn.reaction_energy / 20.0)
            
            # 模板频次权重（常见模板更重要）
            template_weight = template_counts[rxn.name] / len(self.generator.all_reactions)
            
            # 产物重要性
            product_weight = sum(self.species_importance.get(p.name, 0) 
                               for p in rxn.products) / len(rxn.products)
            
            importance = energy_weight + template_weight * 5 + product_weight / 10
            self.reaction_importance[f"R{i}"] = importance
    
    def _identify_dominant_pathways(self):
        """识别主导反应路径（通用算法）"""
        
        # 自动选择起点和终点
        initial_species = [mol.name for mol in self.key_species['initial']]
        
        # 终点：最重要的最终产物
        final_candidates = self.key_species['final']
        if not final_candidates:
            # 如果没有明确的最终产物，选择重要性最高的物种
            final_candidates = sorted(self.generator.all_molecules.values(),
                                    key=lambda x: self.species_importance.get(x[0].name, 0),
                                    reverse=True)[:3]
            final_candidates = [mol for mol, gen in final_candidates]
        
        target_species = [mol.name for mol in final_candidates[:3]]  # 取前3个
        
        self.dominant_pathways = []
        
        # 构建反应图
        G = nx.DiGraph()
        
        # 添加节点
        for smiles, (mol, gen) in self.generator.all_molecules.items():
            G.add_node(mol.name, type='species', generation=gen)
        
        for i, rxn in enumerate(self.generator.all_reactions):
            rxn_id = f"R{i}"
            G.add_node(rxn_id, type='reaction', 
                      template=rxn.name,
                      energy=rxn.reaction_energy,
                      barrier=rxn.activation_energy)
            
            # 添加边
            for reactant in rxn.reactants:
                G.add_edge(reactant.name, rxn_id)
            for product in rxn.products:
                G.add_edge(rxn_id, product.name)
        
        # 寻找路径
        for initial in initial_species:
            for target in target_species:
                if initial in G and target in G and initial != target:
                    try:
                        paths = list(nx.all_simple_paths(G, initial, target, cutoff=8))
                        for path in paths[:2]:  # 每对起终点取前2条路径
                            # 计算路径分数
                            path_energy = 0
                            path_reactions = []
                            for i in range(len(path)-1):
                                if G.nodes[path[i+1]].get('type') == 'reaction':
                                    path_energy += G.nodes[path[i+1]].get('energy', 0)
                                    path_reactions.append(path[i+1])
                            
                            self.dominant_pathways.append({
                                'path': path,
                                'reactions': path_reactions,
                                'total_energy': path_energy,
                                'score': -path_energy - len(path_reactions) * 5,
                                'start': initial,
                                'end': target
                            })
                    except nx.NetworkXNoPath:
                        continue
        
        # 按分数排序，保留前5条
        self.dominant_pathways.sort(key=lambda x: x['score'], reverse=True)
        self.dominant_pathways = self.dominant_pathways[:5]
    
    def create_global_network_view(self):
        """创建全局网络视图"""
        
        fig, ax = plt.subplots(figsize=(20, 14))
        
        # 创建双部图
        G = nx.Graph()
        
        # 添加物种节点
        species_nodes = []
        species_colors = []
        species_sizes = []
        
        for smiles, (mol, gen) in self.generator.all_molecules.items():
            G.add_node(mol.name, bipartite=0, type='species')
            species_nodes.append(mol.name)
            
            # 颜色基于化学性质
            if mol in self.key_species['initial']:
                color = '#FF6B6B'  # 红色 - 初始
            elif mol in self.key_species['final']:
                color = '#4ECDC4'  # 青色 - 最终
            elif mol in self.key_species['metal_containing']:
                color = '#FFD93D'  # 黄色 - 含金属
            elif mol in self.key_species['organic']:
                color = '#95E1D3'  # 浅青 - 有机
            else:
                color = '#CCCCCC'  # 灰色 - 其他
            
            species_colors.append(color)
            
            # 大小基于重要性
            importance = self.species_importance.get(mol.name, 1)
            species_sizes.append(importance * 100 + 200)
        
        # 添加反应节点
        reaction_nodes = []
        reaction_colors = []
        reaction_sizes = []
        
        for i, rxn in enumerate(self.generator.all_reactions):
            rxn_id = f"R{i}"
            G.add_node(rxn_id, bipartite=1, type='reaction')
            reaction_nodes.append(rxn_id)
            
            # 颜色基于模板类型
            template = rxn.name
            reaction_colors.append(self.template_colors.get(template, '#CCCCCC'))
            
            # 大小基于重要性
            importance = self.reaction_importance.get(rxn_id, 1)
            reaction_sizes.append(importance * 80 + 150)
        
        # 添加边
        for rxn in self.generator.all_reactions:
            rxn_id = f"R{self.generator.all_reactions.index(rxn)}"
            for reactant in rxn.reactants:
                G.add_edge(reactant.name, rxn_id)
            for product in rxn.products:
                G.add_edge(rxn_id, product.name)
        
        # 布局
        pos = nx.spring_layout(G, k=3, iterations=50, seed=42)
        
        # 绘制物种节点
        nx.draw_networkx_nodes(G, pos, nodelist=species_nodes,
                              node_color=species_colors,
                              node_size=species_sizes,
                              node_shape='o',
                              alpha=0.8,
                              edgecolors='black',
                              linewidths=2)
        
        # 绘制反应节点
        nx.draw_networkx_nodes(G, pos, nodelist=reaction_nodes,
                              node_color=reaction_colors,
                              node_size=reaction_sizes,
                              node_shape='s',
                              alpha=0.9,
                              edgecolors='black',
                              linewidths=2)
        
        # 绘制边
        nx.draw_networkx_edges(G, pos,
                              edge_color='gray',
                              width=1.5,
                              alpha=0.6)
        
        # 添加标签
        labels = {}
        for node in species_nodes:
            labels[node] = node
        for i, node in enumerate(reaction_nodes):
            rxn = self.generator.all_reactions[i]
            labels[node] = f"{rxn.name}\nΔG={rxn.reaction_energy:.1f}"
        
        nx.draw_networkx_labels(G, pos, labels, font_size=8, font_weight='bold')
        
        # 图例：物种类型
        species_legend = [
            patches.Patch(color='#FF6B6B', label='初始物种'),
            patches.Patch(color='#4ECDC4', label='最终产物'),
            patches.Patch(color='#FFD93D', label='含金属物种'),
            patches.Patch(color='#95E1D3', label='有机物种'),
            patches.Patch(color='#CCCCCC', label='其他物种')
        ]
        
        # 图例：反应模板
        template_legend = []
        for template, color in self.template_colors.items():
            count = self.template_types[template]['count']
            template_legend.append(patches.Patch(color=color, 
                                               label=f'{template} ({count})'))
        
        # 分两列显示图例
        legend1 = ax.legend(handles=species_legend, loc='upper left', 
                           title='物种类型', title_fontsize=12)
        ax.add_artist(legend1)
        
        ax.legend(handles=template_legend, loc='upper right', 
                 title='反应模板类型', title_fontsize=12)
        
        ax.set_title('通用反应网络全局视图\n(物种=圆形, 反应=方形)', 
                    fontsize=16, fontweight='bold', pad=20)
        ax.axis('off')
        
        # 保存
        output_file = self.output_dir / 'universal_global_network.png'
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✅ 通用全局网络图已保存: {output_file}")
        plt.close()

    def create_template_analysis_dashboard(self):
        """创建反应模板分析仪表板"""

        fig = plt.figure(figsize=(20, 16))

        # 创建网格布局
        gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

        # 1. 模板类型分布
        ax1 = fig.add_subplot(gs[0, 0])
        templates = list(self.template_types.keys())
        counts = [self.template_types[t]['count'] for t in templates]
        colors = [self.template_colors[t] for t in templates]

        wedges, texts, autotexts = ax1.pie(counts, labels=templates, colors=colors,
                                          autopct='%1.1f%%', startangle=90)
        ax1.set_title('反应模板分布', fontweight='bold')

        # 2. 模板能量特征
        ax2 = fig.add_subplot(gs[0, 1])
        energies = [self.template_types[t]['avg_energy'] for t in templates]
        bars = ax2.bar(range(len(templates)), energies, color=colors, alpha=0.8)
        ax2.set_xticks(range(len(templates)))
        ax2.set_xticklabels(templates, rotation=45, ha='right')
        ax2.set_ylabel('平均反应能 (kcal/mol)')
        ax2.set_title('模板能量特征', fontweight='bold')
        ax2.axhline(y=0, color='red', linestyle='--', alpha=0.5)
        ax2.grid(True, alpha=0.3)

        # 3. 模板活化能特征
        ax3 = fig.add_subplot(gs[0, 2])
        barriers = [self.template_types[t]['avg_barrier'] for t in templates]
        bars = ax3.bar(range(len(templates)), barriers, color=colors, alpha=0.8)
        ax3.set_xticks(range(len(templates)))
        ax3.set_xticklabels(templates, rotation=45, ha='right')
        ax3.set_ylabel('平均活化能 (kcal/mol)')
        ax3.set_title('模板活化能特征', fontweight='bold')
        ax3.grid(True, alpha=0.3)

        # 4. 模板化学性质分析
        ax4 = fig.add_subplot(gs[1, :2])

        # 统计各模板的化学性质
        metal_involvement = []
        organic_involvement = []

        for template in templates:
            analysis = self.template_types[template]
            metal_involvement.append(1 if analysis['involves_metal'] else 0)
            organic_involvement.append(1 if analysis['involves_organic'] else 0)

        x = np.arange(len(templates))
        width = 0.35

        bars1 = ax4.bar(x - width/2, metal_involvement, width,
                       label='涉及金属', color='#FFD93D', alpha=0.8)
        bars2 = ax4.bar(x + width/2, organic_involvement, width,
                       label='涉及有机物', color='#95E1D3', alpha=0.8)

        ax4.set_xlabel('反应模板')
        ax4.set_ylabel('是否涉及 (1=是, 0=否)')
        ax4.set_title('模板化学性质分析', fontweight='bold')
        ax4.set_xticks(x)
        ax4.set_xticklabels(templates, rotation=45, ha='right')
        ax4.legend()
        ax4.grid(True, alpha=0.3)

        # 5. 模板重要性矩阵
        ax5 = fig.add_subplot(gs[1, 2])

        # 创建重要性矩阵
        importance_matrix = []
        metrics = ['反应数量', '平均|ΔG|', '平均Ea', '涉及金属', '涉及有机物']

        for template in templates:
            analysis = self.template_types[template]
            row = [
                analysis['count'],
                abs(analysis['avg_energy']),
                analysis['avg_barrier'],
                1 if analysis['involves_metal'] else 0,
                1 if analysis['involves_organic'] else 0
            ]
            importance_matrix.append(row)

        # 标准化
        importance_matrix = np.array(importance_matrix)
        for j in range(importance_matrix.shape[1]):
            col = importance_matrix[:, j]
            if col.max() > col.min():
                importance_matrix[:, j] = (col - col.min()) / (col.max() - col.min())

        im = ax5.imshow(importance_matrix.T, cmap='YlOrRd', aspect='auto')
        ax5.set_xticks(range(len(templates)))
        ax5.set_xticklabels(templates, rotation=45, ha='right')
        ax5.set_yticks(range(len(metrics)))
        ax5.set_yticklabels(metrics)
        ax5.set_title('模板重要性热图', fontweight='bold')

        # 添加数值标注
        for i in range(len(metrics)):
            for j in range(len(templates)):
                text = ax5.text(j, i, f'{importance_matrix[j, i]:.2f}',
                               ha="center", va="center", color="black", fontsize=8)

        plt.colorbar(im, ax=ax5, shrink=0.8)

        # 6. 主导模板识别
        ax6 = fig.add_subplot(gs[2, :])

        # 计算综合重要性分数
        template_scores = {}
        for template in templates:
            analysis = self.template_types[template]

            # 综合分数 = 反应数量权重 + 能量权重 + 化学性质权重
            count_score = analysis['count'] / max(counts) * 3
            energy_score = max(0, -analysis['avg_energy'] / 50.0) * 2  # 放热反应更重要
            barrier_score = max(0, (50 - analysis['avg_barrier']) / 50.0) * 1  # 低活化能更重要
            chem_score = (analysis['involves_metal'] + analysis['involves_organic']) * 1

            total_score = count_score + energy_score + barrier_score + chem_score
            template_scores[template] = total_score

        # 排序并绘制
        sorted_templates = sorted(template_scores.items(), key=lambda x: x[1], reverse=True)
        template_names = [t[0] for t in sorted_templates]
        scores = [t[1] for t in sorted_templates]
        template_colors_sorted = [self.template_colors[t] for t in template_names]

        bars = ax6.bar(range(len(template_names)), scores,
                      color=template_colors_sorted, alpha=0.8, edgecolor='black')
        ax6.set_xticks(range(len(template_names)))
        ax6.set_xticklabels(template_names, rotation=45, ha='right')
        ax6.set_ylabel('综合重要性分数')
        ax6.set_title('主导反应模板排名 (分数 = 频次×3 + 放热性×2 + 易反应性×1 + 化学性质×1)',
                     fontweight='bold')
        ax6.grid(True, alpha=0.3)

        # 添加分数标签
        for bar, score in zip(bars, scores):
            height = bar.get_height()
            ax6.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                    f'{score:.1f}', ha='center', va='bottom', fontweight='bold')

        # 高亮前3名
        for i in range(min(3, len(bars))):
            bars[i].set_edgecolor('red')
            bars[i].set_linewidth(3)

        plt.suptitle('反应模板分析仪表板', fontsize=18, fontweight='bold')

        output_file = self.output_dir / 'template_analysis_dashboard.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✅ 模板分析仪表板已保存: {output_file}")
        plt.close()

    def create_pathway_analysis(self):
        """创建路径分析"""

        if not self.dominant_pathways:
            print("⚠️  未找到主导路径")
            return

        fig, axes = plt.subplots(len(self.dominant_pathways), 1,
                                figsize=(18, 5 * len(self.dominant_pathways)))

        if len(self.dominant_pathways) == 1:
            axes = [axes]

        for idx, pathway in enumerate(self.dominant_pathways):
            ax = axes[idx]
            path = pathway['path']
            start = pathway['start']
            end = pathway['end']
            total_energy = pathway['total_energy']

            # 绘制路径能量剖面
            x_positions = list(range(len(path)))
            y_positions = []
            colors = []
            shapes = []
            labels = []

            cumulative_energy = 0

            for i, node in enumerate(path):
                if node.startswith('R'):  # 反应节点
                    rxn_idx = int(node[1:])
                    rxn = self.generator.all_reactions[rxn_idx]

                    # 反应的y坐标 = 前一个物种能量 + 活化能
                    barrier_energy = cumulative_energy + rxn.activation_energy / 10.0
                    y_positions.append(barrier_energy)

                    # 反应节点颜色和标签
                    template = rxn.name
                    colors.append(self.template_colors.get(template, '#CCCCCC'))
                    shapes.append('s')
                    labels.append(f"{template}\nΔG={rxn.reaction_energy:.1f}\nEa={rxn.activation_energy:.1f}")

                    # 更新累积能量
                    cumulative_energy += rxn.reaction_energy / 10.0

                else:  # 物种节点
                    y_positions.append(cumulative_energy)

                    # 物种节点颜色
                    mol = None
                    for smiles, (m, gen) in self.generator.all_molecules.items():
                        if m.name == node:
                            mol = m
                            break

                    if mol:
                        if mol in self.key_species['initial']:
                            colors.append('#FF6B6B')
                        elif mol in self.key_species['final']:
                            colors.append('#4ECDC4')
                        elif mol in self.key_species['metal_containing']:
                            colors.append('#FFD93D')
                        else:
                            colors.append('#95E1D3')
                    else:
                        colors.append('#CCCCCC')

                    shapes.append('o')
                    labels.append(node)

            # 绘制节点
            for i, (x, y, color, shape, label) in enumerate(zip(x_positions, y_positions, colors, shapes, labels)):
                if shape == 'o':  # 物种
                    ax.scatter(x, y, s=500, c=color, marker='o',
                              edgecolors='black', linewidths=2, alpha=0.8, zorder=3)
                else:  # 反应
                    ax.scatter(x, y, s=400, c=color, marker='s',
                              edgecolors='black', linewidths=2, alpha=0.9, zorder=3)

                # 添加标签
                ax.annotate(label, (x, y), xytext=(0, 20), textcoords='offset points',
                           ha='center', va='bottom', fontsize=9, fontweight='bold',
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

            # 绘制连线
            ax.plot(x_positions, y_positions, 'k-', linewidth=3, alpha=0.7, zorder=1)

            ax.set_title(f'路径 {idx+1}: {start} → {end} (总ΔG = {total_energy:.1f} kcal/mol)',
                        fontsize=14, fontweight='bold')
            ax.set_xlabel('反应步骤', fontsize=12)
            ax.set_ylabel('相对能量 (kcal/mol)', fontsize=12)
            ax.grid(True, alpha=0.3)
            ax.axhline(y=0, color='red', linestyle='--', alpha=0.5, label='能量基准')
            ax.legend()

        plt.suptitle('主导反应路径分析', fontsize=16, fontweight='bold')
        plt.tight_layout()

        output_file = self.output_dir / 'universal_pathways.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✅ 通用路径分析图已保存: {output_file}")
        plt.close()

    def generate_summary_report(self):
        """生成分析摘要报告"""

        report_lines = []
        report_lines.append("=" * 80)
        report_lines.append("通用反应网络分析报告")
        report_lines.append("=" * 80)

        # 基本统计
        report_lines.append(f"\n📊 基本统计:")
        report_lines.append(f"   总分子数: {len(self.generator.all_molecules)}")
        report_lines.append(f"   总反应数: {len(self.generator.all_reactions)}")
        report_lines.append(f"   反应模板类型数: {len(self.template_types)}")

        # 关键物种分析
        report_lines.append(f"\n🔬 关键物种分析:")
        report_lines.append(f"   初始物种: {[m.name for m in self.key_species['initial']]}")
        report_lines.append(f"   最终产物: {[m.name for m in self.key_species['final']]}")
        report_lines.append(f"   重要中间体: {[m.name for m in self.key_species['intermediate']]}")
        report_lines.append(f"   含金属物种: {[m.name for m in self.key_species['metal_containing']]}")

        # 主导反应模板
        template_scores = {}
        for template in self.template_types:
            analysis = self.template_types[template]
            count_score = analysis['count'] / len(self.generator.all_reactions) * 100
            energy_score = max(0, -analysis['avg_energy'])
            template_scores[template] = count_score + energy_score

        sorted_templates = sorted(template_scores.items(), key=lambda x: x[1], reverse=True)

        report_lines.append(f"\n⚡ 主导反应模板 (前5名):")
        for i, (template, score) in enumerate(sorted_templates[:5]):
            analysis = self.template_types[template]
            report_lines.append(f"   {i+1}. {template}")
            report_lines.append(f"      反应数量: {analysis['count']} ({analysis['count']/len(self.generator.all_reactions)*100:.1f}%)")
            report_lines.append(f"      平均ΔG: {analysis['avg_energy']:.1f} kcal/mol")
            report_lines.append(f"      平均Ea: {analysis['avg_barrier']:.1f} kcal/mol")
            report_lines.append(f"      涉及金属: {'是' if analysis['involves_metal'] else '否'}")
            report_lines.append(f"      涉及有机物: {'是' if analysis['involves_organic'] else '否'}")

        # 主导路径
        if self.dominant_pathways:
            report_lines.append(f"\n🛤️  主导反应路径 (前3条):")
            for i, pathway in enumerate(self.dominant_pathways[:3]):
                report_lines.append(f"   路径 {i+1}: {pathway['start']} → {pathway['end']}")
                report_lines.append(f"      总能量变化: {pathway['total_energy']:.1f} kcal/mol")
                report_lines.append(f"      反应步数: {len(pathway['reactions'])}")
                report_lines.append(f"      涉及模板: {[self.generator.all_reactions[int(r[1:])].name for r in pathway['reactions']]}")

        report_lines.append(f"\n" + "=" * 80)

        # 保存报告
        report_content = "\n".join(report_lines)
        report_file = self.output_dir / 'analysis_report.txt'
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(report_content)

        print(f"✅ 分析报告已保存: {report_file}")
        print("\n" + report_content)

    def generate_all_visualizations(self):
        """生成所有可视化"""

        print(f"\n{'='*80}")
        print(f"通用反应网络可视化分析")
        print(f"{'='*80}")

        print(f"\n🔍 自动识别的反应模板类型:")
        for template, analysis in self.template_types.items():
            print(f"   - {template}: {analysis['count']} 个反应")

        print(f"\n🎨 生成可视化图表...")

        # 1. 全局网络视图
        self.create_global_network_view()

        # 2. 模板分析仪表板
        self.create_template_analysis_dashboard()

        # 3. 路径分析
        self.create_pathway_analysis()

        # 4. 生成摘要报告
        self.generate_summary_report()

        print(f"\n{'='*80}")
        print(f"✅ 通用可视化分析完成!")
        print(f"   输出目录: {self.output_dir.absolute()}")
        print(f"   - universal_global_network.png (全局网络)")
        print(f"   - template_analysis_dashboard.png (模板分析仪表板)")
        print(f"   - universal_pathways.png (主导路径)")
        print(f"   - analysis_report.txt (分析报告)")
        print(f"{'='*80}")


def main():
    """主函数"""
    print("通用反应网络可视化演示")
    print("=" * 80)

    # 1. 生成网络
    print("\n生成反应网络...")

    li_ion = Molecule.from_smiles('[Li+]', name='Li_ion')
    ec = Molecule.from_smiles('C1COC(=O)O1', name='EC')
    seed_molecules = [li_ion, ec]

    env = Environment(
        temperature=300.0,
        electrode_type='anode',
        voltage=0.05,
        li_activity=1.0
    )

    generator = IterativeNetworkGenerator(
        max_generations=3,
        max_molecules=30
    )

    stats = generator.generate_network(seed_molecules, env)

    print(f"\n网络生成完成:")
    print(f"  总分子数: {stats['total_molecules']}")
    print(f"  总反应数: {stats['total_reactions']}")
    print(f"  最大代数: {stats['max_generation']}")

    # 2. 通用可视化
    output_dir = Path('universal_network_output')
    visualizer = UniversalNetworkVisualizer(generator, output_dir)
    visualizer.generate_all_visualizations()

    return 0


if __name__ == '__main__':
    exit(main())
