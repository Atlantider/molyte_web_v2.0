#!/usr/bin/env python3
"""
RSNet Simple API - 简洁的反应网络生成接口
输入SMILES，自动生成反应网络、计算能量、可视化

用法:
    from rsnet_simple_api import generate_reaction_network
    
    result = generate_reaction_network(
        smiles_list=['C1COC(=O)O1', '[Li+]', 'F[P-](F)(F)(F)(F)F'],
        temperature=300.0,
        max_generations=3
    )
    
    # 结果包含：
    # - result['network']: 反应网络对象
    # - result['reactions']: 所有反应列表
    # - result['molecules']: 所有分子列表
    # - result['visualization']: 可视化文件路径
    # - result['statistics']: 统计信息
"""

import os
import json
import time
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import matplotlib.pyplot as plt
import networkx as nx
from datetime import datetime

# 导入rsnet核心模块
from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.network.generator import NetworkGenerator, ReactionNetwork
from rsnet.network.config import NetworkGenerationConfig
from rsnet.operators.registry import OPERATOR_REGISTRY
from rsnet.compute.reaction_screener import ReactionScreener


class RSNetSimpleAPI:
    """简化的RSNet API接口"""
    
    def __init__(self, output_dir: str = "./rsnet_output"):
        """
        初始化API
        
        Args:
            output_dir: 输出目录，用于保存可视化和结果文件
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def generate_reaction_network(
        self,
        smiles_list: List[str],
        temperature: float = 300.0,
        electrode_type: str = 'anode',
        voltage: float = 0.1,
        max_generations: int = 3,
        max_species: int = 50,
        energy_cutoff: float = 80.0,
        visualize: bool = True,
        save_results: bool = True
    ) -> Dict[str, Any]:
        """
        生成反应网络的主要接口
        
        Args:
            smiles_list: SMILES字符串列表
            temperature: 温度 (K)
            electrode_type: 电极类型 ('anode' 或 'cathode')
            voltage: 电压 (V)
            max_generations: 最大代数
            max_species: 最大分子数
            energy_cutoff: 能量截断值 (kcal/mol)
            visualize: 是否生成可视化
            save_results: 是否保存结果到文件
            
        Returns:
            结果字典，包含网络、反应、分子、统计信息和可视化路径
        """
        print(f"\n{'='*60}")
        print(f"RSNet 反应网络生成")
        print(f"{'='*60}")
        print(f"输入分子数: {len(smiles_list)}")
        print(f"温度: {temperature} K")
        print(f"电极类型: {electrode_type}")
        print(f"电压: {voltage} V")
        print(f"最大代数: {max_generations}")
        print(f"{'='*60}\n")
        
        start_time = time.time()
        
        # 1. 创建分子对象
        print("步骤 1/5: 创建分子对象...")
        molecules = []
        for i, smiles in enumerate(smiles_list):
            try:
                mol = Molecule.from_smiles(smiles, name=f"mol_{i+1}")
                molecules.append(mol)
                print(f"  ✓ {mol.name}: {smiles}")
            except Exception as e:
                print(f"  ✗ 无法解析SMILES: {smiles} - {e}")
        
        if not molecules:
            raise ValueError("没有有效的分子，无法生成网络")
        
        # 2. 创建环境
        print("\n步骤 2/5: 设置反应环境...")
        environment = Environment(
            temperature=temperature,
            electrode_type=electrode_type,
            voltage=voltage,
            li_activity=1.0 if electrode_type == 'anode' else 0.1
        )
        print(f"  ✓ 环境已配置")
        
        # 3. 配置网络生成器
        print("\n步骤 3/5: 配置网络生成器...")
        config = NetworkGenerationConfig(
            max_generations=max_generations,
            max_species=max_species,
            energy_cutoff=energy_cutoff,
            use_structure_based_filtering=True,
            driving_force_threshold=0.2
        )
        
        generator = NetworkGenerator(
            operator_registry=OPERATOR_REGISTRY,
            config=config
        )
        print(f"  ✓ 生成器已配置")
        
        # 4. 生成反应网络
        print("\n步骤 4/5: 生成反应网络...")
        print("  (这可能需要一些时间...)")
        
        network = generator.generate_network(
            seed_molecules=molecules,
            environment=environment,
            max_time=300.0  # 5分钟超时
        )
        
        generation_time = time.time() - start_time
        
        # 5. 收集结果
        print("\n步骤 5/5: 收集和分析结果...")
        stats = network.get_statistics()
        
        print(f"\n{'='*60}")
        print(f"网络生成完成！")
        print(f"{'='*60}")
        print(f"总分子数: {stats['num_molecules']}")
        print(f"总反应数: {stats['num_reactions']}")
        print(f"最大代数: {stats['max_generation']}")
        print(f"生成时间: {generation_time:.2f} 秒")
        print(f"{'='*60}\n")
        
        # 准备返回结果
        result = {
            'network': network,
            'molecules': list(network.molecules.values()),
            'reactions': list(network.reactions.values()),
            'statistics': stats,
            'generation_time': generation_time,
            'environment': {
                'temperature': temperature,
                'electrode_type': electrode_type,
                'voltage': voltage
            }
        }
        
        # 6. 可视化
        if visualize:
            print("生成可视化...")
            viz_path = self._visualize_network(network, result)
            result['visualization_path'] = str(viz_path)
            print(f"  ✓ 可视化已保存: {viz_path}")
        
        # 7. 保存结果
        if save_results:
            print("\n保存结果...")
            json_path = self._save_results(result)
            result['json_path'] = str(json_path)
            print(f"  ✓ 结果已保存: {json_path}")
        
        return result
    
    def _visualize_network(
        self, 
        network: ReactionNetwork, 
        result: Dict[str, Any]
    ) -> Path:
        """
        生成网络可视化
        
        Args:
            network: 反应网络对象
            result: 结果字典
            
        Returns:
            可视化文件路径
        """
        # 创建图形
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        
        # 左图: 网络拓扑
        G = nx.DiGraph()
        
        # 添加节点（分子）
        for smiles, mol in network.molecules.items():
            G.add_node(smiles, label=mol.name)
        
        # 添加边（反应）
        for rxn in network.reactions.values():
            for reactant in rxn.reactants:
                for product in rxn.products:
                    G.add_edge(reactant.smiles, product.smiles)
        
        # 绘制网络
        pos = nx.spring_layout(G, k=2, iterations=50)
        
        # 按代数着色
        node_colors = []
        for node in G.nodes():
            gen = network.molecules[node].generation if hasattr(network.molecules[node], 'generation') else 0
            node_colors.append(gen)
        
        nx.draw_networkx_nodes(
            G, pos, 
            node_color=node_colors,
            node_size=500,
            cmap=plt.cm.viridis,
            alpha=0.8,
            ax=ax1
        )
        
        nx.draw_networkx_edges(
            G, pos,
            edge_color='gray',
            arrows=True,
            arrowsize=20,
            alpha=0.5,
            ax=ax1
        )
        
        # 添加标签（简化的分子名称）
        labels = {node: network.molecules[node].name for node in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels, font_size=8, ax=ax1)
        
        ax1.set_title(f"反应网络拓扑\n({len(network.molecules)} 分子, {len(network.reactions)} 反应)", 
                     fontsize=14, fontweight='bold')
        ax1.axis('off')
        
        # 右图: 统计信息
        ax2.axis('off')
        
        stats = result['statistics']
        env = result['environment']
        
        stats_text = f"""
反应网络统计信息
{'='*40}

网络规模:
  • 总分子数: {stats['num_molecules']}
  • 总反应数: {stats['num_reactions']}
  • 最大代数: {stats['max_generation']}

环境条件:
  • 温度: {env['temperature']} K
  • 电极类型: {env['electrode_type']}
  • 电压: {env['voltage']} V

性能:
  • 生成时间: {result['generation_time']:.2f} 秒
  • 平均每代时间: {result['generation_time']/max(1, stats['max_generation']):.2f} 秒

反应类型分布:
"""
        
        # 统计反应类型
        reaction_types = {}
        for rxn in result['reactions']:
            op_name = rxn.operator_name if hasattr(rxn, 'operator_name') else 'unknown'
            reaction_types[op_name] = reaction_types.get(op_name, 0) + 1
        
        for op_name, count in sorted(reaction_types.items(), key=lambda x: x[1], reverse=True)[:5]:
            stats_text += f"  • {op_name}: {count}\n"
        
        ax2.text(0.1, 0.5, stats_text, 
                fontsize=11, 
                verticalalignment='center',
                family='monospace')
        
        plt.tight_layout()
        
        # 保存
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        viz_path = self.output_dir / f"network_{timestamp}.png"
        plt.savefig(viz_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        return viz_path
    
    def _save_results(self, result: Dict[str, Any]) -> Path:
        """
        保存结果到JSON文件
        
        Args:
            result: 结果字典
            
        Returns:
            JSON文件路径
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        json_path = self.output_dir / f"results_{timestamp}.json"
        
        # 准备可序列化的数据
        serializable_result = {
            'statistics': result['statistics'],
            'generation_time': result['generation_time'],
            'environment': result['environment'],
            'molecules': [
                {
                    'name': mol.name,
                    'smiles': mol.smiles,
                    'generation': mol.generation if hasattr(mol, 'generation') else 0
                }
                for mol in result['molecules']
            ],
            'reactions': [
                {
                    'reactants': [r.smiles for r in rxn.reactants],
                    'products': [p.smiles for p in rxn.products],
                    'operator': rxn.operator_name if hasattr(rxn, 'operator_name') else 'unknown',
                    'energy': rxn.reaction_energy if hasattr(rxn, 'reaction_energy') else None
                }
                for rxn in result['reactions']
            ]
        }
        
        if 'visualization_path' in result:
            serializable_result['visualization_path'] = result['visualization_path']
        
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(serializable_result, f, indent=2, ensure_ascii=False)
        
        return json_path


# 便捷函数
def generate_reaction_network(
    smiles_list: List[str],
    temperature: float = 300.0,
    electrode_type: str = 'anode',
    voltage: float = 0.1,
    max_generations: int = 3,
    max_species: int = 50,
    output_dir: str = "./rsnet_output",
    **kwargs
) -> Dict[str, Any]:
    """
    便捷函数：生成反应网络
    
    Args:
        smiles_list: SMILES字符串列表
        temperature: 温度 (K)
        electrode_type: 电极类型
        voltage: 电压 (V)
        max_generations: 最大代数
        max_species: 最大分子数
        output_dir: 输出目录
        **kwargs: 其他参数
        
    Returns:
        结果字典
    """
    api = RSNetSimpleAPI(output_dir=output_dir)
    return api.generate_reaction_network(
        smiles_list=smiles_list,
        temperature=temperature,
        electrode_type=electrode_type,
        voltage=voltage,
        max_generations=max_generations,
        max_species=max_species,
        **kwargs
    )


# 命令行接口
def main():
    """命令行接口"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='RSNet Simple API - 反应网络生成器',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 基本用法
  python rsnet_simple_api.py --smiles "C1COC(=O)O1" "[Li+]" "F[P-](F)(F)(F)(F)F"
  
  # 指定参数
  python rsnet_simple_api.py --smiles "CCO" "C=C" --temperature 400 --max-generations 4
  
  # 从文件读取
  python rsnet_simple_api.py --input molecules.txt --output ./my_results
        """
    )
    
    parser.add_argument('--smiles', nargs='+', help='SMILES字符串列表')
    parser.add_argument('--input', '-i', help='包含SMILES的文本文件（每行一个）')
    parser.add_argument('--output', '-o', default='./rsnet_output', help='输出目录')
    parser.add_argument('--temperature', '-T', type=float, default=300.0, help='温度 (K)')
    parser.add_argument('--electrode', choices=['anode', 'cathode'], default='anode', help='电极类型')
    parser.add_argument('--voltage', '-V', type=float, default=0.1, help='电压 (V)')
    parser.add_argument('--max-generations', '-g', type=int, default=3, help='最大代数')
    parser.add_argument('--max-species', '-s', type=int, default=50, help='最大分子数')
    parser.add_argument('--no-viz', action='store_true', help='不生成可视化')
    
    args = parser.parse_args()
    
    # 获取SMILES列表
    if args.smiles:
        smiles_list = args.smiles
    elif args.input:
        with open(args.input, 'r') as f:
            smiles_list = [line.strip() for line in f if line.strip()]
    else:
        parser.print_help()
        print("\n错误: 必须提供 --smiles 或 --input 参数")
        return 1
    
    # 生成网络
    try:
        result = generate_reaction_network(
            smiles_list=smiles_list,
            temperature=args.temperature,
            electrode_type=args.electrode,
            voltage=args.voltage,
            max_generations=args.max_generations,
            max_species=args.max_species,
            output_dir=args.output,
            visualize=not args.no_viz
        )
        
        print(f"\n✅ 成功！")
        if 'visualization_path' in result:
            print(f"📊 可视化: {result['visualization_path']}")
        if 'json_path' in result:
            print(f"💾 结果: {result['json_path']}")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ 错误: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    exit(main())
