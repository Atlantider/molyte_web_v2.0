#!/usr/bin/env python3
"""
多分子SMILES处理演示

展示如何同时处理多个分子的SMILES，生成复杂的反应网络。
包括：
1. 多分子输入处理
2. 混合反应网络生成
3. 分子间反应探索
4. 复杂网络可视化
"""

import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # 使用无图形后端

# 添加rsnet到路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.network.generator import NetworkGenerator
from rsnet.network.config import NetworkGenerationConfig, NetworkGenerationPresets
from rsnet.utils.visualization import (
    plot_network_graph, 
    plot_energy_distribution,
    plot_generation_statistics
)
from rsnet.utils.io import (
    create_analysis_report,
    save_molecules_to_smiles,
    load_molecules_from_smiles
)


def demo_multiple_smiles_input():
    """演示多种SMILES输入方式"""
    print("=== 多分子SMILES输入演示 ===")
    
    # 方式1: 直接列表输入
    smiles_list = [
        'CCO',           # 乙醇
        'CC',            # 乙烷  
        'C=C',           # 乙烯
        'C#C',           # 乙炔
        'CC(=O)O',       # 乙酸
        'CCO[Li]',       # 锂乙醇络合物
        'C1CC1',         # 环丙烷
        'c1ccccc1'       # 苯
    ]
    
    molecules = []
    names = ['ethanol', 'ethane', 'ethylene', 'acetylene', 'acetic_acid', 'li_ethanol', 'cyclopropane', 'benzene']
    
    for i, smiles in enumerate(smiles_list):
        try:
            mol = Molecule.from_smiles(smiles, name=names[i])
            molecules.append(mol)
            print(f"✓ 加载分子: {mol.name} ({mol.smiles})")
        except Exception as e:
            print(f"✗ 无法加载 {smiles}: {e}")
    
    print(f"成功加载 {len(molecules)} 个分子")
    return molecules


def demo_file_based_input():
    """演示基于文件的多分子输入"""
    print("\n=== 文件输入演示 ===")
    
    # 创建示例分子文件
    molecules_data = [
        {'smiles': 'CCO', 'name': 'ethanol'},
        {'smiles': 'CC(=O)O', 'name': 'acetic_acid'},
        {'smiles': 'C=C', 'name': 'ethylene'},
        {'smiles': 'CC', 'name': 'ethane'},
        {'smiles': 'C1CC1', 'name': 'cyclopropane'},
        {'smiles': 'CCO[Li]', 'name': 'li_ethanol_complex'},
        {'smiles': 'CC(C)O', 'name': 'isopropanol'},
        {'smiles': 'C=CC=C', 'name': 'butadiene'}
    ]
    
    # 保存到CSV文件
    import pandas as pd
    df = pd.DataFrame(molecules_data)
    csv_file = 'example_molecules.csv'
    df.to_csv(csv_file, index=False)
    print(f"创建分子文件: {csv_file}")
    
    # 从文件加载
    molecules = load_molecules_from_smiles(csv_file, name_column='name')
    print(f"从文件加载了 {len(molecules)} 个分子:")
    for mol in molecules:
        print(f"  - {mol.name}: {mol.smiles}")
    
    return molecules


def demo_mixed_reaction_network():
    """演示多分子混合反应网络生成"""
    print("\n=== 混合反应网络生成演示 ===")
    
    # 选择有代表性的分子组合
    mixed_smiles = [
        'CCO',           # 乙醇 - 有机溶剂
        'C=C',           # 乙烯 - 不饱和烃
        'CC(=O)O',       # 乙酸 - 有机酸
        'C1CC1',         # 环丙烷 - 环状化合物
        'CCO[Li]'        # 锂络合物 - 电池相关
    ]
    
    molecules = []
    names = ['ethanol', 'ethylene', 'acetic_acid', 'cyclopropane', 'li_ethanol']
    
    for i, smiles in enumerate(mixed_smiles):
        mol = Molecule.from_smiles(smiles, name=names[i])
        molecules.append(mol)
    
    print(f"混合分子体系包含 {len(molecules)} 个分子:")
    for mol in molecules:
        print(f"  - {mol.name}: {mol.smiles}")
    
    # 设置复杂环境 - 模拟电池电解液环境
    env = Environment(
        temperature=450.0,      # 中等温度
        electrode_type='cathode',
        voltage=4.0,
        li_activity=0.4,
        interface_type='SEI'
    )
    
    print(f"\n环境条件:")
    print(f"  温度: {env.temperature} K")
    print(f"  电极: {env.electrode_type}")
    print(f"  电压: {env.voltage} V")
    print(f"  Li活度: {env.li_activity}")
    
    # 使用智能配置生成网络
    config = NetworkGenerationConfig(
        max_generations=3,
        max_species=30,
        energy_cutoff=55.0,
        generation_strategy='intelligent',
        operator_selection_strategy='adaptive',
        driving_force_threshold=0.3,
        max_operators_per_generation=6
    )
    
    print(f"\n网络生成配置:")
    print(f"  最大代数: {config.max_generations}")
    print(f"  最大分子数: {config.max_species}")
    print(f"  能量阈值: {config.energy_cutoff} kcal/mol")
    print(f"  生成策略: {config.generation_strategy.value}")
    
    # 生成网络
    generator = NetworkGenerator(config=config)
    print(f"\n开始生成混合反应网络...")
    
    network = generator.generate_network(molecules, env, max_time=45.0)
    
    # 显示结果
    stats = network.get_statistics()
    print(f"\n网络生成完成!")
    print(f"  最终分子数: {stats['num_molecules']}")
    print(f"  反应数: {stats['num_reactions']}")
    print(f"  网络边数: {stats['num_edges']}")
    print(f"  最大代数: {stats['max_generation']}")
    
    # 按代数显示分子分布
    print(f"\n各代分子分布:")
    for gen, count in stats['molecules_by_generation'].items():
        print(f"  第{gen}代: {count} 个分子")
    
    return network, molecules


def demo_intermolecular_reactions():
    """演示分子间反应探索"""
    print("\n=== 分子间反应探索演示 ===")
    
    # 选择容易发生分子间反应的分子
    reactive_smiles = [
        'C=C',           # 乙烯 - 可加成
        'CC(=O)O',       # 乙酸 - 可质子化
        'CCO',           # 乙醇 - 可脱水
        '[Li+]',         # 锂离子 - 可配位
        'C1CC1'          # 环丙烷 - 可开环
    ]
    
    molecules = []
    names = ['ethylene', 'acetic_acid', 'ethanol', 'lithium_ion', 'cyclopropane']
    
    for i, smiles in enumerate(reactive_smiles):
        try:
            mol = Molecule.from_smiles(smiles, name=names[i])
            molecules.append(mol)
        except:
            print(f"跳过无效SMILES: {smiles}")
    
    print(f"反应性分子体系: {len(molecules)} 个分子")
    
    # 高反应性环境
    env = Environment(
        temperature=600.0,      # 高温促进反应
        electrode_type='anode',
        voltage=0.5,
        li_activity=0.8,        # 高Li活度
        interface_type='bulk'
    )
    
    # 专注于分子间反应的配置
    config = NetworkGenerationConfig(
        max_generations=2,
        max_species=20,
        energy_cutoff=50.0,
        generation_strategy='focused',
        operator_selection_strategy='adaptive',
        driving_force_threshold=0.4,  # 较高阈值确保反应活性
        max_operators_per_generation=8
    )
    
    generator = NetworkGenerator(config=config)
    network = generator.generate_network(molecules, env, max_time=30.0)
    
    stats = network.get_statistics()
    print(f"分子间反应网络: {stats['num_molecules']} 分子, {stats['num_reactions']} 反应")
    
    return network


def demo_complex_visualization():
    """演示复杂网络可视化"""
    print("\n=== 复杂网络可视化演示 ===")
    
    # 生成一个较大的网络用于可视化
    molecules = demo_multiple_smiles_input()[:6]  # 取前6个分子
    
    env = Environment(temperature=500.0, electrode_type='cathode', voltage=4.2)
    config = NetworkGenerationConfig(
        max_generations=3,
        max_species=25,
        energy_cutoff=60.0
    )
    
    generator = NetworkGenerator(config=config)
    network = generator.generate_network(molecules, env, max_time=40.0)
    
    stats = network.get_statistics()
    print(f"可视化网络: {stats['num_molecules']} 分子, {stats['num_reactions']} 反应")
    
    if stats['num_molecules'] > 1:
        # 创建多种可视化
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # 网络图 - Spring布局
        plt.sca(axes[0, 0])
        plot_network_graph(network, layout='spring', figsize=(8, 6))
        plt.title('Spring Layout Network', fontsize=14, fontweight='bold')
        
        # 网络图 - Circular布局  
        plt.sca(axes[0, 1])
        plot_network_graph(network, layout='circular', figsize=(8, 6))
        plt.title('Circular Layout Network', fontsize=14, fontweight='bold')
        
        # 能量分布
        if stats['num_reactions'] > 0:
            plt.sca(axes[1, 0])
            plot_energy_distribution(network, figsize=(8, 6))
            plt.title('Energy Distribution', fontsize=14, fontweight='bold')
        else:
            axes[1, 0].text(0.5, 0.5, 'No Energy Data', ha='center', va='center')
            axes[1, 0].set_title('Energy Distribution', fontsize=14, fontweight='bold')
        
        # 生成统计
        plt.sca(axes[1, 1])
        plot_generation_statistics(network, figsize=(8, 6))
        plt.title('Generation Statistics', fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig('multi_molecule_network.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("复杂网络可视化保存为: multi_molecule_network.png")
    else:
        print("网络太小，跳过可视化")
    
    return network


def demo_batch_processing():
    """演示批量处理多个分子集合"""
    print("\n=== 批量处理演示 ===")
    
    # 定义不同的分子集合
    molecule_sets = {
        'alcohols': ['CCO', 'CC(C)O', 'CCC(C)O', 'CC(C)(C)O'],
        'alkenes': ['C=C', 'CC=C', 'C=CC=C', 'CC=CC'],
        'aromatics': ['c1ccccc1', 'Cc1ccccc1', 'c1ccc(O)cc1'],
        'li_complexes': ['CCO[Li]', 'CC(=O)O[Li]', '[Li+]']
    }
    
    results = {}
    
    for set_name, smiles_list in molecule_sets.items():
        print(f"\n处理分子集合: {set_name}")
        
        # 加载分子
        molecules = []
        for i, smiles in enumerate(smiles_list):
            try:
                mol = Molecule.from_smiles(smiles, name=f"{set_name}_{i}")
                molecules.append(mol)
            except:
                print(f"  跳过无效SMILES: {smiles}")
        
        if not molecules:
            print(f"  {set_name} 集合无有效分子，跳过")
            continue
            
        print(f"  加载了 {len(molecules)} 个分子")
        
        # 根据分子类型选择环境
        if 'li' in set_name:
            env = Environment(temperature=400.0, electrode_type='cathode', voltage=4.0, li_activity=0.6)
        elif 'aromatic' in set_name:
            env = Environment(temperature=600.0, pressure=1.0)
        else:
            env = Environment(temperature=450.0, pressure=1.0)
        
        # 快速配置
        config = NetworkGenerationPresets.fast_screening()
        config.max_species = 15
        config.max_generations = 2
        
        # 生成网络
        generator = NetworkGenerator(config=config)
        network = generator.generate_network(molecules, env, max_time=20.0)
        
        stats = network.get_statistics()
        results[set_name] = stats
        
        print(f"  结果: {stats['num_molecules']} 分子, {stats['num_reactions']} 反应")
    
    # 总结批量处理结果
    print(f"\n=== 批量处理总结 ===")
    total_molecules = sum(r['num_molecules'] for r in results.values())
    total_reactions = sum(r['num_reactions'] for r in results.values())
    
    print(f"处理了 {len(results)} 个分子集合")
    print(f"总计生成: {total_molecules} 个分子, {total_reactions} 个反应")
    
    for set_name, stats in results.items():
        print(f"  {set_name}: {stats['num_molecules']} 分子, {stats['num_reactions']} 反应")
    
    return results


def demo_comprehensive_analysis():
    """演示多分子体系的综合分析"""
    print("\n=== 综合分析演示 ===")
    
    # 创建一个复杂的多分子体系
    complex_smiles = [
        'CCO',           # 乙醇
        'CC(=O)O',       # 乙酸  
        'C=C',           # 乙烯
        'CC',            # 乙烷
        'C1CC1',         # 环丙烷
        'CCO[Li]',       # 锂络合物
        'CC(C)O',        # 异丙醇
        'C=CC=C'         # 丁二烯
    ]
    
    molecules = []
    for i, smiles in enumerate(complex_smiles):
        try:
            mol = Molecule.from_smiles(smiles, name=f'mol_{i+1}')
            molecules.append(mol)
        except:
            continue
    
    print(f"复杂体系包含 {len(molecules)} 个分子")
    
    # 电池电解液环境
    env = Environment(
        temperature=500.0,
        electrode_type='cathode', 
        voltage=4.2,
        li_activity=0.5,
        interface_type='SEI'
    )
    
    # 全面配置
    config = NetworkGenerationConfig(
        max_generations=4,
        max_species=40,
        energy_cutoff=55.0,
        generation_strategy='intelligent',
        operator_selection_strategy='adaptive',
        driving_force_threshold=0.25,
        max_operators_per_generation=8
    )
    
    # 生成网络
    generator = NetworkGenerator(config=config)
    print("开始生成复杂反应网络...")
    network = generator.generate_network(molecules, env, max_time=60.0)
    
    # 创建综合分析报告
    output_dir = Path('multi_molecule_analysis')
    create_analysis_report(network, output_dir, include_plots=True)
    
    stats = network.get_statistics()
    print(f"\n综合分析完成!")
    print(f"  网络规模: {stats['num_molecules']} 分子, {stats['num_reactions']} 反应")
    print(f"  分析报告保存在: {output_dir}")
    print(f"  包含文件:")
    
    if output_dir.exists():
        for file in output_dir.iterdir():
            if file.is_file():
                print(f"    - {file.name}")
    
    return network


def main():
    """主函数 - 运行所有多分子演示"""
    print("RSNet 多分子SMILES处理演示")
    print("=" * 50)
    
    try:
        # 创建输出目录
        output_dir = Path('multi_molecule_output')
        output_dir.mkdir(exist_ok=True)
        os.chdir(output_dir)
        
        # 运行各种演示
        molecules = demo_multiple_smiles_input()
        
        file_molecules = demo_file_based_input()
        
        mixed_network, _ = demo_mixed_reaction_network()
        
        inter_network = demo_intermolecular_reactions()
        
        vis_network = demo_complex_visualization()
        
        batch_results = demo_batch_processing()
        
        final_network = demo_comprehensive_analysis()
        
        # 最终总结
        print("\n" + "=" * 50)
        print("多分子处理演示完成!")
        print(f"输出目录: {output_dir.absolute()}")
        
        print("\n演示内容:")
        print("  ✓ 多种SMILES输入方式")
        print("  ✓ 文件批量加载")
        print("  ✓ 混合反应网络生成")
        print("  ✓ 分子间反应探索")
        print("  ✓ 复杂网络可视化")
        print("  ✓ 批量处理不同分子集合")
        print("  ✓ 综合分析报告生成")
        
        print("\n支持的输入格式:")
        print("  - Python列表: ['CCO', 'CC', 'C=C']")
        print("  - CSV文件: smiles,name 格式")
        print("  - JSON文件: [{'smiles': 'CCO', 'name': 'ethanol'}]")
        print("  - 文本文件: 每行一个SMILES")
        
        print("\n生成的文件:")
        for file in output_dir.glob('*'):
            if file.is_file():
                print(f"  - {file.name}")
            elif file.is_dir():
                print(f"  - {file.name}/ (目录)")
        
    except Exception as e:
        print(f"\n演示过程中出现错误: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == '__main__':
    exit(main())
