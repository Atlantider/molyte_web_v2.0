#!/usr/bin/env python3
"""
快速多分子SMILES处理示例

展示如何快速处理多个SMILES并生成反应网络。
"""

import sys
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

# 添加rsnet到路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.network.generator import NetworkGenerator
from rsnet.network.config import NetworkGenerationPresets
from rsnet.utils.io import load_molecules_from_smiles_list, create_analysis_report
from rsnet.utils.visualization import plot_network_graph


def quick_multi_molecule_example():
    """快速多分子处理示例"""
    print("=== 快速多分子SMILES处理 ===")
    
    # 1. 定义多个SMILES
    smiles_list = [
        'CCO',           # 乙醇
        'CC',            # 乙烷
        'C=C',           # 乙烯
        'CC(=O)O',       # 乙酸
        'C1CC1',         # 环丙烷
        'CCO[Li]'        # 锂乙醇络合物
    ]
    
    # 分子名称（可选）
    names = ['ethanol', 'ethane', 'ethylene', 'acetic_acid', 'cyclopropane', 'li_ethanol']
    
    print(f"输入SMILES: {smiles_list}")
    
    # 2. 批量加载分子
    molecules = load_molecules_from_smiles_list(smiles_list, names)
    print(f"成功加载 {len(molecules)} 个分子:")
    for mol in molecules:
        print(f"  - {mol.name}: {mol.smiles}")
    
    # 3. 设置环境（电池化学环境）
    env = Environment(
        temperature=500.0,
        electrode_type='cathode',
        voltage=4.2,
        li_activity=0.4
    )
    print(f"\n环境: {env.temperature}K, {env.electrode_type}, {env.voltage}V")
    
    # 4. 使用预设配置
    config = NetworkGenerationPresets.battery_chemistry()
    config.max_generations = 3
    config.max_species = 20
    print(f"配置: {config.max_generations}代, 最多{config.max_species}个分子")
    
    # 5. 生成网络
    generator = NetworkGenerator(config=config)
    print("\n开始生成反应网络...")
    
    network = generator.generate_network(molecules, env, max_time=30.0)
    
    # 6. 显示结果
    stats = network.get_statistics()
    print(f"\n网络生成完成!")
    print(f"  分子数: {stats['num_molecules']}")
    print(f"  反应数: {stats['num_reactions']}")
    print(f"  最大代数: {stats['max_generation']}")
    
    # 按代数显示
    print(f"\n各代分子分布:")
    for gen in sorted(stats['molecules_by_generation'].keys()):
        count = stats['molecules_by_generation'][gen]
        print(f"  第{gen}代: {count} 个分子")
    
    # 7. 可视化（如果网络足够大）
    if stats['num_molecules'] > 1:
        print(f"\n生成网络可视化...")
        fig = plot_network_graph(network, layout='spring', figsize=(10, 8))
        plt.title(f'Multi-Molecule Reaction Network\n{stats["num_molecules"]} molecules, {stats["num_reactions"]} reactions')
        plt.savefig('quick_multi_network.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("网络图保存为: quick_multi_network.png")
    
    # 8. 生成分析报告
    print(f"\n生成分析报告...")
    create_analysis_report(network, 'quick_analysis/', include_plots=True)
    print("分析报告保存在: quick_analysis/")
    
    return network


def batch_smiles_processing():
    """批量SMILES处理示例"""
    print("\n=== 批量SMILES处理 ===")
    
    # 不同类型的分子组
    molecule_groups = {
        'alcohols': {
            'smiles': ['CCO', 'CC(C)O', 'CCC(C)O'],
            'names': ['ethanol', 'isopropanol', 'sec-butanol']
        },
        'alkenes': {
            'smiles': ['C=C', 'CC=C', 'C=CC=C'],
            'names': ['ethylene', 'propene', 'butadiene']
        },
        'acids': {
            'smiles': ['CC(=O)O', 'CCC(=O)O', 'C(=O)O'],
            'names': ['acetic_acid', 'butyric_acid', 'formic_acid']
        }
    }
    
    all_results = {}
    
    for group_name, group_data in molecule_groups.items():
        print(f"\n处理分子组: {group_name}")
        
        # 加载分子
        molecules = load_molecules_from_smiles_list(
            group_data['smiles'], 
            group_data['names']
        )
        print(f"  加载了 {len(molecules)} 个分子")
        
        # 设置环境
        env = Environment(temperature=450.0, pressure=1.0)
        
        # 快速配置
        config = NetworkGenerationPresets.fast_screening()
        config.max_species = 15
        
        # 生成网络
        generator = NetworkGenerator(config=config)
        network = generator.generate_network(molecules, env, max_time=15.0)
        
        stats = network.get_statistics()
        all_results[group_name] = stats
        
        print(f"  结果: {stats['num_molecules']} 分子, {stats['num_reactions']} 反应")
    
    # 总结
    print(f"\n=== 批量处理总结 ===")
    total_molecules = sum(r['num_molecules'] for r in all_results.values())
    total_reactions = sum(r['num_reactions'] for r in all_results.values())
    
    print(f"处理了 {len(all_results)} 个分子组")
    print(f"总计: {total_molecules} 个分子, {total_reactions} 个反应")
    
    for group, stats in all_results.items():
        print(f"  {group}: {stats['num_molecules']} 分子, {stats['num_reactions']} 反应")
    
    return all_results


def main():
    """主函数"""
    print("RSNet 快速多分子处理示例")
    print("=" * 40)
    
    try:
        # 创建输出目录
        output_dir = Path('quick_multi_output')
        output_dir.mkdir(exist_ok=True)
        
        import os
        os.chdir(output_dir)
        
        # 运行示例
        network = quick_multi_molecule_example()
        
        batch_results = batch_smiles_processing()
        
        print("\n" + "=" * 40)
        print("快速多分子处理完成!")
        print(f"输出目录: {output_dir.absolute()}")
        
        print("\n使用方法总结:")
        print("1. 直接列表输入:")
        print("   smiles_list = ['CCO', 'CC', 'C=C']")
        print("   molecules = load_molecules_from_smiles_list(smiles_list)")
        
        print("\n2. 带名称输入:")
        print("   names = ['ethanol', 'ethane', 'ethylene']")
        print("   molecules = load_molecules_from_smiles_list(smiles_list, names)")
        
        print("\n3. 生成网络:")
        print("   generator = NetworkGenerator()")
        print("   network = generator.generate_network(molecules, env)")
        
        print("\n4. 分析和可视化:")
        print("   create_analysis_report(network, 'output_dir/')")
        
        print(f"\n生成的文件:")
        for file in Path('.').glob('*'):
            if file.is_file():
                print(f"  - {file.name}")
            elif file.is_dir():
                print(f"  - {file.name}/ (目录)")
        
    except Exception as e:
        print(f"\n错误: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == '__main__':
    exit(main())
