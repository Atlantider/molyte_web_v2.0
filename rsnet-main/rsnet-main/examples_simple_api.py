#!/usr/bin/env python3
"""
RSNet Simple API - 使用示例

展示如何使用简化的API生成反应网络
"""

from rsnet_simple_api import generate_reaction_network


def example_1_basic_usage():
    """示例1: 基本用法 - 锂离子电池电解液"""
    print("\n" + "="*60)
    print("示例 1: 锂离子电池电解液 (EC + Li+ + PF6-)")
    print("="*60)
    
    result = generate_reaction_network(
        smiles_list=[
            'C1COC(=O)O1',           # EC (Ethylene Carbonate)
            '[Li+]',                  # 锂离子
            'F[P-](F)(F)(F)(F)F'     # PF6- 阴离子
        ],
        temperature=300.0,
        electrode_type='anode',
        voltage=0.1,
        max_generations=2,
        max_species=30
    )
    
    print(f"\n结果摘要:")
    print(f"  - 生成了 {len(result['molecules'])} 个分子")
    print(f"  - 发现了 {len(result['reactions'])} 个反应")
    print(f"  - 用时 {result['generation_time']:.2f} 秒")
    
    return result


def example_2_organic_molecules():
    """示例2: 有机分子反应"""
    print("\n" + "="*60)
    print("示例 2: 有机分子反应 (乙醇 + 乙烯)")
    print("="*60)
    
    result = generate_reaction_network(
        smiles_list=[
            'CCO',    # 乙醇
            'C=C'     # 乙烯
        ],
        temperature=400.0,
        electrode_type='cathode',
        voltage=3.0,
        max_generations=3,
        max_species=40
    )
    
    print(f"\n结果摘要:")
    print(f"  - 生成了 {len(result['molecules'])} 个分子")
    print(f"  - 发现了 {len(result['reactions'])} 个反应")
    
    # 显示一些反应
    print(f"\n前5个反应:")
    for i, rxn in enumerate(result['reactions'][:5], 1):
        reactants = ' + '.join([r.smiles for r in rxn.reactants])
        products = ' + '.join([p.smiles for p in rxn.products])
        print(f"  {i}. {reactants} → {products}")
    
    return result


def example_3_high_temperature():
    """示例3: 高温反应"""
    print("\n" + "="*60)
    print("示例 3: 高温反应 (甲烷)")
    print("="*60)
    
    result = generate_reaction_network(
        smiles_list=['C'],  # 甲烷
        temperature=800.0,
        electrode_type='anode',
        voltage=0.0,
        max_generations=2,
        max_species=20
    )
    
    print(f"\n结果摘要:")
    print(f"  - 生成了 {len(result['molecules'])} 个分子")
    print(f"  - 发现了 {len(result['reactions'])} 个反应")
    
    return result


def example_4_custom_output():
    """示例4: 自定义输出目录"""
    print("\n" + "="*60)
    print("示例 4: 自定义输出目录")
    print("="*60)
    
    result = generate_reaction_network(
        smiles_list=['CC(=O)O'],  # 乙酸
        temperature=350.0,
        max_generations=2,
        output_dir='./my_custom_results'
    )
    
    print(f"\n输出文件:")
    if 'visualization_path' in result:
        print(f"  - 可视化: {result['visualization_path']}")
    if 'json_path' in result:
        print(f"  - JSON结果: {result['json_path']}")
    
    return result


def example_5_programmatic_access():
    """示例5: 程序化访问结果"""
    print("\n" + "="*60)
    print("示例 5: 程序化访问结果")
    print("="*60)
    
    result = generate_reaction_network(
        smiles_list=['C1CC1'],  # 环丙烷
        temperature=500.0,
        max_generations=2,
        max_species=15,
        visualize=False,  # 不生成可视化
        save_results=False  # 不保存到文件
    )
    
    # 访问网络对象
    network = result['network']
    
    # 分析分子
    print(f"\n分子列表:")
    for mol in result['molecules']:
        print(f"  - {mol.name}: {mol.smiles}")
    
    # 分析反应
    print(f"\n反应列表:")
    for i, rxn in enumerate(result['reactions'], 1):
        reactants = ', '.join([r.name for r in rxn.reactants])
        products = ', '.join([p.name for p in rxn.products])
        operator = rxn.operator_name if hasattr(rxn, 'operator_name') else 'unknown'
        print(f"  {i}. [{operator}] {reactants} → {products}")
    
    # 统计信息
    stats = result['statistics']
    print(f"\n统计信息:")
    for key, value in stats.items():
        print(f"  - {key}: {value}")
    
    return result


def main():
    """运行所有示例"""
    print("\n" + "="*60)
    print("RSNet Simple API - 使用示例")
    print("="*60)
    
    examples = [
        ("基本用法", example_1_basic_usage),
        ("有机分子", example_2_organic_molecules),
        ("高温反应", example_3_high_temperature),
        ("自定义输出", example_4_custom_output),
        ("程序化访问", example_5_programmatic_access)
    ]
    
    print("\n可用示例:")
    for i, (name, _) in enumerate(examples, 1):
        print(f"  {i}. {name}")
    
    print("\n选择要运行的示例 (1-5, 或按回车运行示例1):")
    choice = input("> ").strip()
    
    if not choice:
        choice = "1"
    
    try:
        idx = int(choice) - 1
        if 0 <= idx < len(examples):
            name, func = examples[idx]
            print(f"\n运行示例: {name}")
            result = func()
            print(f"\n✅ 示例完成！")
        else:
            print("无效的选择")
    except ValueError:
        print("无效的输入")
    except Exception as e:
        print(f"\n❌ 错误: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()
