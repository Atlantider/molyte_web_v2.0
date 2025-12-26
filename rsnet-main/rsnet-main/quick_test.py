#!/usr/bin/env python3
"""
简单测试脚本 - 验证RSNet Simple API是否正常工作
"""

import sys
import os

# 添加当前目录到Python路径
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def test_import():
    """测试导入"""
    print("测试 1/4: 导入模块...")
    try:
        from rsnet_simple_api import generate_reaction_network, RSNetSimpleAPI
        print("  ✓ 导入成功")
        return True
    except Exception as e:
        print(f"  ✗ 导入失败: {e}")
        return False


def test_basic_generation():
    """测试基本网络生成"""
    print("\n测试 2/4: 基本网络生成...")
    try:
        from rsnet_simple_api import generate_reaction_network
        
        result = generate_reaction_network(
            smiles_list=['C'],  # 甲烷
            temperature=300.0,
            max_generations=1,
            max_species=5,
            visualize=False,
            save_results=False
        )
        
        assert 'network' in result
        assert 'molecules' in result
        assert 'reactions' in result
        assert len(result['molecules']) >= 1
        
        print(f"  ✓ 生成成功")
        print(f"    - 分子数: {len(result['molecules'])}")
        print(f"    - 反应数: {len(result['reactions'])}")
        return True
        
    except Exception as e:
        print(f"  ✗ 生成失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_multiple_molecules():
    """测试多分子输入"""
    print("\n测试 3/4: 多分子输入...")
    try:
        from rsnet_simple_api import generate_reaction_network
        
        result = generate_reaction_network(
            smiles_list=['CCO', 'C=C'],  # 乙醇 + 乙烯
            temperature=300.0,
            max_generations=1,
            max_species=10,
            visualize=False,
            save_results=False
        )
        
        assert len(result['molecules']) >= 2
        print(f"  ✓ 测试通过")
        print(f"    - 分子数: {len(result['molecules'])}")
        return True
        
    except Exception as e:
        print(f"  ✗ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_visualization():
    """测试可视化生成"""
    print("\n测试 4/4: 可视化生成...")
    try:
        from rsnet_simple_api import generate_reaction_network
        import os
        
        result = generate_reaction_network(
            smiles_list=['C'],
            temperature=300.0,
            max_generations=1,
            max_species=5,
            visualize=True,
            save_results=True,
            output_dir='./test_output'
        )
        
        assert 'visualization_path' in result
        assert os.path.exists(result['visualization_path'])
        assert 'json_path' in result
        assert os.path.exists(result['json_path'])
        
        print(f"  ✓ 测试通过")
        print(f"    - 可视化: {result['visualization_path']}")
        print(f"    - JSON: {result['json_path']}")
        return True
        
    except Exception as e:
        print(f"  ✗ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """运行所有测试"""
    print("="*60)
    print("RSNet Simple API - 快速测试")
    print("="*60)
    
    tests = [
        test_import,
        test_basic_generation,
        test_multiple_molecules,
        test_visualization
    ]
    
    results = []
    for test in tests:
        results.append(test())
    
    print("\n" + "="*60)
    print("测试总结")
    print("="*60)
    passed = sum(results)
    total = len(results)
    print(f"通过: {passed}/{total}")
    
    if passed == total:
        print("\n✅ 所有测试通过！")
        return 0
    else:
        print(f"\n❌ {total - passed} 个测试失败")
        return 1


if __name__ == '__main__':
    exit(main())
