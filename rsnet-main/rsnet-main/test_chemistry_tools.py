"""
化学工具函数测试

验证chemistry_tools.py中的各项功能
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rdkit import Chem
# 直接导入chemistry_tools模块，避免通过__init__.py
import rsnet.utils.chemistry_tools as chem_tools
from rsnet.utils.chemistry_tools import (
    calculate_huckel_energy, find_pi_system, estimate_homo_lumo_gap,
    calculate_pauling_electronegativity, calculate_bond_polarity,
    estimate_bond_dissociation_energy, get_base_bond_energy, calculate_radical_stabilization,
    calculate_ring_strain, get_empirical_ring_strain,
    detect_radicals, is_radical, find_conjugated_systems
)


def test_huckel_theory():
    """测试Hückel理论计算"""
    print("\n" + "="*60)
    print("测试 1: Hückel理论 - HOMO/LUMO计算")
    print("="*60)
    
    # 测试1: 乙烯（最简单的π体系）
    print("\n1.1 乙烯 (C=C)")
    ethene = Chem.MolFromSmiles('C=C')
    result = calculate_huckel_energy(ethene)
    print(f"  π原子数: {len(result['pi_atoms'])}")
    print(f"  π电子数: {result['num_pi_electrons']}")
    print(f"  HOMO能级: {result['homo_energy']}")
    print(f"  LUMO能级: {result['lumo_energy']}")
    if result['homo_energy'] and result['lumo_energy']:
        gap = result['lumo_energy'] - result['homo_energy']
        print(f"  能隙: {gap:.3f} β")
    
    # 测试2: 苯（芳香体系）
    print("\n1.2 苯 (芳香环)")
    benzene = Chem.MolFromSmiles('c1ccccc1')
    result = calculate_huckel_energy(benzene)
    print(f"  π原子数: {len(result['pi_atoms'])}")
    print(f"  π电子数: {result['num_pi_electrons']}")
    print(f"  HOMO能级: {result['homo_energy']}")
    print(f"  LUMO能级: {result['lumo_energy']}")
    if result['homo_energy'] and result['lumo_energy']:
        gap = result['lumo_energy'] - result['homo_energy']
        print(f"  能隙: {gap:.3f} β")
    
    # 测试3: 丁二烯（共轭体系）
    print("\n1.3 丁二烯 (C=C-C=C)")
    butadiene = Chem.MolFromSmiles('C=CC=C')
    result = calculate_huckel_energy(butadiene)
    print(f"  π原子数: {len(result['pi_atoms'])}")
    print(f"  π电子数: {result['num_pi_electrons']}")
    print(f"  HOMO能级: {result['homo_energy']}")
    print(f"  LUMO能级: {result['lumo_energy']}")
    if result['homo_energy'] and result['lumo_energy']:
        gap = result['lumo_energy'] - result['homo_energy']
        print(f"  能隙: {gap:.3f} β")


def test_bde_calculation():
    """测试键解离能计算"""
    print("\n" + "="*60)
    print("测试 2: 键解离能（BDE）计算")
    print("="*60)
    
    # 测试不同类型的C-H键
    molecules = [
        ('甲烷 CH4', 'C'),
        ('乙烷 CH3-CH3', 'CC'),
        ('丙烷 CH3-CH2-CH3', 'CCC'),
        ('异丁烷 (CH3)3CH', 'CC(C)C'),
        ('甲苯 PhCH3', 'Cc1ccccc1'),
        ('丙烯 CH2=CH-CH3', 'C=CC'),
    ]
    
    for name, smiles in molecules:
        mol = Chem.MolFromSmiles(smiles)
        print(f"\n{name}:")
        
        # 找到C-H键并计算BDE
        for bond in mol.GetBonds():
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            
            if (atom1.GetSymbol() == 'C' and atom2.GetSymbol() == 'H') or \
               (atom1.GetSymbol() == 'H' and atom2.GetSymbol() == 'C'):
                bde = estimate_bond_dissociation_energy(bond, mol)
                c_atom = atom1 if atom1.GetSymbol() == 'C' else atom2
                
                # 判断碳的类型
                num_c_neighbors = sum(1 for n in c_atom.GetNeighbors() if n.GetSymbol() == 'C')
                if num_c_neighbors == 0:
                    c_type = "伯碳"
                elif num_c_neighbors == 1:
                    c_type = "伯碳"
                elif num_c_neighbors == 2:
                    c_type = "仲碳"
                elif num_c_neighbors == 3:
                    c_type = "叔碳"
                else:
                    c_type = "未知"
                
                print(f"  C-H键 ({c_type}): BDE = {bde:.1f} kcal/mol")
                break  # 只显示第一个C-H键


def test_ring_strain():
    """测试环张力计算"""
    print("\n" + "="*60)
    print("测试 3: 环张力计算")
    print("="*60)
    
    rings = [
        ('环丙烷', 'C1CC1', 3),
        ('环丁烷', 'C1CCC1', 4),
        ('环戊烷', 'C1CCCC1', 5),
        ('环己烷', 'C1CCCCC1', 6),
    ]
    
    for name, smiles, size in rings:
        mol = Chem.MolFromSmiles(smiles)
        
        # 获取环原子
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() > 0:
            ring_atoms = list(ring_info.AtomRings()[0])
            
            calculated_strain = calculate_ring_strain(ring_atoms, mol)
            empirical_strain = get_empirical_ring_strain(size)
            
            print(f"\n{name} ({size}元环):")
            print(f"  计算张力: {calculated_strain:.1f} kcal/mol")
            print(f"  经验张力: {empirical_strain:.1f} kcal/mol")


def test_radical_detection():
    """测试自由基检测"""
    print("\n" + "="*60)
    print("测试 4: 自由基检测")
    print("="*60)
    
    # 创建一个自由基分子（甲基自由基）
    print("\n4.1 甲基自由基 [CH3]·")
    methyl_radical = Chem.MolFromSmiles('[CH3]')
    radicals = detect_radicals(methyl_radical)
    print(f"  自由基位点: {radicals}")
    print(f"  是自由基: {is_radical(methyl_radical)}")
    
    # 普通分子
    print("\n4.2 甲烷 CH4")
    methane = Chem.MolFromSmiles('C')
    radicals = detect_radicals(methane)
    print(f"  自由基位点: {radicals}")
    print(f"  是自由基: {is_radical(methane)}")


def test_conjugation():
    """测试共轭体系检测"""
    print("\n" + "="*60)
    print("测试 5: 共轭体系检测")
    print("="*60)
    
    molecules = [
        ('乙烯', 'C=C'),
        ('丁二烯', 'C=CC=C'),
        ('苯', 'c1ccccc1'),
        ('萘', 'c1ccc2ccccc2c1'),
    ]
    
    for name, smiles in molecules:
        mol = Chem.MolFromSmiles(smiles)
        systems = find_conjugated_systems(mol)
        
        print(f"\n{name}:")
        print(f"  共轭体系数: {len(systems)}")
        for i, system in enumerate(systems, 1):
            print(f"  体系{i}: {len(system)}个原子")


def main():
    """运行所有测试"""
    print("\n" + "="*60)
    print("化学工具函数测试套件")
    print("="*60)
    
    try:
        test_huckel_theory()
        test_bde_calculation()
        test_ring_strain()
        test_radical_detection()
        test_conjugation()
        
        print("\n" + "="*60)
        print("✅ 所有测试完成！")
        print("="*60)
        
    except Exception as e:
        print(f"\n❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()
