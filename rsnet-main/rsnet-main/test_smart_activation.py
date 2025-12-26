"""
智能算符激活系统测试

测试基于化学和电化学原理的算符激活系统
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# 直接导入核心模块避免循环依赖
from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
# 延迟导入registry
import rsnet.operators.registry as registry_module



def test_high_temperature_decomposition():
    """测试高温分解场景"""
    print("\n" + "="*70)
    print("场景1：高温分解反应")
    print("="*70)
    
    # 创建过氧化物分子（弱O-O键）
    peroxide = Molecule.from_smiles('CCOOCC', name='diethyl_peroxide')
    
    # 高温环境
    env = Environment(
        temperature=600.0,  # 高温
        pressure=1.0,
        solvent='gas'
    )
    
    print(f"\n分子: {peroxide.name}")
    print(f"SMILES: {peroxide.smiles}")
    print(f"环境: T={env.temperature}K, P={env.pressure}atm, 溶剂={env.solvent}")
    
    # 获取激活的算符
    operators = registry_module.get_active_operators_smart([peroxide], env, min_score=0.1)

    
    print(f"\n激活的算符 (共{len(operators)}个):")
    for i, (op, score) in enumerate(operators, 1):
        print(f"  {i}. {op.name:20s} - 分数: {score:.3f}")


def test_electrochemical_reduction():
    """测试电化学还原场景"""
    print("\n" + "="*70)
    print("场景2：电化学还原反应（电池阳极）")
    print("="*70)
    
    # 碳酸乙烯酯（EC）
    ec = Molecule.from_smiles('C1COC(=O)O1', name='ethylene_carbonate')
    
    # 阳极还原环境
    env = Environment(
        temperature=298.15,
        electrode_type='anode',
        voltage=0.5,  # vs Li/Li+，低电压
        li_activity=1.0,
        interface_type='SEI'
    )
    
    print(f"\n分子: {ec.name}")
    print(f"SMILES: {ec.smiles}")
    print(f"环境: T={env.temperature}K, 电极={env.electrode_type}, V={env.voltage}V")
    print(f"      界面={env.interface_type}")
    
    # 获取激活的算符
    operators = registry_module.get_active_operators_smart([ec], env, min_score=0.1)
    
    print(f"\n激活的算符 (共{len(operators)}个):")
    for i, (op, score) in enumerate(operators, 1):
        print(f"  {i}. {op.name:20s} - 分数: {score:.3f}")


def test_ring_strain_release():
    """测试环张力释放场景"""
    print("\n" + "="*70)
    print("场景3：环张力释放反应")
    print("="*70)
    
    # 环丙烷（高张力）
    cyclopropane = Molecule.from_smiles('C1CC1', name='cyclopropane')
    
    # 中等温度环境
    env = Environment(
        temperature=450.0,
        pressure=1.0,
        solvent='gas'
    )
    
    print(f"\n分子: {cyclopropane.name}")
    print(f"SMILES: {cyclopropane.smiles}")
    print(f"环境: T={env.temperature}K, P={env.pressure}atm")
    
    # 获取激活的算符
    operators = registry_module.get_active_operators_smart([cyclopropane], env, min_score=0.1)
    
    print(f"\n激活的算符 (共{len(operators)}个):")
    for i, (op, score) in enumerate(operators, 1):
        print(f"  {i}. {op.name:20s} - 分数: {score:.3f}")


def test_radical_environment():
    """测试自由基环境场景"""
    print("\n" + "="*70)
    print("场景4：自由基反应环境")
    print("="*70)
    
    # 乙烯（可聚合）
    ethylene = Molecule.from_smiles('C=C', name='ethylene')
    
    # 自由基环境（高温 + SEI界面）
    env = Environment(
        temperature=500.0,
        interface_type='SEI',
        electrode_type='anode',
        voltage=0.8
    )
    
    print(f"\n分子: {ethylene.name}")
    print(f"SMILES: {ethylene.smiles}")
    print(f"环境: T={env.temperature}K, 界面={env.interface_type}")
    
    # 获取激活的算符
    operators = registry_module.get_active_operators_smart([ethylene], env, min_score=0.1)
    
    print(f"\n激活的算符 (共{len(operators)}个):")
    for i, (op, score) in enumerate(operators, 1):
        print(f"  {i}. {op.name:20s} - 分数: {score:.3f}")


def test_oxidation_reaction():
    """测试氧化反应场景"""
    print("\n" + "="*70)
    print("场景5：电化学氧化反应（电池阴极）")
    print("="*70)
    
    # 甲苯（有苄位）
    toluene = Molecule.from_smiles('Cc1ccccc1', name='toluene')
    
    # 阴极氧化环境
    env = Environment(
        temperature=298.15,
        electrode_type='cathode',
        voltage=4.2,  # vs Li/Li+，高电压
        interface_type='CEI'
    )
    
    print(f"\n分子: {toluene.name}")
    print(f"SMILES: {toluene.smiles}")
    print(f"环境: T={env.temperature}K, 电极={env.electrode_type}, V={env.voltage}V")
    
    # 获取激活的算符
    operators = registry_module.get_active_operators_smart([toluene], env, min_score=0.1)

    
    print(f"\n激活的算符 (共{len(operators)}个):")
    for i, (op, score) in enumerate(operators, 1):
        print(f"  {i}. {op.name:20s} - 分数: {score:.3f}")


def main():
    """运行所有测试"""
    print("\n" + "="*70)
    print("智能算符激活系统测试")
    print("基于化学和电化学原理")
    print("="*70)
    
    try:
        test_high_temperature_decomposition()
        test_electrochemical_reduction()
        test_ring_strain_release()
        test_radical_environment()
        test_oxidation_reaction()
        
        print("\n" + "="*70)
        print("✅ 所有测试完成！")
        print("="*70)
        print("\n说明：")
        print("- 分数越高，算符越适合当前条件")
        print("- 分数基于环境驱动力、分子特征和化学可行性")
        print("- 系统自动选择最合适的反应路径")
        
    except Exception as e:
        print(f"\n❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()
