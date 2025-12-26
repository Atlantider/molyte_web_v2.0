"""
智能算符激活系统 - 简化演示

展示激活规则配置和评分机制
"""

# 演示激活规则配置
from rsnet.operators.activation_rules import (
    OPERATOR_ACTIVATION_RULES,
    DRIVE_WEIGHTS,
    get_activation_rule,
    get_drive_weight
)

def demo_activation_rules():
    """演示激活规则配置"""
    print("\n" + "="*70)
    print("智能算符激活系统 - 配置演示")
    print("="*70)
    
    print("\n1. 电子转移算符 (ElectronTransferOperator)")
    print("-" * 70)
    rule = get_activation_rule('electron_transfer')
    print(f"化学原理: {rule['chemical_principle']}")
    print(f"必需驱动力: {rule['required_drives']}")
    print(f"增强驱动力: {rule['enhancing_drives']}")
    print(f"分子检查: {list(rule['molecular_checks'].keys())}")
    print(f"权重: {rule['weight']}")
    
    print("\n2. 自由基反应算符 (RadicalReactionOperator)")
    print("-" * 70)
    rule = get_activation_rule('radical_reaction')
    print(f"化学原理: {rule['chemical_principle']}")
    print(f"必需驱动力: {rule['required_drives']}")
    print(f"增强驱动力: {rule['enhancing_drives']}")
    print(f"分子检查: {list(rule['molecular_checks'].keys())}")
    
    print("\n3. 开环反应算符 (RingOpeningOperator)")
    print("-" * 70)
    rule = get_activation_rule('ring_opening')
    print(f"化学原理: {rule['chemical_principle']}")
    print(f"必需驱动力: {rule['required_drives']}")
    print(f"增强驱动力: {rule['enhancing_drives']}")
    print(f"必需特征: {rule['required_features']}")
    
    print("\n4. 驱动力权重配置")
    print("-" * 70)
    important_drives = [
        'electrochemical', 'oxidation', 'reduction',
        'ring_strain', 'high_temperature', 'sei_formation'
    ]
    for drive in important_drives:
        weight = get_drive_weight(drive)
        print(f"  {drive:25s}: {weight:.2f}")
    
    print("\n5. 所有算符概览")
    print("-" * 70)
    for op_name, rule in OPERATOR_ACTIVATION_RULES.items():
        print(f"\n  {op_name:20s} - {rule['description']}")
        print(f"  {'':20s}   权重: {rule['weight']:.2f}")


def demo_scoring_mechanism():
    """演示评分机制"""
    print("\n" + "="*70)
    print("评分机制说明")
    print("="*70)
    
    print("\n总分计算公式:")
    print("  score = 环境匹配度 × 0.4")
    print("        + 分子匹配度 × 0.3")
    print("        + 化学可行性 × 0.2")
    print("        + 算符权重 × 0.1")
    
    print("\n环境匹配度计算:")
    print("  1. 检查必需驱动力是否激活")
    print("  2. 计算驱动力强度加权和")
    print("  3. 必需驱动力权重 × 1.5")
    print("  4. 增强驱动力正常权重")
    
    print("\n分子匹配度计算:")
    print("  1. 检测分子结构特征")
    print("  2. 评估每个特征的匹配度")
    print("  3. 计算加权平均分")
    print("  4. 多分子取最高分")
    
    print("\n化学可行性:")
    print("  - 调用operator.can_apply()")
    print("  - 基于具体算符的化学规则")
    print("  - 返回0.0或1.0")


def demo_usage_example():
    """演示使用示例"""
    print("\n" + "="*70)
    print("使用示例")
    print("="*70)
    
    print("\n代码示例:")
    print("""
from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.operators.registry import OPERATOR_REGISTRY

# 创建分子
mol = Molecule.from_smiles('C1CC1', name='cyclopropane')

# 创建环境
env = Environment(
    temperature=450.0,
    pressure=1.0,
    solvent='gas'
)

# 获取智能激活的算符
operators = OPERATOR_REGISTRY.get_active_operators_smart(
    molecules=[mol],
    environment=env,
    min_score=0.1
)

# 查看结果
for operator, score in operators:
    print(f"{operator.name}: {score:.3f}")
""")
    
    print("\n预期输出:")
    print("  RingOpening: 0.850  # 高张力环优先")
    print("  Decomposition: 0.620  # 高温分解")
    print("  RadicalReaction: 0.450  # 自由基环境")


if __name__ == '__main__':
    demo_activation_rules()
    demo_scoring_mechanism()
    demo_usage_example()
    
    print("\n" + "="*70)
    print("✅ 演示完成！")
    print("="*70)
    print("\n系统特点:")
    print("  ✓ 基于真实化学和电化学原理")
    print("  ✓ 智能评分和优先级排序")
    print("  ✓ 环境-分子-化学三层匹配")
    print("  ✓ 可配置的激活规则")
    print("  ✓ 支持13个反应算符")
