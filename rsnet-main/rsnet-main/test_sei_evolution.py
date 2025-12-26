"""
智能网络演化系统测试 - 锂离子电池SEI膜形成模拟
"""

import sys
import os
import logging

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# 添加项目路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.network.evolution import NetworkEvolver

# Mock xTB Calculator
from unittest.mock import MagicMock
import rsnet.compute.xtb_calculator

class MockXTBCalculator:
    def __init__(self, *args, **kwargs):
        pass
    def calculate_energy(self, mol, *args, **kwargs):
        return -100.0  # Dummy energy
    def optimize_geometry(self, mol, *args, **kwargs):
        return mol  # Return same mol
    def _check_xtb_availability(self):
        pass  # Bypass check

# Patch the real class in BOTH locations
import rsnet.compute.reaction_screener
rsnet.compute.xtb_calculator.XTBCalculator = MockXTBCalculator
rsnet.compute.reaction_screener.XTBCalculator = MockXTBCalculator


def test_sei_evolution():
    print("\n" + "="*70)
    print("SEI膜形成演化模拟 (Intelligent Network Evolution)")
    print("="*70)
    
    # 1. 初始反应物
    # 碳酸乙烯酯 (EC)
    ec = Molecule.from_smiles('C1COC(=O)O1', name='EC')
    # 锂离子 (Li+)
    li_ion = Molecule.from_smiles('[Li+]', name='Li+')
    
    # 2. 反应环境 - 电池阳极界面
    env = Environment(
        temperature=298.15,
        pressure=1.0,
        solvent='EC',
        electrode_type='anode',
        voltage=0.1,  # 0.1V vs Li/Li+ (Strong reduction -> Radicals)
        li_activity=1.0,
        interface_type='SEI'
    )
    
    # 驱动力将根据环境参数自动计算
    # env.get_active_drives() 会返回: reduction=True, electrochemical=True, etc.
    
    print(f"初始物种: {ec.name} ({ec.smiles}), {li_ion.name} ({li_ion.smiles})")
    print(f"环境: {env.electrode_type} (V={env.voltage}V)")
    
    # 3. 初始化演化器
    evolver = NetworkEvolver(environment=env)
    
    # 4. 运行演化 (3代)
    # Gen 0: EC, Li+ -> EC- (Radical Anion)
    # Gen 1: EC- -> Ring Open Radical
    # Gen 2: Radical + Radical -> LEDC / Polymer
    final_species = evolver.evolve([ec, li_ion], max_generations=3)
    
    print("\n" + "="*70)
    print("演化结果统计")
    print("="*70)
    
    radicals = [s for s in final_species if s.is_radical]
    neutrals = [s for s in final_species if not s.is_radical]
    
    print(f"总物种数: {len(final_species)}")
    print(f"自由基: {len(radicals)}")
    print(f"中性分子: {len(neutrals)}")
    
    print("\n关键产物检测:")
    
    # 检查是否有还原产物 (含自由基电子或负电荷)
    reduced_products = [s for s in final_species if 'red' in s.molecule.name or s.is_radical]
    print(f"- 还原产物/自由基: {len(reduced_products)} 个")
    for s in reduced_products[:3]: # 列出前几个
        print(f"  * {s.molecule.name}: {s.molecule.smiles} (Gen {s.generation})")

    # 检查是否有开环产物 (名字可能包含 opened 或 ring_opening 操作后的特征)
    # 由于名字是自动生成的，我们检查大小或名字
    opened_products = [s for s in final_species if 'opened' in s.molecule.name or ('EC' in s.molecule.name and '[' not in s.molecule.smiles and '1' not in s.molecule.smiles)] # 简单的开环检查：没有环标记1
    print(f"- 开环产物: {len(opened_products)} 个")
    for s in opened_products[:3]:
        print(f"  * {s.molecule.name}: {s.molecule.smiles}")
        
    # 检查是否有聚合/偶联产物
    coupled_products = [s for s in final_species if 'coupled' in s.molecule.name or 'polymer' in s.molecule.name]
    print(f"- 聚合/偶联产物: {len(coupled_products)} 个 (SEI主要成分)")
    for s in coupled_products[:3]:
        print(f"  * {s.molecule.name}: {s.molecule.smiles}")
        
    print("\n演化路径验证:")
    if reduced_products:
        print("✅ 检测到电化学还原 (EC -> EC-•)")
    else:
        print("❌ 未检测到还原反应 (检查 ElectronTransfer 算符)")
        
    if opened_products:
        print("✅ 检测到环开环反应 (EC-• -> O-CH2-CH2-O-C•=O)")
    else:
        print("⚠️ 未明显检测到开环 (可能在下一步或条件未满足)")
        
    if coupled_products:
        print("✅ 检测到自由基偶联 (生成 LEDC 前体)")
    else:
        print("⚠️ 未检测到偶联 (可能需要更多代数或自由基浓度)")

if __name__ == "__main__":
    test_sei_evolution()
