"""
Oxidation Reaction Operators for RSNet CEI Chemistry
氧化反应算符 - 用于阴极CEI界面化学反应

包含的反应类型:
1. H抽取 (H-Abstraction)
2. 单电子氧化 (One-electron oxidation)
3. 过氧化物/超氧化物反应
4. 单线态氧攻击
5. 过渡金属配位
6. 碳酸酯开环氧化

作者: RSNet Team
版本: 2.0
"""

from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field
from enum import Enum


class OxidationMechanism(Enum):
    """氧化反应机制"""
    H_ABSTRACTION = "h_abstraction"           # H抽取
    ELECTRON_TRANSFER = "electron_transfer"    # 电子转移
    OXYGEN_INSERTION = "oxygen_insertion"      # 氧插入
    PEROXIDE_ATTACK = "peroxide_attack"        # 过氧化物攻击
    SUPEROXIDE_ATTACK = "superoxide_attack"    # 超氧化物攻击
    SINGLET_O2_ATTACK = "singlet_o2_attack"   # 单线态氧攻击
    RADICAL_COUPLING = "radical_coupling"      # 自由基偶联
    METAL_CATALYZED = "metal_catalyzed"        # 金属催化


@dataclass
class OxidationOperator:
    """氧化反应算符"""
    name: str                                  # 反应名称
    mechanism: OxidationMechanism              # 反应机制
    reactant_patterns: List[str]              # 反应物模式
    oxidant_species: List[str]                # 氧化剂
    product_patterns: List[str]               # 产物模式
    activation_energy_range: Tuple[float, float]  # 活化能范围 (kcal/mol)
    reaction_energy_estimate: float           # 反应能估计 (kcal/mol)
    conditions: List[str]                     # 反应条件
    description: str                          # 描述
    bep_alpha: float = 0.4                    # BEP α参数
    bep_beta: float = 15.0                    # BEP β参数


# =============================================================================
# 预定义的氧化反应算符
# =============================================================================

OXIDATION_OPERATORS: Dict[str, OxidationOperator] = {
    # =========================================================================
    # H抽取反应 (R-H + O·/O⁻ → R· + OH/OH⁻)
    # =========================================================================
    'h_abstraction_by_O-': OxidationOperator(
        name='H-abstraction by O⁻',
        mechanism=OxidationMechanism.H_ABSTRACTION,
        reactant_patterns=['[CH]', 'C[H]', '[CH2]', '[CH3]'],  # C-H键
        oxidant_species=['[O-]'],
        product_patterns=['[C]', '[C·]', '[OH-]'],  # 碳自由基 + 氢氧根
        activation_energy_range=(8.0, 18.0),
        reaction_energy_estimate=-5.0,  # 通常放能
        conditions=['cathode', 'high_voltage'],
        description='氧自由基阴离子从碳上抽取氢原子，生成碳自由基',
        bep_alpha=0.35,
        bep_beta=12.0
    ),
    
    'h_abstraction_by_OH': OxidationOperator(
        name='H-abstraction by ·OH',
        mechanism=OxidationMechanism.H_ABSTRACTION,
        reactant_patterns=['[CH]', 'C[H]'],
        oxidant_species=['[OH]'],
        product_patterns=['[C]', 'O'],  # 碳自由基 + H2O
        activation_energy_range=(3.0, 12.0),
        reaction_energy_estimate=-15.0,  # 强放能
        conditions=['cathode', 'trace_water'],
        description='羟基自由基抽取氢原子，反应非常快',
        bep_alpha=0.30,
        bep_beta=8.0
    ),
    
    'h_abstraction_by_superoxide': OxidationOperator(
        name='H-abstraction by O₂⁻',
        mechanism=OxidationMechanism.SUPEROXIDE_ATTACK,
        reactant_patterns=['[CH]', 'C[H]', '[CH2]'],
        oxidant_species=['[O][O-]'],  # 超氧根
        product_patterns=['[C]', '[O][O]'],  # 碳自由基 + HO2·
        activation_energy_range=(10.0, 22.0),
        reaction_energy_estimate=-3.0,
        conditions=['cathode', 'oxygen_release'],
        description='超氧根从碳上抽取氢，生成HO₂·和碳自由基',
        bep_alpha=0.40,
        bep_beta=15.0
    ),
    
    # =========================================================================
    # 单电子氧化 (R → R⁺· + e⁻)
    # =========================================================================
    'one_electron_oxidation': OxidationOperator(
        name='One-electron Oxidation',
        mechanism=OxidationMechanism.ELECTRON_TRANSFER,
        reactant_patterns=['C', 'O', 'N'],  # 具有孤对电子的原子
        oxidant_species=['h+', 'electrode'],  # 空穴或电极
        product_patterns=['[C+]', '[O+]', '[N+]'],  # 阳离子自由基
        activation_energy_range=(5.0, 15.0),
        reaction_energy_estimate=10.0,  # 取决于HOMO能级
        conditions=['cathode', 'high_voltage'],
        description='分子失去一个电子形成阳离子自由基',
        bep_alpha=0.35,
        bep_beta=10.0
    ),
    
    'carbonate_oxidation': OxidationOperator(
        name='Carbonate Oxidation',
        mechanism=OxidationMechanism.ELECTRON_TRANSFER,
        reactant_patterns=['C1COC(=O)O1', 'COC(=O)OC'],  # EC, DMC
        oxidant_species=['h+'],
        product_patterns=['[C1COC(=O)O1+]', '[COC(=O)OC+]'],  # 阳离子自由基
        activation_energy_range=(8.0, 18.0),
        reaction_energy_estimate=15.0,
        conditions=['cathode', 'high_voltage'],
        description='碳酸酯单电子氧化，是CEI形成的起始步骤',
        bep_alpha=0.38,
        bep_beta=12.0
    ),
    
    # =========================================================================
    # 过氧化物反应 (R + O₂²⁻ → R-OO⁻)
    # =========================================================================
    'peroxide_nucleophilic_attack': OxidationOperator(
        name='Peroxide Nucleophilic Attack',
        mechanism=OxidationMechanism.PEROXIDE_ATTACK,
        reactant_patterns=['C(=O)', '[C+]'],  # 羰基碳或碳正离子
        oxidant_species=['[O-][O-]'],  # 过氧根
        product_patterns=['C(O[O-])'],  # 过氧化物加成产物
        activation_energy_range=(12.0, 25.0),
        reaction_energy_estimate=-8.0,
        conditions=['cathode', 'li2o2_formation'],
        description='过氧根作为亲核试剂攻击缺电子碳中心',
        bep_alpha=0.42,
        bep_beta=18.0
    ),
    
    'peroxide_decomposition': OxidationOperator(
        name='Peroxide Decomposition',
        mechanism=OxidationMechanism.PEROXIDE_ATTACK,
        reactant_patterns=['[O-][O-]', '[Li]OO[Li]'],
        oxidant_species=[],  # 自分解
        product_patterns=['[O-]', '[O]', 'O=O'],  # 氧离子 + 氧气
        activation_energy_range=(18.0, 30.0),
        reaction_energy_estimate=5.0,
        conditions=['cathode', 'high_voltage'],
        description='过氧化物分解释放活性氧',
        bep_alpha=0.45,
        bep_beta=22.0
    ),
    
    # =========================================================================
    # 超氧化物反应 (R + O₂⁻ → R-OO·)
    # =========================================================================
    'superoxide_addition': OxidationOperator(
        name='Superoxide Addition to Double Bond',
        mechanism=OxidationMechanism.SUPEROXIDE_ATTACK,
        reactant_patterns=['C=C', 'C=O'],  # 不饱和键
        oxidant_species=['[O][O-]'],
        product_patterns=['C(O[O])C', 'C(O[O])O'],  # 过氧化物
        activation_energy_range=(8.0, 18.0),
        reaction_energy_estimate=-12.0,
        conditions=['cathode', 'oxygen_release'],
        description='超氧根加成到双键形成过氧化物',
        bep_alpha=0.35,
        bep_beta=12.0
    ),
    
    'superoxide_sn2': OxidationOperator(
        name='Superoxide SN2 Attack',
        mechanism=OxidationMechanism.SUPEROXIDE_ATTACK,
        reactant_patterns=['COC', 'COC(=O)'],  # 醚或酯
        oxidant_species=['[O][O-]'],
        product_patterns=['CO', 'OO'],  # 断裂产物
        activation_energy_range=(15.0, 28.0),
        reaction_energy_estimate=-5.0,
        conditions=['cathode'],
        description='超氧根对C-O键的亲核取代',
        bep_alpha=0.40,
        bep_beta=20.0
    ),
    
    # =========================================================================
    # 单线态氧反应 (R-H + ¹O₂ → R-OOH)
    # =========================================================================
    'singlet_oxygen_ene_reaction': OxidationOperator(
        name='Singlet Oxygen Ene Reaction',
        mechanism=OxidationMechanism.SINGLET_O2_ATTACK,
        reactant_patterns=['C=C', 'CC=C'],  # 烯烃
        oxidant_species=['O=O'],  # 单线态O2 (SMILES相同)
        product_patterns=['C(O)C=C', 'CC(OO)'],  # 烯丙基过氧化物
        activation_energy_range=(3.0, 10.0),
        reaction_energy_estimate=-25.0,  # 强放能
        conditions=['cathode', 'very_high_voltage'],
        description='单线态氧的ene反应，非常快速',
        bep_alpha=0.25,
        bep_beta=6.0
    ),
    
    'singlet_oxygen_cycloaddition': OxidationOperator(
        name='Singlet Oxygen [2+2] Cycloaddition',
        mechanism=OxidationMechanism.SINGLET_O2_ATTACK,
        reactant_patterns=['C=C'],  # 烯烃
        oxidant_species=['O=O'],
        product_patterns=['C1OOC1'],  # 二氧杂环丁烷
        activation_energy_range=(5.0, 15.0),
        reaction_energy_estimate=-20.0,
        conditions=['cathode', 'very_high_voltage'],
        description='单线态氧与双键的[2+2]环加成',
        bep_alpha=0.30,
        bep_beta=8.0
    ),
    
    # =========================================================================
    # 过渡金属催化反应
    # =========================================================================
    'metal_catalyzed_oxidation': OxidationOperator(
        name='Transition Metal Catalyzed Oxidation',
        mechanism=OxidationMechanism.METAL_CATALYZED,
        reactant_patterns=['C[OH]', 'CC(O)'],  # 醇
        oxidant_species=['[Co+3]', '[Ni+3]', '[Mn+3]'],
        product_patterns=['C=O', 'CC(=O)'],  # 醛/酮
        activation_energy_range=(8.0, 18.0),
        reaction_energy_estimate=-10.0,
        conditions=['cathode', 'metal_dissolution'],
        description='溶出过渡金属催化的醇氧化',
        bep_alpha=0.38,
        bep_beta=12.0
    ),
    
    'metal_coordination_complex': OxidationOperator(
        name='Metal-Solvent Complex Formation',
        mechanism=OxidationMechanism.METAL_CATALYZED,
        reactant_patterns=['C1COC(=O)O1', 'COC(=O)OC'],  # 碳酸酯
        oxidant_species=['[Co+2]', '[Ni+2]', '[Mn+2]'],
        product_patterns=['[M]OC'],  # 金属配合物
        activation_energy_range=(2.0, 8.0),
        reaction_energy_estimate=-15.0,  # 配位放能
        conditions=['cathode', 'metal_dissolution'],
        description='溶出金属与溶剂分子配位',
        bep_alpha=0.20,
        bep_beta=5.0
    ),
    
    # =========================================================================
    # 碳酸酯特定反应
    # =========================================================================
    'ec_ring_opening_oxidative': OxidationOperator(
        name='EC Oxidative Ring Opening',
        mechanism=OxidationMechanism.OXYGEN_INSERTION,
        reactant_patterns=['C1COC(=O)O1'],  # EC
        oxidant_species=['[O-]', '[O][O-]'],
        product_patterns=['OCCO', 'OCC(=O)O', 'C(=O)O'],  # 开环产物
        activation_energy_range=(12.0, 22.0),
        reaction_energy_estimate=-8.0,
        conditions=['cathode'],
        description='EC在氧化条件下开环',
        bep_alpha=0.38,
        bep_beta=16.0
    ),
    
    'dmc_ch_oxidation': OxidationOperator(
        name='DMC C-H Oxidation',
        mechanism=OxidationMechanism.H_ABSTRACTION,
        reactant_patterns=['COC(=O)OC'],  # DMC
        oxidant_species=['[O-]', '[OH]'],
        product_patterns=['C(=O)OC(=O)OC', 'OC(=O)OC'],  # 氧化产物
        activation_energy_range=(10.0, 20.0),
        reaction_energy_estimate=-5.0,
        conditions=['cathode'],
        description='DMC甲基上的C-H氧化',
        bep_alpha=0.40,
        bep_beta=15.0
    ),
}


class OxidationOperatorManager:
    """氧化反应算符管理器"""
    
    def __init__(self, 
                 include_peroxide: bool = True,
                 include_superoxide: bool = True,
                 include_singlet_oxygen: bool = True,
                 include_metal_catalysis: bool = True):
        """
        初始化管理器
        
        Args:
            include_peroxide: 包含过氧化物反应
            include_superoxide: 包含超氧化物反应
            include_singlet_oxygen: 包含单线态氧反应
            include_metal_catalysis: 包含金属催化反应
        """
        self.include_peroxide = include_peroxide
        self.include_superoxide = include_superoxide
        self.include_singlet_oxygen = include_singlet_oxygen
        self.include_metal_catalysis = include_metal_catalysis
        
        self.operators = self._filter_operators()
    
    def _filter_operators(self) -> Dict[str, OxidationOperator]:
        """根据设置筛选算符"""
        filtered = {}
        
        for key, op in OXIDATION_OPERATORS.items():
            include = True
            
            if op.mechanism == OxidationMechanism.PEROXIDE_ATTACK and not self.include_peroxide:
                include = False
            elif op.mechanism == OxidationMechanism.SUPEROXIDE_ATTACK and not self.include_superoxide:
                include = False
            elif op.mechanism == OxidationMechanism.SINGLET_O2_ATTACK and not self.include_singlet_oxygen:
                include = False
            elif op.mechanism == OxidationMechanism.METAL_CATALYZED and not self.include_metal_catalysis:
                include = False
            
            if include:
                filtered[key] = op
        
        return filtered
    
    def get_applicable_operators(self, 
                                  smiles: str, 
                                  available_oxidants: List[str]) -> List[OxidationOperator]:
        """
        获取对给定分子适用的氧化算符
        
        Args:
            smiles: 分子SMILES
            available_oxidants: 可用的氧化剂列表
            
        Returns:
            适用的算符列表
        """
        applicable = []
        
        for op in self.operators.values():
            # 检查反应物模式匹配
            reactant_match = any(pattern in smiles for pattern in op.reactant_patterns)
            
            # 检查氧化剂是否可用
            if not op.oxidant_species:
                oxidant_match = True  # 自分解反应
            else:
                oxidant_match = any(ox in available_oxidants for ox in op.oxidant_species)
            
            if reactant_match and oxidant_match:
                applicable.append(op)
        
        return applicable
    
    def estimate_activation_energy(self, 
                                    operator: OxidationOperator,
                                    reaction_energy: float) -> float:
        """
        使用BEP关系估算活化能
        
        Args:
            operator: 氧化算符
            reaction_energy: 反应能 (kcal/mol)
            
        Returns:
            估算的活化能 (kcal/mol)
        """
        alpha = operator.bep_alpha
        beta = operator.bep_beta
        
        if reaction_energy > 0:  # 吸能
            Ea = alpha * reaction_energy + beta
        else:  # 放能 (Hammond假说)
            Ea = beta - (1 - alpha) * abs(reaction_energy)
            Ea = max(2.0, Ea)
        
        # 确保在算符定义的范围内
        Ea_min, Ea_max = operator.activation_energy_range
        Ea = max(Ea_min, min(Ea_max, Ea))
        
        return Ea
    
    def list_operators_by_mechanism(self, mechanism: OxidationMechanism) -> List[OxidationOperator]:
        """按机制分类获取算符"""
        return [op for op in self.operators.values() if op.mechanism == mechanism]
    
    def get_all_oxidant_species(self) -> List[str]:
        """获取所有算符需要的氧化剂物种"""
        species = set()
        for op in self.operators.values():
            species.update(op.oxidant_species)
        return list(species)


# 便捷函数
def get_oxidation_operators_for_cathode(voltage: float = 4.2) -> Dict[str, OxidationOperator]:
    """获取阴极条件下的氧化算符"""
    manager = OxidationOperatorManager()
    return manager.operators


def estimate_oxidation_barrier(reaction_name: str, reaction_energy: float) -> Optional[float]:
    """估算特定氧化反应的活化能"""
    if reaction_name in OXIDATION_OPERATORS:
        op = OXIDATION_OPERATORS[reaction_name]
        manager = OxidationOperatorManager()
        return manager.estimate_activation_energy(op, reaction_energy)
    return None
