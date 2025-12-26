"""
Improved Activation Energy Estimator for RSNet
基于文献的活化能估算器 - 使用反应类型特定的Bell-Evans-Polanyi参数

参考文献:
- Bell, R. P. (1936). 活化能与反应热的关系
- Evans, M. G. & Polanyi, M. (1938). Bell-Evans-Polanyi原理  
- Hammond, G. S. (1955). 过渡态假说
- Marcus, R. A. (1964). 电子转移理论
- Butler, J. A. V. (1924). 电化学动力学

作者: RSNet Team
版本: 2.0 (改进版)
"""

from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
from enum import Enum
import math


class ReactionType(Enum):
    """反应类型枚举"""
    ELECTRON_TRANSFER = "electron_transfer"      # 电子转移
    H_TRANSFER = "h_transfer"                    # 氢/质子转移
    RING_OPENING = "ring_opening"                # 环开裂
    BOND_CLEAVAGE = "bond_cleavage"              # 键断裂
    LI_COORDINATION = "li_coordination"          # 锂配位
    LI_INSERTION = "li_insertion"                # 锂嵌入
    OXIDATION = "oxidation"                      # 氧化
    REDUCTION = "reduction"                      # 还原
    RADICAL_FORMATION = "radical_formation"      # 自由基形成
    SUBSTITUTION = "substitution"                # 取代
    ELIMINATION = "elimination"                  # 消除
    POLYMERIZATION = "polymerization"            # 聚合
    UNKNOWN = "unknown"


@dataclass
class BEPParameters:
    """
    Bell-Evans-Polanyi参数
    
    Ea = α * ΔE + β (吸能反应)
    Ea = β - (1-α) * |ΔE| (放能反应, Hammond假说)
    
    其中:
    - α: Brønsted系数 (0-1), 描述过渡态与产物的相似程度
    - β: 本征势垒 (kcal/mol), 零反应热时的活化能
    """
    alpha: float  # Brønsted系数
    beta: float   # 本征势垒 (kcal/mol)
    source: str   # 参数来源/参考文献


class ImprovedActivationEnergyEstimator:
    """
    基于反应类型的活化能估算器
    
    特点:
    1. 使用文献来源的BEP参数
    2. 考虑反应类型差异
    3. 包含电化学校正
    4. 考虑环应变效应
    """
    
    # 反应类型特定的BEP参数 (α, β, 来源)
    BEP_PARAMETERS: Dict[ReactionType, BEPParameters] = {
        # 电子转移反应 - Marcus理论
        ReactionType.ELECTRON_TRANSFER: BEPParameters(
            alpha=0.35, beta=5.0, 
            source="Marcus, R.A. (1964) J. Chem. Phys."
        ),
        
        # 质子/氢转移 - 典型有机反应
        ReactionType.H_TRANSFER: BEPParameters(
            alpha=0.45, beta=12.0,
            source="Bell, R.P. (1936) Proc. R. Soc."
        ),
        
        # 氧化反应
        ReactionType.OXIDATION: BEPParameters(
            alpha=0.40, beta=15.0,
            source="Typical organic oxidation"
        ),
        
        # 还原反应
        ReactionType.REDUCTION: BEPParameters(
            alpha=0.38, beta=12.0,
            source="Electrochemical reduction"
        ),
        
        # 自由基形成
        ReactionType.RADICAL_FORMATION: BEPParameters(
            alpha=0.30, beta=20.0,
            source="Homolytic bond cleavage"
        ),
        
        # 取代反应
        ReactionType.SUBSTITUTION: BEPParameters(
            alpha=0.42, beta=18.0,
            source="SN2/SN1 average"
        ),
        
        # 消除反应
        ReactionType.ELIMINATION: BEPParameters(
            alpha=0.45, beta=22.0,
            source="E1/E2 average"
        ),
        
        # 聚合反应
        ReactionType.POLYMERIZATION: BEPParameters(
            alpha=0.35, beta=8.0,
            source="Radical polymerization"
        ),
        
        # Li配位 - 低势垒
        ReactionType.LI_COORDINATION: BEPParameters(
            alpha=0.20, beta=3.0,
            source="Li+ solvation kinetics"
        ),
        
        # Li嵌入
        ReactionType.LI_INSERTION: BEPParameters(
            alpha=0.25, beta=8.0,
            source="Li intercalation"
        ),
        
        # 未知类型 - 保守估计
        ReactionType.UNKNOWN: BEPParameters(
            alpha=0.50, beta=15.0,
            source="Default (conservative)"
        ),
    }
    
    # 环开裂的特殊参数 (依赖环大小)
    RING_OPENING_PARAMETERS: Dict[int, BEPParameters] = {
        3: BEPParameters(0.25, 12.0, "Epoxide opening"),      # 三元环 (高应变~27 kcal/mol)
        4: BEPParameters(0.30, 15.0, "Cyclobutane opening"),  # 四元环 (应变~26 kcal/mol)
        5: BEPParameters(0.35, 18.0, "Carbonate 5-ring"),     # 五元环 (应变~6 kcal/mol)
        6: BEPParameters(0.40, 25.0, "6-membered ring"),      # 六元环 (无应变)
    }
    
    # 键断裂的特殊参数 (依赖键类型)
    BOND_CLEAVAGE_PARAMETERS: Dict[str, BEPParameters] = {
        'C-C': BEPParameters(0.30, 35.0, "C-C BDE ~85 kcal/mol"),
        'C-O': BEPParameters(0.35, 30.0, "C-O BDE ~90 kcal/mol"),
        'C-H': BEPParameters(0.40, 38.0, "C-H BDE ~100 kcal/mol"),
        'O-H': BEPParameters(0.42, 40.0, "O-H BDE ~110 kcal/mol"),
        'P-F': BEPParameters(0.25, 25.0, "P-F BDE ~115 kcal/mol"),
        'O-O': BEPParameters(0.20, 15.0, "O-O peroxide ~35 kcal/mol"),
        'C-F': BEPParameters(0.45, 45.0, "C-F BDE ~115 kcal/mol"),
    }
    
    # 法拉第常数 (用于电化学校正)
    FARADAY_KCAL = 23.06  # kcal/mol/V
    
    def __init__(self, temperature: float = 300.0):
        """
        初始化估算器
        
        Args:
            temperature: 温度 (K)
        """
        self.temperature = temperature
        self.R = 0.001987  # 气体常数 (kcal/mol/K)
    
    def estimate(self, 
                 reaction_type: ReactionType,
                 delta_E: float,
                 ring_size: Optional[int] = None,
                 bond_type: Optional[str] = None,
                 electrode_type: Optional[str] = None,
                 voltage: float = 0.0,
                 n_electrons: int = 1) -> Tuple[float, Dict[str, Any]]:
        """
        估算活化能
        
        Args:
            reaction_type: 反应类型
            delta_E: 反应能 (kcal/mol)
            ring_size: 环大小 (如果是环开裂反应)
            bond_type: 键类型 (如果是键断裂反应, 如'C-O')
            electrode_type: 电极类型 ('anode' 或 'cathode')
            voltage: 电极电势 (V)
            n_electrons: 转移电子数
            
        Returns:
            (活化能, 详细信息字典)
        """
        # 获取BEP参数
        params = self._get_bep_params(reaction_type, ring_size, bond_type)
        alpha, beta = params.alpha, params.beta
        
        # 计算基础活化能
        if delta_E > 0:  # 吸能反应
            Ea = alpha * delta_E + beta
            reaction_nature = "endothermic"
        else:  # 放能反应 (Hammond假说: 过渡态更接近反应物)
            Ea = beta - (1 - alpha) * abs(delta_E)
            Ea = max(Ea, 2.0)  # 最小势垒2 kcal/mol
            reaction_nature = "exothermic"
        
        # 记录详细信息
        details = {
            'reaction_type': reaction_type.value,
            'delta_E': delta_E,
            'alpha': alpha,
            'beta': beta,
            'base_Ea': Ea,
            'reaction_nature': reaction_nature,
            'source': params.source,
        }
        
        # 应用电化学校正
        if electrode_type and voltage != 0:
            Ea_before = Ea
            Ea = self._apply_electrochemical_correction(
                Ea, electrode_type, voltage, n_electrons
            )
            details['electrochemical_correction'] = Ea_before - Ea
            details['electrode_type'] = electrode_type
            details['voltage'] = voltage
        
        # 确保非负
        Ea = max(0, Ea)
        details['final_Ea'] = Ea
        
        return Ea, details
    
    def _get_bep_params(self, 
                        reaction_type: ReactionType,
                        ring_size: Optional[int] = None,
                        bond_type: Optional[str] = None) -> BEPParameters:
        """获取BEP参数"""
        
        if reaction_type == ReactionType.RING_OPENING and ring_size:
            return self.RING_OPENING_PARAMETERS.get(
                ring_size, 
                BEPParameters(0.40, 20.0, f"Ring size {ring_size}")
            )
        
        if reaction_type == ReactionType.BOND_CLEAVAGE and bond_type:
            return self.BOND_CLEAVAGE_PARAMETERS.get(
                bond_type,
                BEPParameters(0.35, 30.0, f"Bond type {bond_type}")
            )
        
        return self.BEP_PARAMETERS.get(
            reaction_type, 
            self.BEP_PARAMETERS[ReactionType.UNKNOWN]
        )
    
    def _apply_electrochemical_correction(self,
                                          Ea: float,
                                          electrode_type: str,
                                          voltage: float,
                                          n_electrons: int) -> float:
        """
        Butler-Volmer电化学校正
        
        对于电极反应，活化能受电势影响:
        Ea = Ea0 - α_transfer * n * F * η
        
        其中:
        - α_transfer: 传递系数 (~0.5 for symmetric barrier)
        - n: 转移电子数
        - F: 法拉第常数 (23.06 kcal/mol/V)
        - η: 过电势 (相对于平衡电势)
        """
        alpha_transfer = 0.5  # 对称势垒假设
        
        if electrode_type == 'cathode':
            # 阴极 (还原): 正电势降低还原势垒
            correction = alpha_transfer * n_electrons * self.FARADAY_KCAL * abs(voltage)
            Ea -= correction
        elif electrode_type == 'anode':
            # 阳极 (氧化): 正电势降低氧化势垒
            correction = (1 - alpha_transfer) * n_electrons * self.FARADAY_KCAL * abs(voltage)
            Ea -= correction
        
        return max(2.0, Ea)  # 最小势垒2 kcal/mol
    
    def estimate_rate_constant(self, Ea: float, A: float = 1e13) -> float:
        """
        估算速率常数 (Arrhenius方程)
        
        k = A * exp(-Ea / RT)
        
        Args:
            Ea: 活化能 (kcal/mol)
            A: 指前因子 (s^-1), 默认10^13
            
        Returns:
            速率常数 (s^-1)
        """
        return A * math.exp(-Ea / (self.R * self.temperature))
    
    def estimate_half_life(self, Ea: float, A: float = 1e13) -> float:
        """
        估算反应半衰期
        
        t_1/2 = ln(2) / k
        
        Args:
            Ea: 活化能 (kcal/mol)
            A: 指前因子 (s^-1)
            
        Returns:
            半衰期 (秒)
        """
        k = self.estimate_rate_constant(Ea, A)
        if k > 0:
            return 0.693 / k
        return float('inf')


class EnergyFilterConfig:
    """
    能量筛选配置 - 基于物理化学原理
    
    特点:
    1. 使用物理上合理的截断值
    2. 考虑温度依赖性
    3. 包含电化学驱动修正
    """
    
    # 热力学筛选 (基于平衡常数)
    # K = exp(-ΔG/RT), ΔG < 40 kcal/mol 对应 K > 10^-30
    THERMODYNAMIC_CUTOFF_DEFAULT = 40.0  # kcal/mol
    
    # 动力学筛选 (基于Arrhenius方程)
    # 要求 k > 10^-10 s^-1 (反应时间 < 300年)
    KINETIC_CUTOFF_MAP = {
        250: 22.0,  # kcal/mol @ 250K
        275: 25.0,  # kcal/mol @ 275K
        300: 28.0,  # kcal/mol @ 300K
        325: 31.0,  # kcal/mol @ 325K
        350: 34.0,  # kcal/mol @ 350K
        400: 40.0,  # kcal/mol @ 400K
    }
    
    # 法拉第常数
    FARADAY_KCAL = 23.06  # kcal/mol/V
    
    def __init__(self, 
                 temperature: float = 300.0,
                 voltage: float = 0.0,
                 electrode_type: str = 'anode',
                 custom_cutoff: Optional[float] = None):
        """
        初始化配置
        
        Args:
            temperature: 温度 (K)
            voltage: 电极电势 (V)
            electrode_type: 电极类型
            custom_cutoff: 自定义截断值
        """
        self.temperature = temperature
        self.voltage = voltage
        self.electrode_type = electrode_type
        self.custom_cutoff = custom_cutoff
    
    def get_thermodynamic_cutoff(self) -> float:
        """获取热力学截断值"""
        if self.custom_cutoff:
            return self.custom_cutoff
        
        base = self.THERMODYNAMIC_CUTOFF_DEFAULT
        
        # 电化学驱动修正
        electrochemical_bonus = abs(self.voltage) * self.FARADAY_KCAL
        
        return base + electrochemical_bonus
    
    def get_kinetic_cutoff(self) -> float:
        """获取动力学截断值"""
        if self.custom_cutoff:
            return self.custom_cutoff
        
        # 找到最接近的温度
        temps = sorted(self.KINETIC_CUTOFF_MAP.keys())
        closest_temp = min(temps, key=lambda t: abs(t - self.temperature))
        
        base = self.KINETIC_CUTOFF_MAP[closest_temp]
        
        # 线性插值修正
        RT_factor = self.temperature / 300.0
        base *= RT_factor
        
        # 电化学驱动修正
        electrochemical_bonus = abs(self.voltage) * self.FARADAY_KCAL * 0.5
        
        return base + electrochemical_bonus
    
    def is_thermodynamically_feasible(self, delta_E: float) -> bool:
        """检查热力学可行性"""
        cutoff = self.get_thermodynamic_cutoff()
        return abs(delta_E) <= cutoff
    
    def is_kinetically_accessible(self, Ea: float) -> bool:
        """检查动力学可及性"""
        cutoff = self.get_kinetic_cutoff()
        return Ea <= cutoff
    
    def is_feasible(self, delta_E: float, Ea: float) -> Tuple[bool, str]:
        """
        综合可行性检查
        
        Returns:
            (是否可行, 原因)
        """
        thermo_ok = self.is_thermodynamically_feasible(delta_E)
        kinetic_ok = self.is_kinetically_accessible(Ea)
        
        if thermo_ok and kinetic_ok:
            return True, "feasible"
        elif not thermo_ok and not kinetic_ok:
            return False, "thermodynamically and kinetically unfeasible"
        elif not thermo_ok:
            return False, f"thermodynamically unfeasible (|ΔE|={abs(delta_E):.1f} > {self.get_thermodynamic_cutoff():.1f})"
        else:
            return False, f"kinetically inaccessible (Ea={Ea:.1f} > {self.get_kinetic_cutoff():.1f})"
