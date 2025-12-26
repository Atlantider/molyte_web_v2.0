"""
Cathode Material Models for RSNet
正极材料模型 - 包含各种正极材料的电化学特性和氧释放行为

支持的材料:
- LCO (LiCoO₂): 层状钴酸锂
- NMC (LiNi₀.₈Mn₀.₁Co₀.₁O₂): 三元材料
- NCA (LiNi₀.₈Co₀.₁₅Al₀.₀₅O₂): 镍钴铝
- LFP (LiFePO₄): 磷酸铁锂
- LMO (LiMn₂O₄): 尖晶石锰酸锂
- LNMO (LiNi₀.₅Mn₁.₅O₄): 高电压尖晶石
- LCO_HV: 高电压钴酸锂

氧释放机制:
- Lattice oxygen: 晶格氧脱出
- Oxygen redox: 氧参与电子转移
- Surface reconstruction: 表面重构

作者: RSNet Team
版本: 2.0
"""

from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
from enum import Enum


class CathodeMaterialFamily(Enum):
    """正极材料家族"""
    LAYERED = "layered"           # 层状材料 (LCO, NMC, NCA)
    SPINEL = "spinel"             # 尖晶石 (LMO, LNMO)
    OLIVINE = "olivine"           # 橄榄石 (LFP)
    POLYANION = "polyanion"       # 聚阴离子 (LFMP, LVPF)


class OxygenReleaseType(Enum):
    """氧释放类型"""
    NONE = "none"                 # 无释放 (如LFP)
    LATTICE = "lattice"           # 晶格氧释放
    SURFACE = "surface"           # 表面氧释放
    OXYGEN_REDOX = "oxygen_redox" # 氧参与电子转移


@dataclass
class OxygenReleaseProfile:
    """氧释放特性"""
    onset_voltage: float                      # 起始电压 (V vs Li/Li+)
    peak_voltage: float                       # 峰值电压
    release_type: OxygenReleaseType           # 释放类型
    primary_species: List[str]                # 主要释放物种
    secondary_species: List[str]              # 次要释放物种
    relative_rate: float                      # 相对释放速率 (0-1)


@dataclass
class CathodeMaterial:
    """
    正极材料数据类
    
    包含材料的电化学特性、结构信息和氧释放行为
    """
    code: str                                 # 材料代码 (NMC, LCO等)
    name: str                                 # 完整名称
    formula: str                              # 化学式
    family: CathodeMaterialFamily             # 材料家族
    
    # 电化学窗口
    voltage_min: float                        # 最低工作电压
    voltage_max: float                        # 最高工作电压
    nominal_voltage: float                    # 标称电压
    
    # 容量
    theoretical_capacity: float               # 理论容量 (mAh/g)
    practical_capacity: float                 # 实际容量 (mAh/g)
    
    # 氧释放特性
    oxygen_release: Optional[OxygenReleaseProfile] = None
    
    # 过渡金属溶出
    metal_dissolution_species: List[str] = field(default_factory=list)
    metal_dissolution_threshold: float = 4.5  # 溶出阈值电压
    
    # 稳定性
    thermal_stability: float = 0.5            # 热稳定性 (0-1)
    structural_stability: float = 0.5         # 结构稳定性 (0-1)
    
    # 描述
    description: str = ""
    
    def get_released_species(self, voltage: float) -> List[str]:
        """
        根据电压获取释放的物种
        
        Args:
            voltage: 电极电势 (V vs Li/Li+)
            
        Returns:
            释放物种的key列表
        """
        species = []
        
        # 氧释放
        if self.oxygen_release:
            if voltage >= self.oxygen_release.onset_voltage:
                species.extend(self.oxygen_release.primary_species)
                
                if voltage >= self.oxygen_release.peak_voltage:
                    species.extend(self.oxygen_release.secondary_species)
        
        # 金属溶出
        if voltage >= self.metal_dissolution_threshold:
            species.extend(self.metal_dissolution_species)
        
        return list(set(species))  # 去重
    
    def is_oxygen_releasing(self, voltage: float) -> bool:
        """检查是否在氧释放电压范围"""
        if not self.oxygen_release:
            return False
        return voltage >= self.oxygen_release.onset_voltage
    
    def get_stability_score(self, voltage: float) -> float:
        """
        计算指定电压下的稳定性评分
        
        Args:
            voltage: 电压
            
        Returns:
            0-1的稳定性评分
        """
        base_score = (self.thermal_stability + self.structural_stability) / 2
        
        # 超过最高工作电压时降低稳定性
        if voltage > self.voltage_max:
            overvoltage = voltage - self.voltage_max
            penalty = min(0.5, overvoltage * 0.1)  # 每0.1V降10%
            base_score -= penalty
        
        return max(0, base_score)


# =============================================================================
# 预定义的正极材料库
# =============================================================================

CATHODE_MATERIALS: Dict[str, CathodeMaterial] = {
    'LCO': CathodeMaterial(
        code='LCO',
        name='锂钴氧化物 (Lithium Cobalt Oxide)',
        formula='LiCoO₂',
        family=CathodeMaterialFamily.LAYERED,
        voltage_min=3.0,
        voltage_max=4.2,
        nominal_voltage=3.9,
        theoretical_capacity=274,
        practical_capacity=155,
        oxygen_release=OxygenReleaseProfile(
            onset_voltage=4.3,
            peak_voltage=4.5,
            release_type=OxygenReleaseType.LATTICE,
            primary_species=['O-', 'O2-_radical'],  # 超氧根
            secondary_species=['O2', '1O2', 'O2^2-'],  # 过氧根
            relative_rate=0.7
        ),
        metal_dissolution_species=['Co2+', 'Co3+'],
        metal_dissolution_threshold=4.4,
        thermal_stability=0.4,
        structural_stability=0.5,
        description='最早商业化的正极材料，能量密度高但稳定性一般'
    ),
    
    'NMC': CathodeMaterial(
        code='NMC',
        name='三元材料 (Nickel Manganese Cobalt)',
        formula='LiNi₀.₈Mn₀.₁Co₀.₁O₂',
        family=CathodeMaterialFamily.LAYERED,
        voltage_min=3.0,
        voltage_max=4.3,
        nominal_voltage=3.8,
        theoretical_capacity=275,
        practical_capacity=200,
        oxygen_release=OxygenReleaseProfile(
            onset_voltage=4.4,
            peak_voltage=4.6,
            release_type=OxygenReleaseType.OXYGEN_REDOX,
            primary_species=['O-', 'O2-_radical', 'O2^2-'],  # 超氧根+过氧根
            secondary_species=['1O2', 'O2'],
            relative_rate=0.8
        ),
        metal_dissolution_species=['Ni2+', 'Ni3+', 'Mn2+', 'Co2+'],
        metal_dissolution_threshold=4.5,
        thermal_stability=0.45,
        structural_stability=0.55,
        description='主流动力电池正极，高镍配方性能优异但热稳定性需关注'
    ),
    
    'NMC622': CathodeMaterial(
        code='NMC622',
        name='NMC622',
        formula='LiNi₀.₆Mn₀.₂Co₀.₂O₂',
        family=CathodeMaterialFamily.LAYERED,
        voltage_min=3.0,
        voltage_max=4.3,
        nominal_voltage=3.75,
        theoretical_capacity=270,
        practical_capacity=175,
        oxygen_release=OxygenReleaseProfile(
            onset_voltage=4.5,
            peak_voltage=4.7,
            release_type=OxygenReleaseType.OXYGEN_REDOX,
            primary_species=['O-', 'O2-_radical'],
            secondary_species=['O2^2-', '1O2'],
            relative_rate=0.6
        ),
        metal_dissolution_species=['Ni2+', 'Mn2+', 'Co2+'],
        metal_dissolution_threshold=4.5,
        thermal_stability=0.5,
        structural_stability=0.6,
        description='中镍三元材料，性能与安全性平衡'
    ),
    
    'NMC811': CathodeMaterial(
        code='NMC811',
        name='NMC811 (高镍三元)',
        formula='LiNi₀.₈Mn₀.₁Co₀.₁O₂',
        family=CathodeMaterialFamily.LAYERED,
        voltage_min=3.0,
        voltage_max=4.2,
        nominal_voltage=3.85,
        theoretical_capacity=280,
        practical_capacity=210,
        oxygen_release=OxygenReleaseProfile(
            onset_voltage=4.3,
            peak_voltage=4.5,
            release_type=OxygenReleaseType.OXYGEN_REDOX,
            primary_species=['O-', 'O2-_radical', 'O2^2-'],
            secondary_species=['1O2', 'O2', 'Ni4+'],
            relative_rate=0.9  # 高镍释放更多
        ),
        metal_dissolution_species=['Ni2+', 'Ni3+', 'Ni4+', 'Mn2+', 'Co2+'],
        metal_dissolution_threshold=4.4,
        thermal_stability=0.35,
        structural_stability=0.45,
        description='高镍正极，最高能量密度但热稳定性差'
    ),
    
    'NCA': CathodeMaterial(
        code='NCA',
        name='镍钴铝酸锂 (Nickel Cobalt Aluminum)',
        formula='LiNi₀.₈Co₀.₁₅Al₀.₀₅O₂',
        family=CathodeMaterialFamily.LAYERED,
        voltage_min=3.0,
        voltage_max=4.2,
        nominal_voltage=3.8,
        theoretical_capacity=279,
        practical_capacity=200,
        oxygen_release=OxygenReleaseProfile(
            onset_voltage=4.4,
            peak_voltage=4.6,
            release_type=OxygenReleaseType.OXYGEN_REDOX,
            primary_species=['O-', 'O2-_radical'],
            secondary_species=['1O2', 'O2^2-'],
            relative_rate=0.75
        ),
        metal_dissolution_species=['Ni2+', 'Ni3+', 'Co2+'],
        metal_dissolution_threshold=4.5,
        thermal_stability=0.4,
        structural_stability=0.5,
        description='特斯拉采用的正极材料，Al掺杂提高结构稳定性'
    ),
    
    'LFP': CathodeMaterial(
        code='LFP',
        name='磷酸铁锂 (Lithium Iron Phosphate)',
        formula='LiFePO₄',
        family=CathodeMaterialFamily.OLIVINE,
        voltage_min=2.5,
        voltage_max=3.65,
        nominal_voltage=3.4,
        theoretical_capacity=170,
        practical_capacity=160,
        oxygen_release=None,  # LFP不释放氧！
        metal_dissolution_species=['Fe2+'],  # 极少
        metal_dissolution_threshold=4.0,  # 很高阈值
        thermal_stability=0.95,  # 非常稳定
        structural_stability=0.9,
        description='最安全的正极材料，PO₄骨架锁定氧，无氧释放'
    ),
    
    'LMO': CathodeMaterial(
        code='LMO',
        name='锰酸锂 (Lithium Manganese Oxide)',
        formula='LiMn₂O₄',
        family=CathodeMaterialFamily.SPINEL,
        voltage_min=3.0,
        voltage_max=4.2,
        nominal_voltage=4.0,
        theoretical_capacity=148,
        practical_capacity=120,
        oxygen_release=OxygenReleaseProfile(
            onset_voltage=4.3,
            peak_voltage=4.5,
            release_type=OxygenReleaseType.SURFACE,
            primary_species=['O-'],
            secondary_species=['O2-_radical', 'Mn2+'],
            relative_rate=0.5
        ),
        metal_dissolution_species=['Mn2+', 'Mn3+'],
        metal_dissolution_threshold=4.0,  # Mn歧化导致较早溶出
        thermal_stability=0.6,
        structural_stability=0.55,
        description='尖晶石结构，功率性能好但存在Mn溶出问题'
    ),
    
    'LNMO': CathodeMaterial(
        code='LNMO',
        name='高电压尖晶石 (High-Voltage Spinel)',
        formula='LiNi₀.₅Mn₁.₅O₄',
        family=CathodeMaterialFamily.SPINEL,
        voltage_min=4.4,
        voltage_max=4.9,
        nominal_voltage=4.7,
        theoretical_capacity=147,
        practical_capacity=130,
        oxygen_release=OxygenReleaseProfile(
            onset_voltage=4.85,
            peak_voltage=5.0,
            release_type=OxygenReleaseType.LATTICE,
            primary_species=['O-', 'O2-_radical', 'O2^2-', '1O2'],
            secondary_species=['O2'],
            relative_rate=0.85
        ),
        metal_dissolution_species=['Ni2+', 'Mn2+'],
        metal_dissolution_threshold=4.7,
        thermal_stability=0.5,
        structural_stability=0.6,
        description='高电压正极材料，需要特殊电解液'
    ),
    
    'LRLO': CathodeMaterial(
        code='LRLO',
        name='富锂锰基正极 (Li-Rich Layered Oxide)',
        formula='Li₁.₂Ni₀.₁₃Mn₀.₅₄Co₀.₁₃O₂',
        family=CathodeMaterialFamily.LAYERED,
        voltage_min=2.0,
        voltage_max=4.8,
        nominal_voltage=3.5,
        theoretical_capacity=350,
        practical_capacity=250,
        oxygen_release=OxygenReleaseProfile(
            onset_voltage=4.4,
            peak_voltage=4.6,
            release_type=OxygenReleaseType.OXYGEN_REDOX,
            primary_species=['O-', 'O2-', 'O2^2-', 'O2-_radical'],  # 氧参与电子转移
            secondary_species=['1O2', 'O2'],
            relative_rate=1.0  # 氧活性最高
        ),
        metal_dissolution_species=['Mn2+', 'Ni2+', 'Co2+'],
        metal_dissolution_threshold=4.5,
        thermal_stability=0.3,
        structural_stability=0.35,
        description='超高容量但电压衰减严重，氧参与充放电'
    ),
}


def get_cathode_material(code: str) -> Optional[CathodeMaterial]:
    """
    获取正极材料信息
    
    Args:
        code: 材料代码 (不区分大小写)
        
    Returns:
        CathodeMaterial对象或None
    """
    return CATHODE_MATERIALS.get(code.upper())


def list_available_materials() -> List[str]:
    """列出所有可用的正极材料"""
    return list(CATHODE_MATERIALS.keys())


def get_materials_by_family(family: CathodeMaterialFamily) -> List[CathodeMaterial]:
    """按材料家族获取材料列表"""
    return [m for m in CATHODE_MATERIALS.values() if m.family == family]


def print_material_summary(code: str):
    """打印材料摘要信息"""
    material = get_cathode_material(code)
    if not material:
        print(f"未知材料: {code}")
        return
    
    print(f"\n{'='*60}")
    print(f"材料: {material.name}")
    print(f"化学式: {material.formula}")
    print(f"家族: {material.family.value}")
    print(f"{'='*60}")
    print(f"电压范围: {material.voltage_min} - {material.voltage_max} V")
    print(f"标称电压: {material.nominal_voltage} V")
    print(f"理论容量: {material.theoretical_capacity} mAh/g")
    print(f"实际容量: {material.practical_capacity} mAh/g")
    print(f"热稳定性: {material.thermal_stability:.0%}")
    print(f"结构稳定性: {material.structural_stability:.0%}")
    
    if material.oxygen_release:
        print(f"\n氧释放特性:")
        print(f"  类型: {material.oxygen_release.release_type.value}")
        print(f"  起始电压: {material.oxygen_release.onset_voltage} V")
        print(f"  峰值电压: {material.oxygen_release.peak_voltage} V")
        print(f"  主要物种: {', '.join(material.oxygen_release.primary_species)}")
        print(f"  次要物种: {', '.join(material.oxygen_release.secondary_species)}")
    else:
        print(f"\n氧释放: 无 (材料稳定)")
    
    print(f"\n金属溶出:")
    print(f"  阈值电压: {material.metal_dissolution_threshold} V")
    print(f"  溶出物种: {', '.join(material.metal_dissolution_species)}")
    
    print(f"\n描述: {material.description}")
    print(f"{'='*60}\n")
