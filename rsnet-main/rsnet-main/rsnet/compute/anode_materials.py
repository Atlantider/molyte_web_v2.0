"""
Anode Material Models for RSNet
负极材料模型 - 包含各种负极材料的电化学特性和SEI形成行为

支持的材料:
- Li Metal: 锂金属负极
- Na Metal: 钠金属负极  
- K Metal: 钾金属负极
- Graphite: 石墨 (锂离子电池主流)
- Silicon: 硅负极 (高容量)
- Hard Carbon: 硬碳 (钠离子电池)
- LTO: 钛酸锂 (Li₄Ti₅O₁₂)

SEI形成机制:
- 金属负极: 直接还原，界面不稳定
- 嵌入型负极: 首圈形成SEI后稳定
- 合金化负极: 体积变化大，SEI持续破裂

作者: RSNet Team
版本: 2.0
"""

from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
from enum import Enum


class AnodeMaterialFamily(Enum):
    """负极材料家族"""
    METAL = "metal"                    # 金属 (Li, Na, K)
    INTERCALATION = "intercalation"    # 嵌入型 (石墨, LTO)
    ALLOY = "alloy"                    # 合金化 (Si, Sn)
    CONVERSION = "conversion"           # 转化型
    HARD_CARBON = "hard_carbon"         # 硬碳 (无序碳)


class CarrierIon(Enum):
    """载流子离子类型"""
    LITHIUM = "Li"
    SODIUM = "Na"
    POTASSIUM = "K"


class SEIFormationType(Enum):
    """SEI形成类型"""
    STABLE = "stable"              # 稳定SEI (石墨)
    CONTINUOUS = "continuous"      # 持续形成 (Si)
    MOSAIC = "mosaic"              # 马赛克结构 (金属)
    MINIMAL = "minimal"            # 极少SEI (LTO)


@dataclass
class SEIFormationProfile:
    """SEI形成特性"""
    onset_voltage: float                      # SEI开始形成电压 (V vs M/M+)
    main_formation_range: Tuple[float, float] # 主要形成电压范围
    formation_type: SEIFormationType          # 形成类型
    primary_components: List[str]             # 主要SEI组分
    is_stable: bool                           # SEI是否稳定
    volume_change_impact: float               # 体积变化影响 (0-1)


@dataclass
class AnodeSpeciesProfile:
    """负极物种注入配置"""
    ion_smiles: str                           # 离子SMILES (如 [Li+])
    atom_smiles: str                          # 原子SMILES (如 [Li])
    reduced_smiles: Optional[str] = None      # 过还原态 (如 [Li-])
    solvated_ion_name: str = ""               # 溶剂化离子名称
    additional_species: List[str] = field(default_factory=list)  # 额外物种


@dataclass
class AnodeMaterial:
    """
    负极材料数据类
    
    包含材料的电化学特性、载流子信息和SEI形成行为
    """
    code: str                                 # 材料代码
    name: str                                 # 完整名称
    formula: str                              # 化学式
    family: AnodeMaterialFamily               # 材料家族
    carrier_ion: CarrierIon                   # 载流子离子
    
    # 电化学窗口
    voltage_min: float                        # 最低工作电压 (vs M/M+)
    voltage_max: float                        # 最高工作电压
    nominal_voltage: float                    # 标称电压
    
    # 容量
    theoretical_capacity: float               # 理论容量 (mAh/g)
    practical_capacity: float                 # 实际容量 (mAh/g)
    
    # 体积变化
    volume_change: float                      # 充放电体积变化 %
    
    # SEI特性
    sei_profile: Optional[SEIFormationProfile] = None
    
    # 物种配置
    species_profile: Optional[AnodeSpeciesProfile] = None
    
    # 描述
    description: str = ""
    
    def get_species_to_inject(self, voltage: float) -> List[str]:
        """
        根据电压获取应该注入的物种
        
        Args:
            voltage: 电极电势 (V vs M/M+)
            
        Returns:
            物种SMILES列表
        """
        species = []
        
        if self.species_profile:
            # 总是添加离子和原子
            species.append(self.species_profile.ion_smiles)
            species.append(self.species_profile.atom_smiles)
            
            # 低电压时添加还原态
            if voltage < 0.1 and self.species_profile.reduced_smiles:
                species.append(self.species_profile.reduced_smiles)
            
            # 添加材料特定物种
            species.extend(self.species_profile.additional_species)
        
        return species
    
    def is_in_sei_formation_range(self, voltage: float) -> bool:
        """检查是否在SEI形成电压范围"""
        if not self.sei_profile:
            return False
        low, high = self.sei_profile.main_formation_range
        return low <= voltage <= high


# =============================================================================
# 物种配置
# =============================================================================

LITHIUM_SPECIES = AnodeSpeciesProfile(
    ion_smiles='[Li+]',
    atom_smiles='[Li]',
    reduced_smiles='[Li-]',
    solvated_ion_name='溶剂化锂离子',
    additional_species=[]
)

SODIUM_SPECIES = AnodeSpeciesProfile(
    ion_smiles='[Na+]',
    atom_smiles='[Na]',
    reduced_smiles='[Na-]',
    solvated_ion_name='溶剂化钠离子',
    additional_species=[]
)

POTASSIUM_SPECIES = AnodeSpeciesProfile(
    ion_smiles='[K+]',
    atom_smiles='[K]',
    reduced_smiles='[K-]',
    solvated_ion_name='溶剂化钾离子',
    additional_species=[]
)

SILICON_SPECIES = AnodeSpeciesProfile(
    ion_smiles='[Li+]',
    atom_smiles='[Li]',
    reduced_smiles=None,
    solvated_ion_name='溶剂化锂离子',
    additional_species=['[Si]', '[SiH4]']  # 硅和硅烷
)

GRAPHITE_SPECIES = AnodeSpeciesProfile(
    ion_smiles='[Li+]',
    atom_smiles='[Li]',
    reduced_smiles=None,
    solvated_ion_name='溶剂化锂离子',
    additional_species=[]  # 石墨不释放额外物种
)


# =============================================================================
# 预定义的负极材料库
# =============================================================================

ANODE_MATERIALS: Dict[str, AnodeMaterial] = {
    # =========================================================================
    # 锂系负极
    # =========================================================================
    'LI_METAL': AnodeMaterial(
        code='LI_METAL',
        name='锂金属 (Lithium Metal)',
        formula='Li',
        family=AnodeMaterialFamily.METAL,
        carrier_ion=CarrierIon.LITHIUM,
        voltage_min=0.0,
        voltage_max=0.0,
        nominal_voltage=0.0,
        theoretical_capacity=3860,
        practical_capacity=3000,
        volume_change=100,  # 无限体积变化
        sei_profile=SEIFormationProfile(
            onset_voltage=1.0,
            main_formation_range=(0.0, 0.8),
            formation_type=SEIFormationType.MOSAIC,
            primary_components=['Li2CO3', 'LiF', 'Li2O', 'ROLi'],
            is_stable=False,  # 枝晶导致SEI破裂
            volume_change_impact=1.0
        ),
        species_profile=LITHIUM_SPECIES,
        description='最高容量但枝晶和SEI稳定性问题'
    ),
    
    'GRAPHITE': AnodeMaterial(
        code='GRAPHITE',
        name='石墨 (Graphite)',
        formula='C6 (LiC6)',
        family=AnodeMaterialFamily.INTERCALATION,
        carrier_ion=CarrierIon.LITHIUM,
        voltage_min=0.01,
        voltage_max=0.3,
        nominal_voltage=0.1,
        theoretical_capacity=372,
        practical_capacity=350,
        volume_change=10,  # 10%体积变化
        sei_profile=SEIFormationProfile(
            onset_voltage=0.8,
            main_formation_range=(0.1, 0.8),
            formation_type=SEIFormationType.STABLE,
            primary_components=['Li2CO3', 'LiF', 'LEDC', '(CH2OCO2Li)2'],
            is_stable=True,
            volume_change_impact=0.1
        ),
        species_profile=GRAPHITE_SPECIES,
        description='锂离子电池主流负极，SEI稳定'
    ),
    
    'SILICON': AnodeMaterial(
        code='SILICON',
        name='硅 (Silicon)',
        formula='Si (Li15Si4)',
        family=AnodeMaterialFamily.ALLOY,
        carrier_ion=CarrierIon.LITHIUM,
        voltage_min=0.01,
        voltage_max=0.5,
        nominal_voltage=0.2,
        theoretical_capacity=4200,
        practical_capacity=2000,
        volume_change=300,  # 300%体积变化!
        sei_profile=SEIFormationProfile(
            onset_voltage=1.0,
            main_formation_range=(0.1, 1.0),
            formation_type=SEIFormationType.CONTINUOUS,
            primary_components=['Li2CO3', 'LiF', 'Li2SiO3', 'SiOx'],
            is_stable=False,  # 体积变化导致SEI持续破裂
            volume_change_impact=1.0
        ),
        species_profile=SILICON_SPECIES,
        description='超高容量但循环稳定性差'
    ),
    
    'LTO': AnodeMaterial(
        code='LTO',
        name='钛酸锂 (Lithium Titanate)',
        formula='Li4Ti5O12',
        family=AnodeMaterialFamily.INTERCALATION,
        carrier_ion=CarrierIon.LITHIUM,
        voltage_min=1.5,
        voltage_max=1.6,
        nominal_voltage=1.55,
        theoretical_capacity=175,
        practical_capacity=165,
        volume_change=0,  # 零应变材料
        sei_profile=SEIFormationProfile(
            onset_voltage=0.8,  # 高于工作电压
            main_formation_range=(0.5, 0.8),
            formation_type=SEIFormationType.MINIMAL,
            primary_components=[],  # 几乎无SEI
            is_stable=True,
            volume_change_impact=0.0
        ),
        species_profile=LITHIUM_SPECIES,
        description='零应变材料，几乎无SEI，安全性极高'
    ),
    
    'SIC': AnodeMaterial(
        code='SIC',
        name='硅碳复合 (Si/C Composite)',
        formula='Si-C',
        family=AnodeMaterialFamily.ALLOY,
        carrier_ion=CarrierIon.LITHIUM,
        voltage_min=0.01,
        voltage_max=0.5,
        nominal_voltage=0.15,
        theoretical_capacity=1500,
        practical_capacity=800,
        volume_change=120,
        sei_profile=SEIFormationProfile(
            onset_voltage=1.0,
            main_formation_range=(0.1, 1.0),
            formation_type=SEIFormationType.CONTINUOUS,
            primary_components=['Li2CO3', 'LiF', 'Li2SiO3', 'LEDC'],
            is_stable=False,
            volume_change_impact=0.6
        ),
        species_profile=SILICON_SPECIES,
        description='平衡容量和稳定性的复合材料'
    ),
    
    # =========================================================================
    # 钠系负极
    # =========================================================================
    'NA_METAL': AnodeMaterial(
        code='NA_METAL',
        name='钠金属 (Sodium Metal)',
        formula='Na',
        family=AnodeMaterialFamily.METAL,
        carrier_ion=CarrierIon.SODIUM,
        voltage_min=0.0,
        voltage_max=0.0,
        nominal_voltage=0.0,
        theoretical_capacity=1166,
        practical_capacity=800,
        volume_change=100,
        sei_profile=SEIFormationProfile(
            onset_voltage=1.2,
            main_formation_range=(0.0, 1.0),
            formation_type=SEIFormationType.MOSAIC,
            primary_components=['Na2CO3', 'NaF', 'Na2O', 'RONa'],
            is_stable=False,
            volume_change_impact=1.0
        ),
        species_profile=SODIUM_SPECIES,
        description='钠离子电池金属负极'
    ),
    
    'HARD_CARBON': AnodeMaterial(
        code='HARD_CARBON',
        name='硬碳 (Hard Carbon)',
        formula='C (NaCx)',
        family=AnodeMaterialFamily.HARD_CARBON,
        carrier_ion=CarrierIon.SODIUM,
        voltage_min=0.0,
        voltage_max=0.5,
        nominal_voltage=0.1,
        theoretical_capacity=350,
        practical_capacity=300,
        volume_change=8,
        sei_profile=SEIFormationProfile(
            onset_voltage=1.0,
            main_formation_range=(0.1, 0.8),
            formation_type=SEIFormationType.STABLE,
            primary_components=['Na2CO3', 'NaF', 'LEDC_Na', 'RONa'],
            is_stable=True,
            volume_change_impact=0.1
        ),
        species_profile=SODIUM_SPECIES,
        description='钠离子电池主流负极'
    ),
    
    'SOFT_CARBON': AnodeMaterial(
        code='SOFT_CARBON',
        name='软碳 (Soft Carbon)',
        formula='C',
        family=AnodeMaterialFamily.HARD_CARBON,
        carrier_ion=CarrierIon.SODIUM,
        voltage_min=0.1,
        voltage_max=0.6,
        nominal_voltage=0.2,
        theoretical_capacity=250,
        practical_capacity=200,
        volume_change=10,
        sei_profile=SEIFormationProfile(
            onset_voltage=1.0,
            main_formation_range=(0.2, 0.8),
            formation_type=SEIFormationType.STABLE,
            primary_components=['Na2CO3', 'NaF'],
            is_stable=True,
            volume_change_impact=0.1
        ),
        species_profile=SODIUM_SPECIES,
        description='可石墨化碳材料'
    ),
    
    # =========================================================================
    # 钾系负极
    # =========================================================================
    'K_METAL': AnodeMaterial(
        code='K_METAL',
        name='钾金属 (Potassium Metal)',
        formula='K',
        family=AnodeMaterialFamily.METAL,
        carrier_ion=CarrierIon.POTASSIUM,
        voltage_min=0.0,
        voltage_max=0.0,
        nominal_voltage=0.0,
        theoretical_capacity=685,
        practical_capacity=500,
        volume_change=100,
        sei_profile=SEIFormationProfile(
            onset_voltage=1.5,
            main_formation_range=(0.0, 1.2),
            formation_type=SEIFormationType.MOSAIC,
            primary_components=['K2CO3', 'KF', 'K2O', 'ROK'],
            is_stable=False,
            volume_change_impact=1.0
        ),
        species_profile=POTASSIUM_SPECIES,
        description='钾离子电池金属负极'
    ),
    
    'K_GRAPHITE': AnodeMaterial(
        code='K_GRAPHITE',
        name='钾石墨 (Potassium Graphite)',
        formula='KC8',
        family=AnodeMaterialFamily.INTERCALATION,
        carrier_ion=CarrierIon.POTASSIUM,
        voltage_min=0.1,
        voltage_max=0.5,
        nominal_voltage=0.2,
        theoretical_capacity=279,
        practical_capacity=250,
        volume_change=60,  # 比锂石墨大得多
        sei_profile=SEIFormationProfile(
            onset_voltage=1.2,
            main_formation_range=(0.2, 1.0),
            formation_type=SEIFormationType.STABLE,
            primary_components=['K2CO3', 'KF', 'LEDC_K'],
            is_stable=True,
            volume_change_impact=0.4
        ),
        species_profile=POTASSIUM_SPECIES,
        description='钾可嵌入石墨，但体积变化较大'
    ),
}


def get_anode_material(code: str) -> Optional[AnodeMaterial]:
    """
    获取负极材料信息
    
    Args:
        code: 材料代码 (不区分大小写)
        
    Returns:
        AnodeMaterial对象或None
    """
    return ANODE_MATERIALS.get(code.upper())


def list_available_anode_materials() -> List[str]:
    """列出所有可用的负极材料"""
    return list(ANODE_MATERIALS.keys())


def get_materials_by_carrier(carrier: CarrierIon) -> List[AnodeMaterial]:
    """按载流子离子获取材料列表"""
    return [m for m in ANODE_MATERIALS.values() if m.carrier_ion == carrier]


def print_anode_material_summary(code: str):
    """打印材料摘要信息"""
    material = get_anode_material(code)
    if not material:
        print(f"未知材料: {code}")
        return
    
    print(f"\n{'='*60}")
    print(f"材料: {material.name}")
    print(f"化学式: {material.formula}")
    print(f"家族: {material.family.value}")
    print(f"载流子: {material.carrier_ion.value}⁺")
    print(f"{'='*60}")
    print(f"电压范围: {material.voltage_min} - {material.voltage_max} V")
    print(f"理论容量: {material.theoretical_capacity} mAh/g")
    print(f"实际容量: {material.practical_capacity} mAh/g")
    print(f"体积变化: {material.volume_change}%")
    
    if material.sei_profile:
        print(f"\nSEI形成特性:")
        print(f"  类型: {material.sei_profile.formation_type.value}")
        print(f"  形成范围: {material.sei_profile.main_formation_range[0]} - {material.sei_profile.main_formation_range[1]} V")
        print(f"  稳定性: {'稳定' if material.sei_profile.is_stable else '不稳定'}")
        print(f"  主要组分: {', '.join(material.sei_profile.primary_components)}")
    
    print(f"\n描述: {material.description}")
    print(f"{'='*60}\n")
