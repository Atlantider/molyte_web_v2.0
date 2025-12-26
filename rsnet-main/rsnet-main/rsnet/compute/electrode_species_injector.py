"""
Electrode Species Injector for RSNet
电极物种自动注入器 - 根据电极类型自动添加关键反应物种

包含:
- 阳极物种: Li⁺, Li, e⁻
- 阴极物种: O²⁻, O⁻, O₂²⁻(过氧根), O₂⁻(超氧根), O₂, ¹O₂, ·OH, 过渡金属离子

化学背景:
- 阳极SEI: Li⁺ + e⁻ + 溶剂 → SEI (还原反应)
- 阴极CEI: 溶剂 + O·/O₂⁻ → CEI (氧化反应)

过氧根与超氧根:
- O₂²⁻ (过氧根): 两个额外电子, O-O键弱, 强亲核试剂
- O₂⁻ (超氧根): 一个额外电子, 自由基性质, 非常活泼

作者: RSNet Team
版本: 2.0
"""

from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
from enum import Enum


class ElectrodeType(Enum):
    """电极类型"""
    ANODE = "anode"      # 阳极 (负极) - 还原环境
    CATHODE = "cathode"  # 阴极 (正极) - 氧化环境


@dataclass
class ElectrodeSpecies:
    """电极物种数据类"""
    name: str           # 物种名称
    smiles: str         # SMILES表示
    charge: int         # 电荷
    role: str           # 在界面反应中的角色
    origin: str         # 物种来源
    reactivity: float   # 反应活性 (0-1)
    conditions: List[str] = field(default_factory=list)  # 生成条件


class ElectrodeSpeciesInjector:
    """
    电极物种自动注入器
    
    根据电极类型、材料和电压自动添加关键反应物种，
    确保反应网络包含所有相关的界面化学物种。
    """
    
    # ============================================================
    # 阳极(负极)关键物种 - 按载流子类型分类
    # ============================================================
    
    # 锂系负极物种
    LITHIUM_ANODE_SPECIES: Dict[str, ElectrodeSpecies] = {
        'Li+': ElectrodeSpecies(
            name='锂离子',
            smiles='[Li+]',
            charge=1,
            role='SEI形成的核心物种，与溶剂还原产物反应',
            origin='电解液盐 (LiPF6, LiTFSI等)',
            reactivity=0.9,
            conditions=['always']
        ),
        'Li': ElectrodeSpecies(
            name='锂原子',
            smiles='[Li]',
            charge=0,
            role='Li⁺还原产物，参与SEI形成',
            origin='Li⁺ + e⁻ → Li',
            reactivity=1.0,
            conditions=['low_voltage']
        ),
        'Li-': ElectrodeSpecies(
            name='锂阴离子',
            smiles='[Li-]',
            charge=-1,
            role='极端还原条件下的锂物种',
            origin='Li + e⁻ → Li⁻',
            reactivity=1.0,
            conditions=['very_low_voltage']
        ),
    }
    
    # 钠系负极物种
    SODIUM_ANODE_SPECIES: Dict[str, ElectrodeSpecies] = {
        'Na+': ElectrodeSpecies(
            name='钠离子',
            smiles='[Na+]',
            charge=1,
            role='SEI形成核心物种，钠离子电池',
            origin='电解液盐 (NaPF6, NaClO4等)',
            reactivity=0.85,
            conditions=['always']
        ),
        'Na': ElectrodeSpecies(
            name='钠原子',
            smiles='[Na]',
            charge=0,
            role='Na⁺还原产物',
            origin='Na⁺ + e⁻ → Na',
            reactivity=0.95,
            conditions=['low_voltage']
        ),
        'Na-': ElectrodeSpecies(
            name='钠阴离子',
            smiles='[Na-]',
            charge=-1,
            role='极端还原条件下的钠物种',
            origin='Na + e⁻ → Na⁻',
            reactivity=1.0,
            conditions=['very_low_voltage']
        ),
    }
    
    # 钾系负极物种
    POTASSIUM_ANODE_SPECIES: Dict[str, ElectrodeSpecies] = {
        'K+': ElectrodeSpecies(
            name='钾离子',
            smiles='[K+]',
            charge=1,
            role='SEI形成核心物种，钾离子电池',
            origin='电解液盐 (KPF6等)',
            reactivity=0.8,
            conditions=['always']
        ),
        'K': ElectrodeSpecies(
            name='钾原子',
            smiles='[K]',
            charge=0,
            role='K⁺还原产物',
            origin='K⁺ + e⁻ → K',
            reactivity=0.9,
            conditions=['low_voltage']
        ),
        'K-': ElectrodeSpecies(
            name='钾阴离子',
            smiles='[K-]',
            charge=-1,
            role='极端还原条件下的钾物种',
            origin='K + e⁻ → K⁻',
            reactivity=1.0,
            conditions=['very_low_voltage']
        ),
    }
    
    # 硅负极特有物种
    SILICON_ANODE_SPECIES: Dict[str, ElectrodeSpecies] = {
        'SiH4': ElectrodeSpecies(
            name='硅烷',
            smiles='[SiH4]',
            charge=0,
            role='硅负极还原分解产物',
            origin='Si + 4H⁺ + 4e⁻ → SiH₄',
            reactivity=0.6,
            conditions=['silicon', 'low_voltage']
        ),
        'Si': ElectrodeSpecies(
            name='硅原子',
            smiles='[Si]',
            charge=0,
            role='合金化反应中心',
            origin='硅负极活性材料',
            reactivity=0.7,
            conditions=['silicon']
        ),
        'LixSi': ElectrodeSpecies(
            name='锂硅合金',
            smiles='[Li][Si]',
            charge=0,
            role='锂硅合金中间体',
            origin='xLi + Si → LixSi',
            reactivity=0.8,
            conditions=['silicon', 'low_voltage']
        ),
    }
    
    # 向后兼容: 默认阳极物种 (锂)
    ANODE_SPECIES = LITHIUM_ANODE_SPECIES
    

    # ============================================================
    # 阴极(正极)关键物种 - 氧化性环境
    # ============================================================
    CATHODE_SPECIES: Dict[str, ElectrodeSpecies] = {
        # --- 基础氧物种 ---
        'O2-': ElectrodeSpecies(
            name='氧离子 (晶格氧)',
            smiles='[O-2]',
            charge=-2,
            role='正极材料晶格氧，强还原剂/亲核试剂',
            origin='LiMO₂晶格 (M=Co,Ni,Mn)',
            reactivity=0.8,
            conditions=['always']
        ),
        'O-': ElectrodeSpecies(
            name='氧自由基阴离子',
            smiles='[O-]',
            charge=-1,
            role='单电子氧化产物，活泼亲电试剂',
            origin='O²⁻ - e⁻ → O⁻',
            reactivity=0.95,
            conditions=['high_voltage']  # > 4.0V
        ),
        'O': ElectrodeSpecies(
            name='氧原子/单线态氧',
            smiles='[O]',
            charge=0,
            role='高反应性氧物种',
            origin='高电压下晶格氧释放',
            reactivity=1.0,
            conditions=['very_high_voltage']  # > 4.5V
        ),
        
        # --- 过氧根 O₂²⁻ (Peroxide) ---
        'O2^2-': ElectrodeSpecies(
            name='过氧根 (Peroxide)',
            smiles='[O-][O-]',
            charge=-2,
            role='O-O键弱的二氧化物种，强亲核试剂，可分解为O²⁻',
            origin='2O²⁻ - 2e⁻ → O₂²⁻ 或 晶格缺陷',
            reactivity=0.85,
            conditions=['high_voltage', 'li2o2_formation']
        ),
        
        # --- 超氧根 O₂⁻ (Superoxide) ---
        'O2-_radical': ElectrodeSpecies(
            name='超氧根 (Superoxide)',
            smiles='[O][O-]',
            charge=-1,
            role='自由基性质，高反应性氧化剂，攻击C-H键',
            origin='O₂ + e⁻ → O₂⁻ 或 O₂²⁻ - e⁻ → O₂⁻',
            reactivity=0.98,
            conditions=['high_voltage', 'oxygen_release']
        ),
        
        # --- 分子氧 ---
        'O2': ElectrodeSpecies(
            name='分子氧 (三重态)',
            smiles='O=O',
            charge=0,
            role='正极材料释放的氧气，可进一步活化',
            origin='高电压下材料分解: LiMO₂ → LiM₁₋ₓO₂₋ᵧ + O₂',
            reactivity=0.6,
            conditions=['very_high_voltage', 'degradation']
        ),
        '1O2': ElectrodeSpecies(
            name='单线态氧 (Singlet Oxygen)',
            smiles='O=O',  # 相同SMILES，不同电子态
            charge=0,
            role='高活性氧物种，攻击C=C双键和碳酸酯',
            origin='O⁻ + O⁻ → ¹O₂ 或 超氧根歧化',
            reactivity=1.0,
            conditions=['very_high_voltage']
        ),
        
        # --- 羟基物种 ---
        'OH': ElectrodeSpecies(
            name='羟基自由基',
            smiles='[OH]',
            charge=0,
            role='强氧化剂，H抽取反应',
            origin='O⁻ + H⁺ → ·OH 或 水分解',
            reactivity=1.0,
            conditions=['high_voltage', 'trace_water']
        ),
        'OH-': ElectrodeSpecies(
            name='氢氧根离子',
            smiles='[OH-]',
            charge=-1,
            role='亲核试剂，与碳酸酯反应',
            origin='H₂O + O²⁻ → 2OH⁻',
            reactivity=0.7,
            conditions=['trace_water']
        ),
        
        # --- 过渡金属离子 (材料溶出) ---
        'Co2+': ElectrodeSpecies(
            name='钴离子',
            smiles='[Co+2]',
            charge=2,
            role='LCO分解产物，可与溶剂配位',
            origin='LiCoO₂ + HF → Co²⁺ + LiF + H₂O',
            reactivity=0.5,
            conditions=['lco', 'nmc']
        ),
        'Co3+': ElectrodeSpecies(
            name='三价钴离子',
            smiles='[Co+3]',
            charge=3,
            role='氧化态钴，强氧化剂',
            origin='充电态LiCoO₂',
            reactivity=0.7,
            conditions=['lco', 'high_voltage']
        ),
        'Ni2+': ElectrodeSpecies(
            name='镍离子',
            smiles='[Ni+2]',
            charge=2,
            role='NMC/NCA分解产物',
            origin='LiNiO₂ + HF → Ni²⁺',
            reactivity=0.5,
            conditions=['nmc', 'nca']
        ),
        'Ni3+': ElectrodeSpecies(
            name='三价镍离子',
            smiles='[Ni+3]',
            charge=3,
            role='高氧化态镍',
            origin='充电态NMC',
            reactivity=0.65,
            conditions=['nmc', 'high_voltage']
        ),
        'Ni4+': ElectrodeSpecies(
            name='四价镍离子',
            smiles='[Ni+4]',
            charge=4,
            role='极高氧化态，不稳定',
            origin='过充电NMC',
            reactivity=0.8,
            conditions=['nmc', 'very_high_voltage']
        ),
        'Mn2+': ElectrodeSpecies(
            name='锰离子',
            smiles='[Mn+2]',
            charge=2,
            role='LMO/NMC歧化产物',
            origin='2Mn³⁺ → Mn²⁺ + Mn⁴⁺',
            reactivity=0.5,
            conditions=['lmo', 'nmc']
        ),
        'Mn3+': ElectrodeSpecies(
            name='三价锰离子',
            smiles='[Mn+3]',
            charge=3,
            role='Jahn-Teller活性，导致结构变化',
            origin='放电态LMO',
            reactivity=0.6,
            conditions=['lmo']
        ),
        'Fe2+': ElectrodeSpecies(
            name='亚铁离子',
            smiles='[Fe+2]',
            charge=2,
            role='LFP分解产物(少见)',
            origin='LiFePO₄分解',
            reactivity=0.4,
            conditions=['lfp', 'high_temp']
        ),
    }
    
    def __init__(self, 
                 electrode_type: str = 'anode',
                 cathode_material: str = 'NMC',
                 anode_material: str = 'GRAPHITE',
                 voltage: float = 3.7,
                 include_peroxide: bool = True,
                 include_superoxide: bool = True,
                 include_singlet_oxygen: bool = True,
                 include_metal_ions: bool = True):
        """
        初始化注入器
        
        Args:
            electrode_type: 电极类型 ('anode' 或 'cathode')
            cathode_material: 正极材料 (NMC, LCO, NCA, LFP, LMO等)
            anode_material: 负极材料 (LI_METAL, GRAPHITE, SILICON, NA_METAL, K_METAL, HARD_CARBON等)
            voltage: 电极电势 (V vs M/M+)
            include_peroxide: 是否包含过氧根 O₂²⁻
            include_superoxide: 是否包含超氧根 O₂⁻
            include_singlet_oxygen: 是否包含单线态氧
            include_metal_ions: 是否包含过渡金属离子
        """
        self.electrode_type = electrode_type.lower()
        self.cathode_material = cathode_material.upper()
        self.anode_material = anode_material.upper()
        self.voltage = voltage
        self.include_peroxide = include_peroxide
        self.include_superoxide = include_superoxide
        self.include_singlet_oxygen = include_singlet_oxygen
        self.include_metal_ions = include_metal_ions
        
        # 导入材料模型
        if self.electrode_type == 'cathode':
            from .cathode_materials import CATHODE_MATERIALS
            self.material_info = CATHODE_MATERIALS.get(self.cathode_material)
        else:
            from .anode_materials import ANODE_MATERIALS, CarrierIon
            self.anode_info = ANODE_MATERIALS.get(self.anode_material)
            # 确定载流子类型
            if self.anode_info:
                self.carrier_ion = self.anode_info.carrier_ion
            else:
                self.carrier_ion = CarrierIon.LITHIUM  # 默认锂
    
    def get_voltage_level(self) -> str:
        """根据电压返回电压等级"""
        if self.electrode_type == 'anode':
            if self.voltage < 0.1:
                return 'very_low_voltage'
            elif self.voltage < 1.0:
                return 'low_voltage'
            else:
                return 'normal'
        else:  # cathode
            if self.voltage > 4.5:
                return 'very_high_voltage'
            elif self.voltage > 4.0:
                return 'high_voltage'
            else:
                return 'normal'
    
    def get_species_to_inject(self) -> List[Dict[str, Any]]:
        """
        获取应该注入的物种列表
        
        Returns:
            物种信息字典列表
        """
        species_list = []
        voltage_level = self.get_voltage_level()
        
        if self.electrode_type == 'anode':
            species_list.extend(self._get_anode_species(voltage_level))
        else:
            species_list.extend(self._get_cathode_species(voltage_level))
        
        return species_list
    
    def _get_anode_species(self, voltage_level: str) -> List[Dict]:
        """获取阳极物种 - 根据材料和载流子类型"""
        result = []
        
        # 根据载流子类型选择物种库
        from .anode_materials import CarrierIon
        
        if hasattr(self, 'carrier_ion'):
            if self.carrier_ion == CarrierIon.SODIUM:
                species_dict = self.SODIUM_ANODE_SPECIES
            elif self.carrier_ion == CarrierIon.POTASSIUM:
                species_dict = self.POTASSIUM_ANODE_SPECIES
            else:
                species_dict = self.LITHIUM_ANODE_SPECIES
        else:
            species_dict = self.LITHIUM_ANODE_SPECIES
        
        # 添加载流子物种
        for key, species in species_dict.items():
            if 'always' in species.conditions:
                should_add = True
            else:
                should_add = voltage_level in species.conditions
            
            if should_add:
                result.append({
                    'key': key,
                    'name': species.name,
                    'smiles': species.smiles,
                    'charge': species.charge,
                    'role': species.role,
                    'origin': species.origin,
                    'reactivity': species.reactivity
                })
        
        # 硅负极特有物种
        if self.anode_material in ['SILICON', 'SIC']:
            for key, species in self.SILICON_ANODE_SPECIES.items():
                if 'silicon' in species.conditions:
                    if voltage_level in species.conditions or 'silicon' in species.conditions:
                        result.append({
                            'key': key,
                            'name': species.name,
                            'smiles': species.smiles,
                            'charge': species.charge,
                            'role': species.role,
                            'origin': species.origin,
                            'reactivity': species.reactivity
                        })
        
        return result

    
    def _get_cathode_species(self, voltage_level: str) -> List[Dict]:
        """获取阴极物种"""
        result = []
        material_conditions = [self.cathode_material.lower()]
        
        # 总是需要的基本物种
        always_include = ['O2-', 'O-']
        
        # 根据电压和材料添加物种
        for key, species in self.CATHODE_SPECIES.items():
            should_add = False
            
            # 检查是否是基本物种
            if key in always_include:
                should_add = True
            
            # 检查电压条件
            elif voltage_level in species.conditions:
                should_add = True
            
            # 检查材料条件
            elif any(cond in material_conditions for cond in species.conditions):
                should_add = True
            
            # 特殊条件检查
            if key == 'O2^2-' and not self.include_peroxide:
                should_add = False
            if key == 'O2-_radical' and not self.include_superoxide:
                should_add = False
            if key == '1O2' and not self.include_singlet_oxygen:
                should_add = False
            if species.charge >= 2 and 'Co' in key or 'Ni' in key or 'Mn' in key or 'Fe' in key:
                if not self.include_metal_ions:
                    should_add = False
            
            if should_add:
                result.append({
                    'key': key,
                    'name': species.name,
                    'smiles': species.smiles,
                    'charge': species.charge,
                    'role': species.role,
                    'origin': species.origin,
                    'reactivity': species.reactivity
                })
        
        return result
    
    def inject_species(self, 
                       initial_smiles: List[str],
                       verbose: bool = True) -> Tuple[List[str], List[Dict]]:
        """
        向初始分子列表注入电极特定物种
        
        Args:
            initial_smiles: 用户提供的初始分子SMILES列表
            verbose: 是否打印注入信息
            
        Returns:
            (扩展后的SMILES列表, 注入的物种信息列表)
        """
        extended_smiles = list(initial_smiles)
        injected_species = []
        
        species_to_add = self.get_species_to_inject()
        
        if verbose:
            print(f"\n{'='*60}")
            print(f"电极物种自动注入 - {self.electrode_type.upper()}")
            if self.electrode_type == 'cathode':
                print(f"正极材料: {self.cathode_material}")
            print(f"电压: {self.voltage}V vs Li/Li+")
            print(f"{'='*60}")
        
        for species in species_to_add:
            smiles = species['smiles']
            
            # 检查是否已存在
            if smiles not in extended_smiles:
                extended_smiles.append(smiles)
                injected_species.append(species)
                
                if verbose:
                    print(f"✅ 添加: {species['name']}")
                    print(f"   SMILES: {smiles}")
                    print(f"   角色: {species['role']}")
                    print(f"   来源: {species['origin']}")
                    print()
        
        if verbose:
            print(f"总计注入 {len(injected_species)} 个物种")
            print(f"扩展后分子数: {len(extended_smiles)}")
            print(f"{'='*60}\n")
        
        return extended_smiles, injected_species


def auto_inject_electrode_species(
    initial_smiles: List[str],
    electrode_type: str = 'anode',
    cathode_material: str = 'NMC',
    voltage: float = 3.7,
    include_peroxide: bool = True,
    include_superoxide: bool = True,
    verbose: bool = True
) -> Tuple[List[str], List[Dict]]:
    """
    自动注入电极物种 (便捷函数)
    
    Args:
        initial_smiles: 初始分子SMILES列表
        electrode_type: 电极类型
        cathode_material: 正极材料
        voltage: 电压
        include_peroxide: 包含过氧根
        include_superoxide: 包含超氧根
        verbose: 打印信息
        
    Returns:
        (扩展的SMILES列表, 注入的物种信息)
    """
    injector = ElectrodeSpeciesInjector(
        electrode_type=electrode_type,
        cathode_material=cathode_material,
        voltage=voltage,
        include_peroxide=include_peroxide,
        include_superoxide=include_superoxide
    )
    
    return injector.inject_species(initial_smiles, verbose=verbose)
