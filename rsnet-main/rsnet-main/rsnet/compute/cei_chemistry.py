"""
CEI Chemistry Integration Module for RSNet
CEI化学集成模块 - 统一入口点

本模块提供阴极界面化学的完整支持:
1. 电极物种自动注入 (过氧根、超氧根、单线态氧等)
2. 正极材料特性 (LCO, NMC, NCA, LFP等)
3. CEI组分识别
4. 氧化反应算符

使用方式:
---------
from rsnet.compute.cei_chemistry import integrate_cei_chemistry

# 扩展初始分子列表
extended_smiles, injected_info = integrate_cei_chemistry(
    initial_smiles=['C1COC(=O)O1', '[Li+]'],
    electrode_type='cathode',
    cathode_material='NMC',
    voltage=4.5
)

作者: RSNet Team
版本: 2.0
"""

from typing import Dict, List, Optional, Tuple, Any

# 导入子模块
from .electrode_species_injector import (
    ElectrodeSpeciesInjector, 
    ElectrodeType,
    ElectrodeSpecies,
    auto_inject_electrode_species
)
from .cathode_materials import (
    CathodeMaterial,
    CathodeMaterialFamily,
    OxygenReleaseProfile,
    OxygenReleaseType,
    CATHODE_MATERIALS,
    get_cathode_material,
    list_available_materials,
    print_material_summary
)
from .interface_components import (
    InterfaceComponent,
    CEIComponentType,
    InterphaseFace,
    InterfaceComponentRecognizer,
    SEI_COMPONENTS,
    CEI_COMPONENTS,
    ALL_INTERFACE_COMPONENTS,
    recognize_sei_component,
    recognize_cei_component
)
from .oxidation_operators import (
    OxidationOperator,
    OxidationMechanism,
    OxidationOperatorManager,
    OXIDATION_OPERATORS,
    get_oxidation_operators_for_cathode,
    estimate_oxidation_barrier
)


class CEIChemistryIntegrator:
    """
    CEI化学集成器
    
    一站式解决方案，用于将CEI化学整合到反应网络生成中
    """
    
    def __init__(self,
                 electrode_type: str = 'anode',
                 cathode_material: str = 'NMC',
                 voltage: float = 3.7,
                 include_peroxide: bool = True,
                 include_superoxide: bool = True,
                 include_singlet_oxygen: bool = True,
                 include_metal_ions: bool = True,
                 verbose: bool = True):
        """
        初始化集成器
        
        Args:
            electrode_type: 电极类型 ('anode' 或 'cathode')
            cathode_material: 正极材料 (NMC, LCO, NCA, LFP, LMO等)
            voltage: 电极电势 (V vs Li/Li+)
            include_peroxide: 包含过氧根O₂²⁻
            include_superoxide: 包含超氧根O₂⁻
            include_singlet_oxygen: 包含单线态氧¹O₂
            include_metal_ions: 包含过渡金属离子
            verbose: 打印详细信息
        """
        self.electrode_type = electrode_type.lower()
        self.cathode_material = cathode_material.upper()
        self.voltage = voltage
        self.include_peroxide = include_peroxide
        self.include_superoxide = include_superoxide
        self.include_singlet_oxygen = include_singlet_oxygen
        self.include_metal_ions = include_metal_ions
        self.verbose = verbose
        
        # 初始化子模块
        self.species_injector = ElectrodeSpeciesInjector(
            electrode_type=electrode_type,
            cathode_material=cathode_material,
            voltage=voltage,
            include_peroxide=include_peroxide,
            include_superoxide=include_superoxide,
            include_singlet_oxygen=include_singlet_oxygen,
            include_metal_ions=include_metal_ions
        )
        
        self.component_recognizer = InterfaceComponentRecognizer(
            interphase_type=InterphaseFace.CEI if electrode_type == 'cathode' else InterphaseFace.SEI
        )
        
        self.oxidation_manager = OxidationOperatorManager(
            include_peroxide=include_peroxide,
            include_superoxide=include_superoxide,
            include_singlet_oxygen=include_singlet_oxygen,
            include_metal_catalysis=include_metal_ions
        )
        
        # 获取材料信息
        self.material_info = get_cathode_material(cathode_material)
    
    def integrate(self, initial_smiles: List[str]) -> Dict[str, Any]:
        """
        执行完整的CEI化学集成
        
        Args:
            initial_smiles: 初始分子SMILES列表
            
        Returns:
            包含扩展分子和相关信息的字典
        """
        result = {
            'original_smiles': initial_smiles,
            'extended_smiles': [],
            'injected_species': [],
            'electrode_type': self.electrode_type,
            'cathode_material': self.cathode_material,
            'voltage': self.voltage,
            'material_info': None,
            'active_oxidation_operators': [],
            'warnings': []
        }
        
        # 1. 注入电极特定物种
        extended_smiles, injected = self.species_injector.inject_species(
            initial_smiles, 
            verbose=self.verbose
        )
        result['extended_smiles'] = extended_smiles
        result['injected_species'] = injected
        
        # 2. 获取材料信息
        if self.material_info:
            result['material_info'] = {
                'code': self.material_info.code,
                'name': self.material_info.name,
                'formula': self.material_info.formula,
                'voltage_range': (self.material_info.voltage_min, self.material_info.voltage_max),
                'is_oxygen_releasing': self.material_info.is_oxygen_releasing(self.voltage),
                'stability_score': self.material_info.get_stability_score(self.voltage)
            }
            
            # 检查电压警告
            if self.voltage > self.material_info.voltage_max:
                result['warnings'].append(
                    f"⚠️ 电压 {self.voltage}V 超过 {self.cathode_material} 最高推荐电压 {self.material_info.voltage_max}V"
                )
            
            if self.material_info.oxygen_release:
                if self.voltage >= self.material_info.oxygen_release.onset_voltage:
                    result['warnings'].append(
                        f"⚠️ 电压已进入氧释放区域 (onset: {self.material_info.oxygen_release.onset_voltage}V)"
                    )
        
        # 3. 获取可用的氧化算符
        available_oxidants = [s['smiles'] for s in injected]
        result['active_oxidation_operators'] = []
        
        for smiles in initial_smiles:
            applicable = self.oxidation_manager.get_applicable_operators(smiles, available_oxidants)
            for op in applicable:
                result['active_oxidation_operators'].append({
                    'operator_name': op.name,
                    'mechanism': op.mechanism.value,
                    'target_molecule': smiles,
                    'activation_energy_range': op.activation_energy_range
                })
        
        if self.verbose:
            print(f"\n🔬 CEI化学集成完成")
            print(f"   原始分子: {len(initial_smiles)}个")
            print(f"   扩展后: {len(extended_smiles)}个")
            print(f"   注入物种: {len(injected)}个")
            print(f"   激活算符: {len(result['active_oxidation_operators'])}个")
            if result['warnings']:
                for w in result['warnings']:
                    print(f"   {w}")
        
        return result
    
    def score_product_as_interface_component(self, product_smiles: str) -> Dict[str, Any]:
        """
        评估产物作为界面层组分的可能性
        
        Args:
            product_smiles: 产物SMILES
            
        Returns:
            组分评分信息
        """
        components = self.component_recognizer.recognize(product_smiles)
        
        if components:
            best = max(components, key=lambda c: c.importance)
            return {
                'is_interface_component': True,
                'component_name': best.name,
                'component_type': best.component_type.value,
                'importance': best.importance,
                'origin': best.origin,
                'all_matches': [c.name for c in components]
            }
        
        return {
            'is_interface_component': False,
            'component_name': None,
            'importance': 0.0
        }


def integrate_cei_chemistry(
    initial_smiles: List[str],
    electrode_type: str = 'anode',
    cathode_material: str = 'NMC',
    voltage: float = 3.7,
    include_peroxide: bool = True,
    include_superoxide: bool = True,
    verbose: bool = True
) -> Tuple[List[str], Dict[str, Any]]:
    """
    便捷函数: 将CEI化学集成到反应网络
    
    Args:
        initial_smiles: 初始分子SMILES列表
        electrode_type: 电极类型
        cathode_material: 正极材料
        voltage: 电压
        include_peroxide: 包含过氧根
        include_superoxide: 包含超氧根
        verbose: 打印信息
        
    Returns:
        (扩展的SMILES列表, 集成信息字典)
    """
    integrator = CEIChemistryIntegrator(
        electrode_type=electrode_type,
        cathode_material=cathode_material,
        voltage=voltage,
        include_peroxide=include_peroxide,
        include_superoxide=include_superoxide,
        verbose=verbose
    )
    
    result = integrator.integrate(initial_smiles)
    
    return result['extended_smiles'], result


# 导出
__all__ = [
    # 主要类
    'CEIChemistryIntegrator',
    'ElectrodeSpeciesInjector',
    'InterfaceComponentRecognizer',
    'OxidationOperatorManager',
    
    # 材料
    'CathodeMaterial',
    'CATHODE_MATERIALS',
    'get_cathode_material',
    'list_available_materials',
    
    # 组分
    'InterfaceComponent',
    'SEI_COMPONENTS',
    'CEI_COMPONENTS',
    
    # 算符
    'OXIDATION_OPERATORS',
    'OxidationOperator',
    
    # 便捷函数
    'integrate_cei_chemistry',
    'auto_inject_electrode_species',
    'recognize_sei_component',
    'recognize_cei_component',
]
