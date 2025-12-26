"""
CEI (Cathode Electrolyte Interphase) Component Library for RSNet
CEI组分库 - 阴极/正极界面层化学组分

CEI与SEI的区别:
- SEI (阳极): 主要由还原反应产生，含Li₂CO₃, LiF, 有机锂盐
- CEI (阴极): 主要由氧化反应产生，含过渡金属氧化物、聚碳酸酯、氧化产物

CEI主要组分:
1. 无机组分: Li₂CO₃, LiF, Li₃PO₄, MxOy
2. 有机组分: 聚碳酸酯, 羧酸锂, 氧化聚合物
3. 过渡金属配合物: M-溶剂配合物

作者: RSNet Team
版本: 2.0
"""

from typing import Dict, List, Optional
from dataclasses import dataclass, field
from enum import Enum


class CEIComponentType(Enum):
    """CEI组分类型"""
    INORGANIC_LITHIUM = "inorganic_li"      # 无机锂盐
    INORGANIC_METAL_OXIDE = "metal_oxide"    # 过渡金属氧化物
    ORGANIC_LITHIUM = "organic_li"           # 有机锂盐
    POLYMER = "polymer"                       # 聚合物
    OLIGOMER = "oligomer"                     # 低聚物
    SMALL_MOLECULE = "small_molecule"         # 小分子
    METAL_COMPLEX = "metal_complex"           # 金属配合物
    PEROXIDE = "peroxide"                     # 过氧化物


class InterphaseFace(Enum):
    """界面层类型"""
    SEI = "sei"            # 阳极界面 (Solid Electrolyte Interphase)
    CEI = "cei"            # 阴极界面 (Cathode Electrolyte Interphase)
    BOTH = "both"          # 两者都有


@dataclass
class InterfaceComponent:
    """界面层组分"""
    name: str                               # 组分名称
    formula: str                            # 化学式
    smiles_patterns: List[str]              # SMILES模式匹配
    component_type: CEIComponentType         # 组分类型
    interphase: InterphaseFace              # 所属界面
    importance: float                        # 重要性评分 (0-1)
    origin: str                             # 形成来源
    description: str                        # 描述
    parent_reactions: List[str] = field(default_factory=list)  # 生成反应


# =============================================================================
# SEI组分 (阳极界面层)
# =============================================================================

SEI_COMPONENTS: List[InterfaceComponent] = [
    InterfaceComponent(
        name='Li₂CO₃',
        formula='Li2CO3',
        smiles_patterns=['[Li]OC(=O)O[Li]', 'O=C([O-])[O-].[Li+].[Li+]'],
        component_type=CEIComponentType.INORGANIC_LITHIUM,
        interphase=InterphaseFace.BOTH,
        importance=1.0,
        origin='EC/PC还原分解',
        description='主要无机SEI组分，提供Li⁺传导性',
        parent_reactions=['EC + 2e⁻ + 2Li⁺ → Li₂CO₃ + C₂H₄']
    ),
    InterfaceComponent(
        name='LiF',
        formula='LiF',
        smiles_patterns=['[Li]F', '[Li+].[F-]'],
        component_type=CEIComponentType.INORGANIC_LITHIUM,
        interphase=InterphaseFace.BOTH,
        importance=1.0,
        origin='LiPF₆分解 / HF与Li₂O反应',
        description='关键SEI/CEI组分，电子绝缘但Li⁺导通',
        parent_reactions=['LiPF₆ → LiF + PF₅', 'Li₂O + 2HF → 2LiF + H₂O']
    ),
    InterfaceComponent(
        name='Li₂O',
        formula='Li2O',
        smiles_patterns=['[Li]O[Li]', '[Li+].[O-2].[Li+]'],
        component_type=CEIComponentType.INORGANIC_LITHIUM,
        interphase=InterphaseFace.SEI,
        importance=0.8,
        origin='Li与微量O₂反应 / 碳酸锂分解',
        description='内层SEI组分',
        parent_reactions=['4Li + O₂ → 2Li₂O']
    ),
    InterfaceComponent(
        name='LEDC (Li₂EDC)',
        formula='(CH2OCO2Li)2',
        smiles_patterns=['[Li]OC(=O)OCC', 'CCOC(=O)O[Li]', '[Li]OC(=O)OCCOC(=O)O[Li]'],
        component_type=CEIComponentType.ORGANIC_LITHIUM,
        interphase=InterphaseFace.SEI,
        importance=0.9,
        origin='EC单电子还原',
        description='锂乙基二碳酸酯，主要有机SEI组分',
        parent_reactions=['2EC + 2e⁻ + 2Li⁺ → LEDC + C₂H₄']
    ),
    InterfaceComponent(
        name='ROLi (烷氧锂)',
        formula='ROLi',
        smiles_patterns=['[Li]OC', '[Li]OCC', 'C[O-].[Li+]'],
        component_type=CEIComponentType.ORGANIC_LITHIUM,
        interphase=InterphaseFace.SEI,
        importance=0.7,
        origin='线性碳酸酯还原',
        description='烷氧锂盐',
        parent_reactions=['DMC + e⁻ + Li⁺ → CH₃OLi + CH₃O·']
    ),
    InterfaceComponent(
        name='LiₓPFᵧOᵤ',
        formula='LixPFyOz',
        smiles_patterns=['[Li]OP', 'FP(F)O[Li]', 'O=P([O-])(F)F'],
        component_type=CEIComponentType.INORGANIC_LITHIUM,
        interphase=InterphaseFace.BOTH,
        importance=0.6,
        origin='LiPF₆水解',
        description='含磷氟化物',
        parent_reactions=['LiPF₆ + H₂O → LiPOF₂ + 2HF']
    ),
]

# =============================================================================
# CEI组分 (阴极界面层) - 氧化反应产物
# =============================================================================

CEI_COMPONENTS: List[InterfaceComponent] = [
    InterfaceComponent(
        name='Li₂CO₃ (CEI)',
        formula='Li2CO3',
        smiles_patterns=['[Li]OC(=O)O[Li]'],
        component_type=CEIComponentType.INORGANIC_LITHIUM,
        interphase=InterphaseFace.CEI,
        importance=1.0,
        origin='碳酸酯氧化分解',
        description='阴极表面主要无机组分',
        parent_reactions=['EC + O⁻ → Li₂CO₃ + oxidized organics']
    ),
    InterfaceComponent(
        name='聚碳酸酯 (Polycarbonates)',
        formula='(OC(=O)O)n',
        smiles_patterns=['OC(=O)OC(=O)O', 'C(=O)OC(=O)', 'OC(=O)OC'],
        component_type=CEIComponentType.POLYMER,
        interphase=InterphaseFace.CEI,
        importance=0.85,
        origin='EC/PC氧化聚合',
        description='碳酸酯开环聚合产物',
        parent_reactions=['nEC + O⁻ → 聚碳酸酯']
    ),
    InterfaceComponent(
        name='Li₂O₂ (过氧化锂)',
        formula='Li2O2',
        smiles_patterns=['[Li]OO[Li]', '[Li+].[O-][O-].[Li+]'],
        component_type=CEIComponentType.PEROXIDE,
        interphase=InterphaseFace.CEI,
        importance=0.8,
        origin='Li与O₂反应 (高电压)',
        description='过氧化物中间体，可进一步分解',
        parent_reactions=['2Li + O₂ → Li₂O₂']
    ),
    InterfaceComponent(
        name='LiO₂ (超氧化锂)',
        formula='LiO2',
        smiles_patterns=['[Li]O[O]', '[Li+].[O-][O]'],
        component_type=CEIComponentType.PEROXIDE,
        interphase=InterphaseFace.CEI,
        importance=0.75,
        origin='Li + O₂⁻ 反应',
        description='超氧化物，高反应性',
        parent_reactions=['Li + O₂⁻ → LiO₂']
    ),
    InterfaceComponent(
        name='羧酸锂 (RCO₂Li)',
        formula='RCO2Li',
        smiles_patterns=['C(=O)O[Li]', 'CC(=O)O[Li]', '[Li]OC(=O)C'],
        component_type=CEIComponentType.ORGANIC_LITHIUM,
        interphase=InterphaseFace.CEI,
        importance=0.7,
        origin='溶剂C-H氧化',
        description='碳酸酯氧化降解产物',
        parent_reactions=['R-H + O⁻ → R-OH → R-CO₂Li']
    ),
    InterfaceComponent(
        name='醛类 (RCHO)',
        formula='RCHO',
        smiles_patterns=['C=O', 'CC=O', 'OCC=O'],
        component_type=CEIComponentType.SMALL_MOLECULE,
        interphase=InterphaseFace.CEI,
        importance=0.5,
        origin='醇氧化',
        description='醇类氧化中间体',
        parent_reactions=['R-CH₂OH + O⁻ → R-CHO + OH⁻']
    ),
    InterfaceComponent(
        name='CO₂',
        formula='CO2',
        smiles_patterns=['O=C=O', 'C(=O)=O'],
        component_type=CEIComponentType.SMALL_MOLECULE,
        interphase=InterphaseFace.BOTH,
        importance=0.6,
        origin='碳酸酯完全分解',
        description='碳酸酯氧化最终产物',
        parent_reactions=['Li₂CO₃ → Li₂O + CO₂']
    ),
    InterfaceComponent(
        name='过渡金属碳酸盐',
        formula='MCO3',
        smiles_patterns=['[Co]OC(=O)O', '[Ni]OC(=O)O', '[Mn]OC(=O)O'],
        component_type=CEIComponentType.INORGANIC_METAL_OXIDE,
        interphase=InterphaseFace.CEI,
        importance=0.65,
        origin='溶出金属与碳酸根反应',
        description='过渡金属碳酸盐沉淀',
        parent_reactions=['M²⁺ + CO₃²⁻ → MCO₃']
    ),
    InterfaceComponent(
        name='过渡金属氟化物',
        formula='MF2',
        smiles_patterns=['[Co]F', '[Ni]F', '[Mn]F', 'F[Co]F', 'F[Ni]F'],
        component_type=CEIComponentType.INORGANIC_METAL_OXIDE,
        interphase=InterphaseFace.CEI,
        importance=0.6,
        origin='溶出金属与F⁻反应',
        description='过渡金属氟化物',
        parent_reactions=['M²⁺ + 2F⁻ → MF₂']
    ),
    InterfaceComponent(
        name='过渡金属-溶剂配合物',
        formula='M(solvent)n',
        smiles_patterns=['[Co]OC', '[Ni]OC', '[Mn]OC'],
        component_type=CEIComponentType.METAL_COMPLEX,
        interphase=InterphaseFace.CEI,
        importance=0.5,
        origin='溶出金属与溶剂配位',
        description='过渡金属配位化合物',
        parent_reactions=['M²⁺ + nSolvent → M(Solvent)ₙ²⁺']
    ),
    InterfaceComponent(
        name='Li₃PO₄',
        formula='Li3PO4',
        smiles_patterns=['[Li]OP(=O)([O-])[O-]', 'O=P([O-])([O-])[O-].[Li+].[Li+].[Li+]'],
        component_type=CEIComponentType.INORGANIC_LITHIUM,
        interphase=InterphaseFace.BOTH,
        importance=0.55,
        origin='LiPF₆完全水解',
        description='磷酸锂，稳定的界面组分',
        parent_reactions=['LiPF₆ + 4H₂O → Li₃PO₄ + 6HF + LiOH']
    ),
]

# =============================================================================
# 合并组分库
# =============================================================================

ALL_INTERFACE_COMPONENTS: List[InterfaceComponent] = SEI_COMPONENTS + CEI_COMPONENTS


class InterfaceComponentRecognizer:
    """界面组分识别器"""
    
    def __init__(self, interphase_type: InterphaseFace = InterphaseFace.BOTH):
        """
        初始化识别器
        
        Args:
            interphase_type: 要识别的界面类型
        """
        self.interphase_type = interphase_type
        self.components = self._filter_components()
    
    def _filter_components(self) -> List[InterfaceComponent]:
        """根据界面类型筛选组分"""
        if self.interphase_type == InterphaseFace.BOTH:
            return ALL_INTERFACE_COMPONENTS
        return [c for c in ALL_INTERFACE_COMPONENTS 
                if c.interphase == self.interphase_type or c.interphase == InterphaseFace.BOTH]
    
    def recognize(self, smiles: str) -> List[InterfaceComponent]:
        """
        识别SMILES对应的界面组分
        
        Args:
            smiles: 分子SMILES
            
        Returns:
            匹配的组分列表
        """
        matches = []
        for component in self.components:
            for pattern in component.smiles_patterns:
                if pattern in smiles or smiles in pattern:
                    matches.append(component)
                    break
        return matches
    
    def is_interface_component(self, smiles: str) -> bool:
        """检查是否是界面层组分"""
        return len(self.recognize(smiles)) > 0
    
    def get_importance_score(self, smiles: str) -> float:
        """获取组分重要性评分"""
        components = self.recognize(smiles)
        if not components:
            return 0.0
        return max(c.importance for c in components)
    
    def classify_products(self, products: List[str]) -> Dict[str, List[str]]:
        """
        将产物分类
        
        Args:
            products: 产物SMILES列表
            
        Returns:
            {组分类型: [SMILES列表]}
        """
        classification = {t.value: [] for t in CEIComponentType}
        
        for smiles in products:
            components = self.recognize(smiles)
            for comp in components:
                classification[comp.component_type.value].append(smiles)
        
        return {k: list(set(v)) for k, v in classification.items() if v}


# 便捷函数
def recognize_sei_component(smiles: str) -> List[InterfaceComponent]:
    """识别SEI组分"""
    recognizer = InterfaceComponentRecognizer(InterphaseFace.SEI)
    return recognizer.recognize(smiles)


def recognize_cei_component(smiles: str) -> List[InterfaceComponent]:
    """识别CEI组分"""
    recognizer = InterfaceComponentRecognizer(InterphaseFace.CEI)
    return recognizer.recognize(smiles)


def get_all_sei_patterns() -> List[str]:
    """获取所有SEI组分的SMILES模式"""
    patterns = []
    for comp in SEI_COMPONENTS:
        patterns.extend(comp.smiles_patterns)
    return list(set(patterns))


def get_all_cei_patterns() -> List[str]:
    """获取所有CEI组分的SMILES模式"""
    patterns = []
    for comp in CEI_COMPONENTS:
        patterns.extend(comp.smiles_patterns)
    return list(set(patterns))
