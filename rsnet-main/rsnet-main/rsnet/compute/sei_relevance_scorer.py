"""
SEI Relevance Scorer and Reaction Prioritization for RSNet
SEI相关反应评分器 - 聚焦电池界面化学相关反应

特点:
1. 基于化学直觉的相关性评分
2. 考虑Li参与、碳酸酯反应、电极匹配
3. 支持已知SEI组分识别
4. 提供综合优先级排序

作者: RSNet Team
版本: 2.0 (改进版)
"""

from typing import Dict, List, Optional, Tuple, Set
from dataclasses import dataclass, field
from enum import Enum
import re


class SEIComponentType(Enum):
    """SEI组分类型"""
    INORGANIC_LITHIUM_SALT = "inorganic_li_salt"    # Li2CO3, LiF, Li2O
    ORGANIC_LITHIUM_SALT = "organic_li_salt"         # ROCO2Li, ROLi
    POLYMER = "polymer"                               # 聚合物
    OLIGOMER = "oligomer"                             # 低聚物
    SMALL_MOLECULE = "small_molecule"                 # 小分子产物


@dataclass
class SEIComponent:
    """已知SEI组分"""
    name: str
    smiles_patterns: List[str]
    component_type: SEIComponentType
    importance: float  # 0-1, 重要性评分
    description: str


class SEIRelevanceScorer:
    """
    SEI相关反应评分器
    
    评分维度:
    1. 锂参与 (30%)
    2. 碳酸酯/电解液反应 (20%)
    3. 电极类型匹配 (20%)
    4. 产物稳定性 (15%)
    5. 已知SEI组分形成 (15%)
    """
    
    # 已知SEI组分
    KNOWN_SEI_COMPONENTS: List[SEIComponent] = [
        SEIComponent(
            name="Li2CO3",
            smiles_patterns=["[Li]OC(=O)O[Li]", "O=C([O-])[O-]"],
            component_type=SEIComponentType.INORGANIC_LITHIUM_SALT,
            importance=1.0,
            description="碳酸锂 - 主要无机SEI组分"
        ),
        SEIComponent(
            name="LiF",
            smiles_patterns=["[Li]F", "[Li+].[F-]"],
            component_type=SEIComponentType.INORGANIC_LITHIUM_SALT,
            importance=1.0,
            description="氟化锂 - 关键SEI组分，来自LiPF6分解"
        ),
        SEIComponent(
            name="Li2O",
            smiles_patterns=["[Li]O[Li]", "[Li+].[O-2].[Li+]"],
            component_type=SEIComponentType.INORGANIC_LITHIUM_SALT,
            importance=0.8,
            description="氧化锂 - 内层SEI组分"
        ),
        SEIComponent(
            name="LiOH",
            smiles_patterns=["[Li]O", "[Li+].[OH-]"],
            component_type=SEIComponentType.INORGANIC_LITHIUM_SALT,
            importance=0.7,
            description="氢氧化锂"
        ),
        SEIComponent(
            name="LEDC",
            smiles_patterns=["[Li]OC(=O)OCC", "CCOC(=O)O[Li]"],
            component_type=SEIComponentType.ORGANIC_LITHIUM_SALT,
            importance=0.9,
            description="锂乙基二碳酸酯 - EC分解产物"
        ),
        SEIComponent(
            name="LDMC",
            smiles_patterns=["COC(=O)O[Li]", "[Li]OC(=O)OC"],
            component_type=SEIComponentType.ORGANIC_LITHIUM_SALT,
            importance=0.8,
            description="锂二甲基碳酸酯 - DMC分解产物"
        ),
        SEIComponent(
            name="Li ethoxide",
            smiles_patterns=["CCO[Li]", "[Li]OCC"],
            component_type=SEIComponentType.ORGANIC_LITHIUM_SALT,
            importance=0.7,
            description="乙醇锂"
        ),
        SEIComponent(
            name="Li methoxide",
            smiles_patterns=["CO[Li]", "[Li]OC"],
            component_type=SEIComponentType.ORGANIC_LITHIUM_SALT,
            importance=0.6,
            description="甲醇锂"
        ),
        SEIComponent(
            name="PF5",
            smiles_patterns=["FP(F)(F)(F)F", "F[P](F)(F)(F)F"],
            component_type=SEIComponentType.SMALL_MOLECULE,
            importance=0.5,
            description="五氟化磷 - LiPF6分解中间体"
        ),
        SEIComponent(
            name="POF3",
            smiles_patterns=["FP(=O)(F)F", "O=P(F)(F)F"],
            component_type=SEIComponentType.SMALL_MOLECULE,
            importance=0.5,
            description="三氟氧化磷 - 水解产物"
        ),
        SEIComponent(
            name="CO2",
            smiles_patterns=["O=C=O", "C(=O)=O"],
            component_type=SEIComponentType.SMALL_MOLECULE,
            importance=0.4,
            description="二氧化碳 - 碳酸酯分解产物"
        ),
        SEIComponent(
            name="C2H4",
            smiles_patterns=["C=C", "[CH2]=[CH2]"],
            component_type=SEIComponentType.SMALL_MOLECULE,
            importance=0.4,
            description="乙烯 - EC开环产物"
        ),
    ]
    
    # 碳酸酯电解液SMILES模式
    CARBONATE_PATTERNS = [
        r'C1COC\(=O\)O1',      # EC (环状碳酸酯)
        r'COC\(=O\)OC',         # DMC
        r'CCOC\(=O\)OCC',       # DEC
        r'COC\(=O\)OCC',        # EMC
        r'C1CCOC\(=O\)O1',      # PC
        r'O=C\(O',              # 通用碳酸酯
    ]
    
    # 电解液盐SMILES模式
    SALT_PATTERNS = [
        r'\[Li\+\]',            # Li+
        r'F\[P-\]\(F\)\(F\)\(F\)\(F\)F',  # PF6-
        r'FC\(F\)\(F\)S',       # CF3- (TFSI相关)
        r'\[BF4-\]',            # BF4-
    ]
    
    def __init__(self, electrode_type: str = 'anode'):
        """
        初始化评分器
        
        Args:
            electrode_type: 电极类型 ('anode' 或 'cathode')
        """
        self.electrode_type = electrode_type
        
        # 编译正则表达式
        self._carbonate_re = [re.compile(p) for p in self.CARBONATE_PATTERNS]
        self._salt_re = [re.compile(p) for p in self.SALT_PATTERNS]
    
    def score_reaction(self, reaction) -> Tuple[float, Dict[str, float]]:
        """
        评分反应对SEI形成的相关性
        
        Args:
            reaction: 反应对象 (需要有reactants和products属性)
            
        Returns:
            (总分0.0-1.0, 各维度分数字典)
        """
        scores = {}
        
        # 1. Li参与 (30%)
        scores['lithium_involvement'] = self._score_lithium_involvement(reaction)
        
        # 2. 碳酸酯/电解液反应 (20%)
        scores['electrolyte_relevance'] = self._score_electrolyte_relevance(reaction)
        
        # 3. 电极类型匹配 (20%)
        scores['electrode_match'] = self._score_electrode_match(reaction)
        
        # 4. 产物稳定性 (15%)
        scores['product_stability'] = self._score_product_stability(reaction)
        
        # 5. 已知SEI组分形成 (15%)
        scores['sei_component_formation'] = self._score_sei_component_formation(reaction)
        
        # 加权总分
        total_score = (
            0.30 * scores['lithium_involvement'] +
            0.20 * scores['electrolyte_relevance'] +
            0.20 * scores['electrode_match'] +
            0.15 * scores['product_stability'] +
            0.15 * scores['sei_component_formation']
        )
        
        return min(1.0, total_score), scores
    
    def _score_lithium_involvement(self, reaction) -> float:
        """评估锂参与程度"""
        score = 0.0
        
        # 检查反应物和产物中的Li
        reactant_li = self._count_lithium(reaction.reactants)
        product_li = self._count_lithium(reaction.products)
        
        if reactant_li > 0 or product_li > 0:
            score += 0.5  # Li参与基础分
            
            # Li作为反应物（如Li+还原）
            if reactant_li > 0:
                score += 0.25
            
            # Li作为产物（如形成Li盐）
            if product_li > 0:
                score += 0.25
        
        return min(1.0, score)
    
    def _score_electrolyte_relevance(self, reaction) -> float:
        """评估电解液相关性"""
        score = 0.0
        
        # 检查碳酸酯参与
        for mol in reaction.reactants:
            smiles = getattr(mol, 'smiles', '')
            for pattern in self._carbonate_re:
                if pattern.search(smiles):
                    score += 0.5
                    break
        
        # 检查电解液盐参与
        for mol in reaction.reactants:
            smiles = getattr(mol, 'smiles', '')
            for pattern in self._salt_re:
                if pattern.search(smiles):
                    score += 0.5
                    break
        
        return min(1.0, score)
    
    def _score_electrode_match(self, reaction) -> float:
        """评估反应类型与电极匹配程度"""
        reaction_type = getattr(reaction, 'reaction_type', None)
        reaction_name = getattr(reaction, 'name', '').lower()
        
        anode_reactions = {
            'reduction', 'electron_gain', 'ring_opening', 
            'li_insertion', 'reductive_cleavage'
        }
        cathode_reactions = {
            'oxidation', 'electron_loss', 'dehydrogenation',
            'li_extraction', 'oxidative_cleavage'
        }
        
        if self.electrode_type == 'anode':
            # 阳极主要发生还原反应
            if reaction_type in anode_reactions:
                return 1.0
            for keyword in ['reduct', 'ring_open', 'electron']:
                if keyword in reaction_name:
                    return 0.8
            return 0.3  # 未明确的反应给予低分
        else:
            # 阴极主要发生氧化反应
            if reaction_type in cathode_reactions:
                return 1.0
            for keyword in ['oxid', 'dehydrogen']:
                if keyword in reaction_name:
                    return 0.8
            return 0.3
    
    def _score_product_stability(self, reaction) -> float:
        """评估产物稳定性"""
        score = 0.0
        total_products = len(reaction.products)
        
        for product in reaction.products:
            # 检查是否是已知稳定产物
            smiles = getattr(product, 'smiles', '')
            
            # 小分子稳定产物
            small_stable = ['O=C=O', 'C=C', '[Li]F', '[Li]O[Li]']
            for stable in small_stable:
                if stable in smiles:
                    score += 1.0 / total_products
                    break
            
            # 无机锂盐 (高稳定性)
            if re.search(r'\[Li\]O|Li[OF]', smiles):
                score += 0.5 / total_products
        
        return min(1.0, score)
    
    def _score_sei_component_formation(self, reaction) -> float:
        """评估是否形成已知SEI组分"""
        max_importance = 0.0
        
        for product in reaction.products:
            smiles = getattr(product, 'smiles', '')
            
            for component in self.KNOWN_SEI_COMPONENTS:
                for pattern in component.smiles_patterns:
                    if pattern in smiles or smiles in pattern:
                        max_importance = max(max_importance, component.importance)
                        break
        
        return max_importance
    
    def _count_lithium(self, molecules) -> int:
        """计算分子列表中锂原子数"""
        count = 0
        for mol in molecules:
            smiles = getattr(mol, 'smiles', '')
            count += smiles.count('Li') + smiles.count('[Li+]')
        return count


class ReactionPrioritizer:
    """
    反应优先级排序器
    
    综合考虑:
    1. SEI相关性 (40%)
    2. 热力学可行性 (30%)
    3. 动力学可及性 (30%)
    """
    
    def __init__(self, 
                 sei_scorer: SEIRelevanceScorer,
                 thermodynamic_cutoff: float = 40.0,
                 kinetic_cutoff: float = 35.0):
        """
        初始化排序器
        
        Args:
            sei_scorer: SEI相关性评分器
            thermodynamic_cutoff: 热力学截断值 (kcal/mol)
            kinetic_cutoff: 动力学截断值 (kcal/mol)
        """
        self.sei_scorer = sei_scorer
        self.thermodynamic_cutoff = thermodynamic_cutoff
        self.kinetic_cutoff = kinetic_cutoff
    
    def prioritize_reactions(self, reactions) -> List[Tuple[any, float, Dict[str, float]]]:
        """
        按优先级排序反应
        
        Args:
            reactions: 反应列表
            
        Returns:
            [(反应, 总分, 详细分数)] 按分数降序排列
        """
        scored_reactions = []
        
        for reaction in reactions:
            # 获取SEI相关性分数
            relevance_score, relevance_details = self.sei_scorer.score_reaction(reaction)
            
            # 热力学得分
            delta_E = abs(getattr(reaction, 'reaction_energy', 0) or 0)
            thermo_score = max(0, 1 - delta_E / self.thermodynamic_cutoff)
            
            # 动力学得分
            Ea = getattr(reaction, 'activation_energy', 30.0) or 30.0
            kinetic_score = max(0, 1 - Ea / self.kinetic_cutoff)
            
            # 综合得分
            total_score = (
                0.40 * relevance_score +
                0.30 * thermo_score +
                0.30 * kinetic_score
            )
            
            details = {
                **relevance_details,
                'relevance_total': relevance_score,
                'thermodynamic_score': thermo_score,
                'kinetic_score': kinetic_score,
                'delta_E': delta_E,
                'Ea': Ea,
                'total_score': total_score
            }
            
            scored_reactions.append((reaction, total_score, details))
        
        # 按分数降序排列
        scored_reactions.sort(key=lambda x: x[1], reverse=True)
        
        return scored_reactions
    
    def filter_and_prioritize(self, 
                              reactions,
                              min_relevance: float = 0.3,
                              min_total_score: float = 0.2) -> List:
        """
        筛选并排序反应
        
        Args:
            reactions: 反应列表
            min_relevance: 最小相关性阈值
            min_total_score: 最小总分阈值
            
        Returns:
            筛选后的反应列表
        """
        prioritized = self.prioritize_reactions(reactions)
        
        filtered = [
            reaction for reaction, score, details in prioritized
            if details['relevance_total'] >= min_relevance and score >= min_total_score
        ]
        
        return filtered
    
    def get_top_reactions(self, reactions, top_n: int = 20) -> List:
        """
        获取优先级最高的N个反应
        
        Args:
            reactions: 反应列表
            top_n: 返回数量
            
        Returns:
            优先级最高的N个反应
        """
        prioritized = self.prioritize_reactions(reactions)
        return [reaction for reaction, _, _ in prioritized[:top_n]]
    
    def categorize_reactions(self, reactions) -> Dict[str, List]:
        """
        按类别分类反应
        
        Returns:
            {
                'high_priority': [...],   # 高优先级 (分数 > 0.7)
                'medium_priority': [...], # 中优先级 (0.4 < 分数 <= 0.7)
                'low_priority': [...],    # 低优先级 (分数 <= 0.4)
            }
        """
        prioritized = self.prioritize_reactions(reactions)
        
        categories = {
            'high_priority': [],
            'medium_priority': [],
            'low_priority': []
        }
        
        for reaction, score, _ in prioritized:
            if score > 0.7:
                categories['high_priority'].append(reaction)
            elif score > 0.4:
                categories['medium_priority'].append(reaction)
            else:
                categories['low_priority'].append(reaction)
        
        return categories


def create_battery_reaction_filter(electrode_type: str = 'anode',
                                   temperature: float = 300.0,
                                   voltage: float = 0.0):
    """
    创建电池反应筛选器 (便捷函数)
    
    Args:
        electrode_type: 电极类型
        temperature: 温度 (K)
        voltage: 电极电势 (V)
        
    Returns:
        配置好的ReactionPrioritizer实例
    """
    from .improved_activation_energy import EnergyFilterConfig
    
    # 创建能量筛选配置
    energy_config = EnergyFilterConfig(
        temperature=temperature,
        voltage=voltage,
        electrode_type=electrode_type
    )
    
    # 创建SEI评分器
    sei_scorer = SEIRelevanceScorer(electrode_type=electrode_type)
    
    # 创建优先级排序器
    prioritizer = ReactionPrioritizer(
        sei_scorer=sei_scorer,
        thermodynamic_cutoff=energy_config.get_thermodynamic_cutoff(),
        kinetic_cutoff=energy_config.get_kinetic_cutoff()
    )
    
    return prioritizer
