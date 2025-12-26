"""
RSNet Compute Module Improvements
改进模块初始化 - 导出所有改进功能
"""

from .improved_activation_energy import (
    ImprovedActivationEnergyEstimator,
    EnergyFilterConfig,
    ReactionType,
    BEPParameters
)

from .sei_relevance_scorer import (
    SEIRelevanceScorer,
    ReactionPrioritizer,
    SEIComponentType,
    SEIComponent,
    create_battery_reaction_filter
)

# 方便导入
__all__ = [
    # 活化能估算
    'ImprovedActivationEnergyEstimator',
    'EnergyFilterConfig',
    'ReactionType',
    'BEPParameters',
    
    # SEI相关性评分
    'SEIRelevanceScorer',
    'ReactionPrioritizer',
    'SEIComponentType',
    'SEIComponent',
    'create_battery_reaction_filter',
]


def get_improved_reaction_screener(electrode_type: str = 'anode',
                                   temperature: float = 300.0,
                                   voltage: float = 0.0):
    """
    获取改进的反应筛选器
    
    便捷函数，返回完整配置的筛选器组合
    
    Args:
        electrode_type: 电极类型 ('anode' 或 'cathode')
        temperature: 温度 (K)
        voltage: 电极电势 (V)
        
    Returns:
        dict: 包含所有改进组件的字典
    """
    # 活化能估算器
    ea_estimator = ImprovedActivationEnergyEstimator(temperature=temperature)
    
    # 能量筛选配置
    energy_config = EnergyFilterConfig(
        temperature=temperature,
        voltage=voltage,
        electrode_type=electrode_type
    )
    
    # SEI相关性评分器
    sei_scorer = SEIRelevanceScorer(electrode_type=electrode_type)
    
    # 反应优先级排序器
    prioritizer = ReactionPrioritizer(
        sei_scorer=sei_scorer,
        thermodynamic_cutoff=energy_config.get_thermodynamic_cutoff(),
        kinetic_cutoff=energy_config.get_kinetic_cutoff()
    )
    
    return {
        'ea_estimator': ea_estimator,
        'energy_config': energy_config,
        'sei_scorer': sei_scorer,
        'prioritizer': prioritizer,
        
        # 便捷方法
        'estimate_activation_energy': ea_estimator.estimate,
        'check_feasibility': energy_config.is_feasible,
        'score_relevance': sei_scorer.score_reaction,
        'prioritize_reactions': prioritizer.prioritize_reactions,
    }
