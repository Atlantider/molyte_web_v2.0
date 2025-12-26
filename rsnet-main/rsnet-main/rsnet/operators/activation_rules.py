"""
算符激活规则配置

基于化学和电化学原理的算符激活规则，包括：
- 必需驱动力
- 增强驱动力
- 分子特征要求
- 激活权重
"""

from typing import Dict, List, Any


# 算符激活规则配置
OPERATOR_ACTIVATION_RULES = {
    'electron_transfer': {
        'required_drives': ['electrochemical'],
        'enhancing_drives': ['oxidation', 'reduction'],
        'required_features': [],  # π体系在运行时检测
        'molecular_checks': {
            'has_pi_systems': 0.5,  # 有π体系增加0.5分
            'has_heteroatoms': 0.3,  # 有杂原子增加0.3分
            'has_conjugated_systems': 0.4,  # 有共轭体系增加0.4分
        },
        'weight': 1.0,
        'description': '电子转移反应（氧化还原）',
        'chemical_principle': 'HOMO/LUMO能级匹配，电化学电位驱动'
    },
    
    'radical_reaction': {
        'required_drives': [],  # 可以在多种条件下发生
        'enhancing_drives': ['radical_environment', 'high_temperature', 'surface_reaction'],
        'required_features': [],
        'molecular_checks': {
            'has_ch_bonds': 0.4,  # 有C-H键用于抽取
            'is_radical': 0.8,  # 已经是自由基
            'has_double_bonds': 0.3,  # 有双键用于加成
        },
        'weight': 0.9,
        'description': '自由基反应（H抽取、偶联）',
        'chemical_principle': 'BDE驱动，自由基稳定性'
    },
    
    'ring_opening': {
        'required_drives': [],
        'enhancing_drives': ['thermal', 'high_temperature', 'ring_strain'],
        'required_features': ['small_rings'],
        'molecular_checks': {
            'has_small_rings': 1.0,  # 必须有小环
            'ring_strain_high': 0.5,  # 高张力环
        },
        'weight': 0.8,
        'description': '开环反应（环张力释放）',
        'chemical_principle': 'Baeyer张力理论，环张力能量释放'
    },
    
    'elimination': {
        'required_drives': ['thermal', 'high_temperature'],
        'enhancing_drives': [],
        'required_features': [],
        'molecular_checks': {
            'has_leaving_group': 0.8,  # 有离去基团
            'has_beta_h': 0.6,  # 有β-H
            'has_polar_bonds': 0.3,  # 有极性键
        },
        'weight': 0.7,
        'description': '消除反应（E1/E2机理）',
        'chemical_principle': 'Zaitsev规则，离去基团能力'
    },
    
    'coordination': {
        'required_drives': [],
        'enhancing_drives': ['li_coordination', 'solution_phase', 'surface_reaction'],
        'required_features': [],
        'molecular_checks': {
            'has_metal_ion': 0.9,  # 有金属离子
            'has_ligand_sites': 0.7,  # 有配体位点
        },
        'weight': 0.75,
        'description': '配位反应（金属-配体）',
        'chemical_principle': 'HSAB理论，配位场理论'
    },
    
    'polymerization': {
        'required_drives': [],
        'enhancing_drives': ['radical_environment', 'sei_formation', 'surface_reaction'],
        'required_features': [],
        'molecular_checks': {
            'is_polymerizable': 0.8,  # 可聚合（有C=C）
            'is_initiator': 0.9,  # 是引发剂
        },
        'weight': 0.65,
        'description': '聚合反应（自由基聚合）',
        'chemical_principle': '自由基链式反应，SEI膜形成'
    },
    
    'decomposition': {
        'required_drives': ['high_temperature'],
        'enhancing_drives': [],
        'required_features': [],
        'molecular_checks': {
            'has_weak_bonds': 0.8,  # 有弱键
            'temperature_sufficient': 0.5,  # 温度足够
        },
        'weight': 0.7,
        'description': '热分解反应',
        'chemical_principle': 'BDE，温度依赖的键断裂'
    },
    
    # 原有算符
    'hydrogen_transfer': {
        'required_drives': [],
        'enhancing_drives': ['thermal', 'acidic_environment', 'basic_environment'],
        'required_features': [],
        'molecular_checks': {
            'has_acidic_h': 0.6,
            'has_h_acceptor': 0.6,
        },
        'weight': 0.6,
        'description': '氢转移反应',
        'chemical_principle': '酸碱质子转移'
    },
    
    'bond_breaking': {
        'required_drives': [],
        'enhancing_drives': ['thermal', 'high_temperature', 'weak_bond_cleavage'],
        'required_features': [],
        'molecular_checks': {
            'has_weak_bonds': 0.7,
            'has_small_rings': 0.4,
        },
        'weight': 0.65,
        'description': '键断裂反应',
        'chemical_principle': '键能，热力学驱动'
    },
    
    'cyclization': {
        'required_drives': [],
        'enhancing_drives': ['thermal', 'ring_closure_favorable'],
        'required_features': [],
        'molecular_checks': {
            'has_reactive_ends': 0.7,
            'chain_length_suitable': 0.5,
        },
        'weight': 0.6,
        'description': '环化反应',
        'chemical_principle': '熵驱动，环张力平衡'
    },
    
    'addition': {
        'required_drives': [],
        'enhancing_drives': ['thermal', 'pi_system_reaction'],
        'required_features': [],
        'molecular_checks': {
            'has_double_bonds': 0.7,
            'has_nucleophile': 0.5,
            'has_electrophile': 0.5,
        },
        'weight': 0.7,
        'description': '加成反应',
        'chemical_principle': 'π键反应性，亲核/亲电加成'
    },
    
    'rearrangement': {
        'required_drives': [],
        'enhancing_drives': ['thermal', 'high_temperature'],
        'required_features': [],
        'molecular_checks': {
            'has_rearrangeable_structure': 0.6,
            'thermodynamically_favorable': 0.4,
        },
        'weight': 0.5,
        'description': '重排反应',
        'chemical_principle': '热力学稳定性，分子内重组'
    },
    
    'redox': {
        'required_drives': [],
        'enhancing_drives': ['electrochemical', 'oxidation', 'reduction'],
        'required_features': [],
        'molecular_checks': {
            'has_redox_center': 0.8,
            'voltage_appropriate': 0.6,
        },
        'weight': 0.85,
        'description': '氧化还原反应',
        'chemical_principle': '电化学电位，电子转移'
    },
}


# 驱动力权重配置（用于综合评分）
DRIVE_WEIGHTS = {
    # 电化学驱动力（高权重）
    'electrochemical': 1.0,
    'oxidation': 0.9,
    'reduction': 0.9,
    
    # 热力学驱动力
    'thermal': 0.7,
    'high_temperature': 0.8,
    'cryogenic': 0.6,
    
    # 分子特征驱动力
    'ring_strain': 0.85,
    'weak_bond_cleavage': 0.75,
    'polar_bond_cleavage': 0.7,
    'pi_system_reaction': 0.7,
    
    # 界面和环境驱动力
    'surface_reaction': 0.8,
    'sei_formation': 0.85,
    'cei_formation': 0.8,
    'radical_environment': 0.75,
    
    # 配位驱动力
    'li_coordination': 0.8,
    
    # 其他
    'high_pressure': 0.5,
    'gas_phase': 0.4,
    'solution_phase': 0.5,
}


# 分子特征检测函数映射
FEATURE_DETECTORS = {
    'has_pi_systems': lambda tags: tags.get('has_pi_systems', False),
    'has_heteroatoms': lambda tags: tags.get('has_heteroatoms', False),
    'has_conjugated_systems': lambda tags: len(tags.get('pi_systems', [])) > 1,
    'has_ch_bonds': lambda tags: True,  # 大多数有机分子都有C-H键
    'is_radical': lambda tags: False,  # 需要特殊检测
    'has_double_bonds': lambda tags: tags.get('has_pi_systems', False),
    'has_small_rings': lambda tags: tags.get('has_small_rings', False),
    'ring_strain_high': lambda tags: any(
        len(ring) <= 4 for ring in tags.get('small_rings', [])
    ),
    'has_leaving_group': lambda tags: tags.get('has_polar_bonds', False),
    'has_beta_h': lambda tags: True,  # 简化检测
    'has_polar_bonds': lambda tags: tags.get('has_polar_bonds', False),
    'has_metal_ion': lambda tags: False,  # 需要特殊检测
    'has_ligand_sites': lambda tags: tags.get('has_heteroatoms', False),
    'is_polymerizable': lambda tags: tags.get('has_pi_systems', False),
    'is_initiator': lambda tags: False,  # 需要特殊检测
    'has_weak_bonds': lambda tags: tags.get('has_weak_bonds', False),
    'temperature_sufficient': lambda tags: True,  # 在环境检查中处理
    'has_acidic_h': lambda tags: tags.get('has_acidic_hydrogens', False),
    'has_h_acceptor': lambda tags: len(tags.get('h_bond_acceptors', [])) > 0,
    'has_reactive_ends': lambda tags: tags.get('has_pi_systems', False) or tags.get('has_polar_bonds', False),
    'chain_length_suitable': lambda tags: True,  # 简化
    'has_nucleophile': lambda tags: tags.get('has_heteroatoms', False),
    'has_electrophile': lambda tags: tags.get('has_pi_systems', False),
    'has_rearrangeable_structure': lambda tags: tags.get('has_small_rings', False),
    'thermodynamically_favorable': lambda tags: True,  # 需要能量计算
    'has_redox_center': lambda tags: tags.get('has_heteroatoms', False) or tags.get('has_pi_systems', False),
    'voltage_appropriate': lambda tags: True,  # 在环境检查中处理
}


def get_activation_rule(operator_name: str) -> Dict[str, Any]:
    """获取算符的激活规则"""
    return OPERATOR_ACTIVATION_RULES.get(operator_name, {})


def get_drive_weight(drive_name: str) -> float:
    """获取驱动力的权重"""
    return DRIVE_WEIGHTS.get(drive_name, 0.5)
