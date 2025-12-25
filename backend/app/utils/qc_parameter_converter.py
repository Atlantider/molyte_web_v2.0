"""
QC引擎参数转换器

提供Gaussian和PySCF之间的参数兼容性转换
"""
from typing import Dict, Optional
import logging

logger = logging.getLogger(__name__)


# Gaussian → PySCF 基组映射
GAUSSIAN_TO_PYSCF_BASIS = {
    # Pople基组
    '6-31G*': '6-31G(d)',
    '6-31G**': '6-31G(d,p)',
    '6-31+G*': '6-31+G(d)',
    '6-31+G**': '6-31+G(d,p)',
    '6-311G*': '6-311G(d)',
    '6-311G**': '6-311G(d,p)',
    '6-311+G*': '6-311+G(d)',
    '6-311+G**': '6-311+G(d,p)',
    
    # Correlation-consistent基组 (PySCF使用小写)
    'cc-pVDZ': 'cc-pvdz',
    'cc-pVTZ': 'cc-pvtz',
    'cc-pVQZ': 'cc-pvqz',
    'aug-cc-pVDZ': 'aug-cc-pvdz',
    'aug-cc-pVTZ': 'aug-cc-pvtz',
    'aug-cc-pVQZ': 'aug-cc-pvqz',
    
    # Def2基组
    'def2-SVP': 'def2-svp',
    'def2-SVPD': 'def2-svpd',
    'def2-TZVP': 'def2-tzvp',
    'def2-TZVPD': 'def2-tzvpd',
    'def2-QZVP': 'def2-qzvp',
}


# 泛函名称映射 (Gaussian → PySCF)
FUNCTIONAL_MAP = {
    # 杂化泛函
    'B3LYP': 'b3lyp',
    'B3LYP5': 'b3lyp5',
    'PBE0': 'pbe0',
    'PBE1PBE': 'pbe0',
    'M06': 'm06',
    'M06-2X': 'm062x',
    'M06-HF': 'm06hf',
    'M06-L': 'm06l',
    
    # 远程校正泛函
    'wB97': 'wb97',
    'wB97X': 'wb97x',
    'wB97X-D': 'wb97x-d',
    'wB97X-D3': 'wb97x-d3',
    'CAM-B3LYP': 'camb3lyp',
    
    # GGA泛函
    'PBE': 'pbe',
    'BLYP': 'blyp',
    'BP86': 'bp86',
    
    # Meta-GGA
    'TPSS': 'tpss',
    'M11-L': 'm11l',
}


# 溶剂介电常数映射
SOLVENT_DIELECTRIC = {
    'water': 78.3553,
    'h2o': 78.3553,
    'acetonitrile': 35.688,
    'ch3cn': 35.688,
    'dmso': 46.826,
    'dimethylsulfoxide': 46.826,
    'ethanol': 24.852,
    'etoh': 24.852,
    'methanol': 32.613,
    'meoh': 32.613,
    'dichloromethane': 8.93,
    'dcm': 8.93,
    'chloroform': 4.7113,
    'chcl3': 4.7113,
    'thf': 7.4257,
    'tetrahydrofuran': 7.4257,
    'toluene': 2.374,
    'benzene': 2.2706,
}


def convert_basis_set(gaussian_basis: str) -> str:
    """
    将Gaussian基组名称转换为PySCF格式
    
    Args:
        gaussian_basis: Gaussian格式的基组名称
        
    Returns:
        PySCF格式的基组名称
    """
    # 直接查找映射
    if gaussian_basis in GAUSSIAN_TO_PYSCF_BASIS:
        pyscf_basis = GAUSSIAN_TO_PYSCF_BASIS[gaussian_basis]
        logger.debug(f"Converted basis set: {gaussian_basis} → {pyscf_basis}")
        return pyscf_basis
    
    # 如果不在映射中，尝试转换为小写
    lower_basis = gaussian_basis.lower()
    if lower_basis != gaussian_basis:
        logger.info(f"Using lowercase basis set: {gaussian_basis} → {lower_basis}")
        return lower_basis
    
    # 原样返回
    logger.warning(f"No mapping found for basis set: {gaussian_basis}, using as-is")
    return gaussian_basis


def convert_functional(gaussian_functional: str) -> str:
    """
    将Gaussian泛函名称转换为PySCF格式
    
    Args:
        gaussian_functional: Gaussian格式的泛函名称
        
    Returns:
        PySCF格式的泛函名称
    """
    # 直接查找映射
    if gaussian_functional in FUNCTIONAL_MAP:
        pyscf_functional = FUNCTIONAL_MAP[gaussian_functional]
        logger.debug(f"Converted functional: {gaussian_functional} → {pyscf_functional}")
        return pyscf_functional
    
    # 尝试小写
    lower_func = gaussian_functional.lower()
    if lower_func != gaussian_functional:
        logger.info(f"Using lowercase functional: {gaussian_functional} → {lower_func}")
        return lower_func
    
    # 原样返回
    logger.warning(f"No mapping found for functional: {gaussian_functional}, using as-is")
    return gaussian_functional


def get_solvent_dielectric(solvent_name: str) -> float:
    """
    获取溶剂的介电常数
    
    Args:
        solvent_name: 溶剂名称
        
    Returns:
        介电常数，默认为水的介电常数
    """
    if not solvent_name:
        return 78.3553  # 默认水
    
    solvent_lower = solvent_name.lower().strip()
    
    if solvent_lower in SOLVENT_DIELECTRIC:
        eps = SOLVENT_DIELECTRIC[solvent_lower]
        logger.debug(f"Solvent {solvent_name} → ε={eps}")
        return eps
    
    # 默认水
    logger.warning(f"Unknown solvent: {solvent_name}, using water ε=78.3553")
    return 78.3553


def get_minimal_basis(basis_set: str) -> str:
    """
    获取基组的最小版本(用于渐进式计算)
    
    Args:
        basis_set: 原始基组
        
    Returns:
        最小基组
    """
    # 基组降级映射
    minimal_map = {
        '6-311G(d,p)': '6-31G',
        '6-311G(d)': '6-31G',
        '6-31G(d,p)': '6-31G',
        '6-31G(d)': '6-31G',
        'cc-pvtz': 'cc-pvdz',
        'cc-pvqz': 'cc-pvdz',
        'aug-cc-pvtz': 'cc-pvdz',
        'aug-cc-pvqz': 'cc-pvdz',
        'def2-tzvp': 'def2-svp',
        'def2-qzvp': 'def2-svp',
    }
    
    basis_lower = basis_set.lower()
    if basis_lower in minimal_map:
        minimal = minimal_map[basis_lower]
        logger.debug(f"Minimal basis for {basis_set}: {minimal}")
        return minimal
    
    # 默认STO-3G
    return 'sto-3g'


def should_use_xtb_preopt(molecule_info: Dict) -> bool:
    """
    判断是否应该使用XTB预优化
    
    Args:
        molecule_info: 分子信息字典
        
    Returns:
        是否使用XTB预优化
    """
    # 大分子建议使用
    num_atoms = molecule_info.get('num_atoms', 0)
    if num_atoms > 20:
        return True
    
    # 复杂结构(环、多重键)
    smiles = molecule_info.get('smiles', '')
    if '=' in smiles or '#' in smiles or any(c in smiles for c in '123456789'):
        return True
    
    # 带电体系
    if abs(molecule_info.get('charge', 0)) > 0:
        return True
    
    # 默认使用(除非是非常小的分子)
    return num_atoms > 5


def get_scf_convergence_params(molecule_info: Dict) -> Dict:
    """
    根据分子特征获取SCF收敛参数
    
    Args:
        molecule_info: 分子信息
        
    Returns:
        SCF参数字典
    """
    params = {
        'max_cycle': 150,
        'conv_tol': 1e-8,
    }
    
    num_atoms = molecule_info.get('num_atoms', 0)
    charge = molecule_info.get('charge', 0)
    multiplicity = molecule_info.get('multiplicity', 1)
    
    # 大分子: 增加循环数
    if num_atoms > 50:
        params['max_cycle'] = 300
        params['level_shift'] = 0.1
    elif num_atoms > 100:
        params['max_cycle'] = 500
        params['level_shift'] = 0.2
    
    # 带电体系: 使用level shift
    if abs(charge) > 0:
        params['level_shift'] = params.get('level_shift', 0) + 0.2
    
    # 开壳层: 调整参数
    if multiplicity > 1:
        params['max_cycle'] = max(params['max_cycle'], 200)
        # DIIS对开壳层更重要
        params['diis'] = True
    
    logger.debug(f"SCF params for molecule (N={num_atoms}, Q={charge}, M={multiplicity}): {params}")
    return params


class ParameterConverter:
    """统一的参数转换器"""
    
    @staticmethod
    def gaussian_to_pyscf_params(params: Dict) -> Dict:
        """
        Gaussian参数完整转换为PySCF参数
        
        Args:
            params: Gaussian风格的参数字典
            
        Returns:
            PySCF风格的参数字典
        """
        pyscf_params = params.copy()
        
        # 基组转换
        if 'basis_set' in pyscf_params:
            pyscf_params['basis_set'] = convert_basis_set(pyscf_params['basis_set'])
        
        # 泛函转换
        if 'functional' in pyscf_params:
            pyscf_params['functional'] = convert_functional(pyscf_params['functional'])
        
        # multiplicity → spin
        if 'multiplicity' in pyscf_params:
            pyscf_params['spin'] = pyscf_params['multiplicity'] - 1
        
        # 溶剂处理
        if pyscf_params.get('solvent_model') in ['pcm', 'smd']:
            solvent_name = pyscf_params.get('solvent_name', 'water')
            pyscf_params['solvent_dielectric'] = get_solvent_dielectric(solvent_name)
        
        return pyscf_params
    
    @staticmethod
    def prepare_for_engine(params: Dict, engine: str) -> Dict:
        """
        为指定引擎准备参数
        
        Args:
            params: 原始参数
            engine: 目标引擎 ('gaussian' or 'pyscf')
            
        Returns:
            转换后的参数
        """
        if engine == 'pyscf':
            return ParameterConverter.gaussian_to_pyscf_params(params)
        
        # Gaussian保持原样
        return params.copy()
