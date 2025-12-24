"""
分子信息服务
通过 SMILES 获取分子名称、CAS 号等信息
"""

import logging
import requests
from typing import Dict, Optional
from rdkit import Chem
from rdkit.Chem import Descriptors
import time

logger = logging.getLogger(__name__)

# 缓存分子信息，避免重复查询
_molecule_info_cache: Dict[str, Dict] = {}

def get_molecule_info(smiles: str) -> Dict:
    """
    通过 SMILES 获取分子信息
    
    Args:
        smiles: 分子的 SMILES 表示
        
    Returns:
        包含分子信息的字典
    """
    if smiles in _molecule_info_cache:
        return _molecule_info_cache[smiles]
    
    # 基本信息
    info = {
        'smiles': smiles,
        'name': None,
        'cas_number': None,
        'molecular_formula': None,
        'molecular_weight': None,
        'iupac_name': None
    }
    
    try:
        # 使用 RDKit 计算基本信息
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            info['molecular_formula'] = Chem.rdMolDescriptors.CalcMolFormula(mol)
            info['molecular_weight'] = round(Descriptors.MolWt(mol), 2)
            
        # 尝试通过 PubChem 获取更多信息
        pubchem_info = get_pubchem_info(smiles)
        if pubchem_info:
            info.update(pubchem_info)
            
    except Exception as e:
        logger.warning(f"Failed to get molecule info for {smiles}: {e}")
    
    # 缓存结果
    _molecule_info_cache[smiles] = info
    return info

def get_pubchem_info(smiles: str, timeout: int = 5) -> Optional[Dict]:
    """
    通过 PubChem API 获取分子信息
    
    Args:
        smiles: 分子的 SMILES 表示
        timeout: 请求超时时间（秒）
        
    Returns:
        分子信息字典或 None
    """
    try:
        # 首先通过 SMILES 搜索 CID
        search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/JSON"
        
        response = requests.get(search_url, timeout=timeout)
        if response.status_code != 200:
            return None
            
        search_data = response.json()
        if 'IdentifierList' not in search_data or 'CID' not in search_data['IdentifierList']:
            return None
            
        cid = search_data['IdentifierList']['CID'][0]
        
        # 获取详细信息
        detail_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,IUPACName/JSON"
        
        response = requests.get(detail_url, timeout=timeout)
        if response.status_code != 200:
            return None
            
        detail_data = response.json()
        if 'PropertyTable' not in detail_data or 'Properties' not in detail_data['PropertyTable']:
            return None
            
        properties = detail_data['PropertyTable']['Properties'][0]
        
        # 获取同义词（包括常用名称和 CAS 号）
        synonyms_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
        
        synonyms = []
        cas_number = None
        common_name = None
        
        try:
            response = requests.get(synonyms_url, timeout=timeout)
            if response.status_code == 200:
                synonyms_data = response.json()
                if 'InformationList' in synonyms_data and 'Information' in synonyms_data['InformationList']:
                    synonyms = synonyms_data['InformationList']['Information'][0].get('Synonym', [])
                    
                    # 查找 CAS 号（格式：数字-数字-数字）
                    import re
                    cas_pattern = re.compile(r'^\d+-\d+-\d+$')
                    for synonym in synonyms:
                        if cas_pattern.match(synonym):
                            cas_number = synonym
                            break
                    
                    # 选择最短的非数字名称作为常用名称
                    for synonym in synonyms:
                        if not synonym.isdigit() and not cas_pattern.match(synonym):
                            if not common_name or len(synonym) < len(common_name):
                                common_name = synonym
                                
        except Exception as e:
            logger.debug(f"Failed to get synonyms for CID {cid}: {e}")
        
        return {
            'name': common_name,
            'cas_number': cas_number,
            'iupac_name': properties.get('IUPACName'),
            'molecular_formula': properties.get('MolecularFormula'),
            'molecular_weight': properties.get('MolecularWeight')
        }
        
    except Exception as e:
        logger.debug(f"PubChem lookup failed for {smiles}: {e}")
        return None

def get_simple_molecule_name(smiles: str) -> str:
    """
    获取简单的分子名称，优先使用常用名称
    
    Args:
        smiles: 分子的 SMILES 表示
        
    Returns:
        分子名称字符串
    """
    info = get_molecule_info(smiles)
    
    # 优先级：常用名称 > 分子式 > SMILES
    if info.get('name'):
        return info['name']
    elif info.get('molecular_formula'):
        return info['molecular_formula']
    else:
        return smiles

# 常见分子的硬编码映射（用于快速查找）
COMMON_MOLECULES = {
    'CCO': {'name': 'Ethanol', 'cas_number': '64-17-5', 'molecular_formula': 'C2H6O'},
    'CC(=O)O': {'name': 'Acetic acid', 'cas_number': '64-19-7', 'molecular_formula': 'C2H4O2'},
    'C1COC(=O)O1': {'name': 'Ethylene carbonate', 'cas_number': '96-49-1', 'molecular_formula': 'C3H4O3'},
    'O': {'name': 'Water', 'cas_number': '7732-18-5', 'molecular_formula': 'H2O'},
    'CO': {'name': 'Methanol', 'cas_number': '67-56-1', 'molecular_formula': 'CH4O'},
    'CC': {'name': 'Ethane', 'cas_number': '74-84-0', 'molecular_formula': 'C2H6'},
    'C': {'name': 'Methane', 'cas_number': '74-82-8', 'molecular_formula': 'CH4'},
}

def get_molecule_info_fast(smiles: str) -> Dict:
    """
    快速获取分子信息，优先使用硬编码映射
    
    Args:
        smiles: 分子的 SMILES 表示
        
    Returns:
        包含分子信息的字典
    """
    # 首先检查硬编码映射
    if smiles in COMMON_MOLECULES:
        info = COMMON_MOLECULES[smiles].copy()
        info['smiles'] = smiles
        return info
    
    # 否则使用完整的查找逻辑
    return get_molecule_info(smiles)
