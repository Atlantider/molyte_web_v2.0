"""
化学工具函数模块

提供基于化学物理原理的智能计算工具，包括：
- 量子化学近似（Hückel理论）
- 键能和键解离能计算
- 环张力分析
- 自由基稳定性评估
- 电负性和极性分析

所有方法都基于真实的化学原理，避免硬编码。
"""

import numpy as np
from typing import List, Dict, Any, Optional, Tuple, Set
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import math


# ==================== 量子化学近似 ====================

def calculate_huckel_energy(mol: Chem.Mol, atom_indices: Optional[List[int]] = None) -> Dict[str, Any]:
    """
    使用Hückel理论计算π电子体系的能级
    
    原理：
    - 只考虑π电子（sp2杂化的碳原子）
    - 使用简化的哈密顿矩阵
    - α = 库仑积分（设为0作为参考）
    - β = 共振积分（相邻原子间）
    
    Args:
        mol: RDKit分子对象
        atom_indices: 要考虑的原子索引列表（如果为None，自动检测π体系）
    
    Returns:
        {
            'homo_energy': HOMO能级,
            'lumo_energy': LUMO能级,
            'homo_index': HOMO轨道索引,
            'lumo_index': LUMO轨道索引,
            'all_energies': 所有能级列表,
            'pi_atoms': π电子原子列表
        }
    """
    # 1. 找到π电子体系
    if atom_indices is None:
        pi_atoms = find_pi_system(mol)
    else:
        pi_atoms = atom_indices
    
    if len(pi_atoms) == 0:
        return {
            'homo_energy': None,
            'lumo_energy': None,
            'homo_index': None,
            'lumo_index': None,
            'all_energies': [],
            'pi_atoms': []
        }
    
    # 2. 构建Hückel矩阵（邻接矩阵）
    n = len(pi_atoms)
    H = np.zeros((n, n))
    
    # 对角元素 = α（设为0）
    # 非对角元素 = β（相邻原子）或 0（不相邻）
    
    atom_idx_map = {atom_idx: i for i, atom_idx in enumerate(pi_atoms)}
    
    for i, atom_idx in enumerate(pi_atoms):
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in atom_idx_map:
                j = atom_idx_map[neighbor_idx]
                # 检查是否有π键
                bond = mol.GetBondBetweenAtoms(atom_idx, neighbor_idx)
                if bond and (bond.GetBondType() == Chem.BondType.DOUBLE or 
                           bond.GetBondType() == Chem.BondType.AROMATIC):
                    H[i, j] = 1.0  # β的单位（归一化）
    
    # 3. 求解本征值（能级）
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    
    # 能级从低到高排序
    sorted_indices = np.argsort(eigenvalues)
    energies = eigenvalues[sorted_indices]
    
    # 4. 计算π电子数
    num_pi_electrons = sum(1 for atom_idx in pi_atoms 
                          if mol.GetAtomWithIdx(atom_idx).GetHybridization() == Chem.HybridizationType.SP2)
    
    # 每个轨道最多容纳2个电子
    num_occupied = (num_pi_electrons + 1) // 2
    
    # 5. 确定HOMO和LUMO
    if num_occupied > 0 and num_occupied <= len(energies):
        homo_index = num_occupied - 1
        homo_energy = energies[homo_index]
    else:
        homo_index = None
        homo_energy = None
    
    if num_occupied < len(energies):
        lumo_index = num_occupied
        lumo_energy = energies[lumo_index]
    else:
        lumo_index = None
        lumo_energy = None
    
    return {
        'homo_energy': float(homo_energy) if homo_energy is not None else None,
        'lumo_energy': float(lumo_energy) if lumo_energy is not None else None,
        'homo_index': homo_index,
        'lumo_index': lumo_index,
        'all_energies': energies.tolist(),
        'pi_atoms': pi_atoms,
        'num_pi_electrons': num_pi_electrons
    }


def find_pi_system(mol: Chem.Mol) -> List[int]:
    """
    找到分子中的π电子体系
    
    Returns:
        π电子原子的索引列表
    """
    pi_atoms = []
    
    for atom in mol.GetAtoms():
        # sp2杂化的碳、氮、氧
        if atom.GetHybridization() in [Chem.HybridizationType.SP2, Chem.HybridizationType.SP]:
            pi_atoms.append(atom.GetIdx())
        # 芳香原子
        elif atom.GetIsAromatic():
            pi_atoms.append(atom.GetIdx())
    
    return pi_atoms


def estimate_homo_lumo_gap(mol: Chem.Mol) -> float:
    """
    估算HOMO-LUMO能隙
    
    较小的能隙表示：
    - 更容易发生电子转移
    - 更高的反应活性
    
    Returns:
        能隙（β单位），如果无法计算则返回None
    """
    huckel_result = calculate_huckel_energy(mol)
    
    if huckel_result['homo_energy'] is not None and huckel_result['lumo_energy'] is not None:
        return huckel_result['lumo_energy'] - huckel_result['homo_energy']
    
    return None


# ==================== 电负性和极性 ====================

def calculate_pauling_electronegativity(atom: Chem.Atom) -> float:
    """
    获取Pauling电负性
    
    基于元素周期表的经验值
    """
    electronegativity_table = {
        'H': 2.20, 'C': 2.55, 'N': 3.04, 'O': 3.44, 'F': 3.98,
        'P': 2.19, 'S': 2.58, 'Cl': 3.16, 'Br': 2.96, 'I': 2.66,
        'Li': 0.98, 'Na': 0.93, 'K': 0.82, 'Mg': 1.31, 'Ca': 1.00,
        'Al': 1.61, 'Si': 1.90, 'B': 2.04
    }
    
    symbol = atom.GetSymbol()
    return electronegativity_table.get(symbol, 2.0)  # 默认值


def calculate_bond_polarity(bond: Chem.Bond, mol: Chem.Mol) -> float:
    """
    计算键的极性
    
    极性 = |电负性差|
    
    Returns:
        极性值（0-4范围）
    """
    atom1 = bond.GetBeginAtom()
    atom2 = bond.GetEndAtom()
    
    en1 = calculate_pauling_electronegativity(atom1)
    en2 = calculate_pauling_electronegativity(atom2)
    
    return abs(en1 - en2)


# ==================== 键能计算 ====================

def estimate_bond_dissociation_energy(bond: Chem.Bond, mol: Chem.Mol) -> float:
    """
    估算键解离能（BDE）
    
    BDE = 基础键能 - 自由基稳定化能
    
    原理：
    - 不同类型的键有不同的基础键能
    - 形成的自由基越稳定，BDE越低
    - 考虑共轭、超共轭、共振效应
    
    Args:
        bond: 要计算的键
        mol: 分子对象
    
    Returns:
        BDE（kcal/mol）
    """
    # 1. 基础键能
    base_bde = get_base_bond_energy(bond)
    
    # 2. 计算两端原子断裂后的自由基稳定化能
    atom1 = bond.GetBeginAtom()
    atom2 = bond.GetEndAtom()
    
    stabilization1 = calculate_radical_stabilization(atom1, mol, bond)
    stabilization2 = calculate_radical_stabilization(atom2, mol, bond)
    
    # 3. BDE = 基础能量 - 稳定化能
    bde = base_bde - (stabilization1 + stabilization2)
    
    return bde


def get_base_bond_energy(bond: Chem.Bond) -> float:
    """
    获取基础键能（kcal/mol）
    
    基于键类型和原子类型
    """
    bond_type = bond.GetBondType()
    atom1 = bond.GetBeginAtom().GetSymbol()
    atom2 = bond.GetEndAtom().GetSymbol()
    
    # 键能表（kcal/mol）
    bond_energies = {
        # C-X键
        ('C', 'C', Chem.BondType.SINGLE): 83,
        ('C', 'C', Chem.BondType.DOUBLE): 146,
        ('C', 'C', Chem.BondType.TRIPLE): 200,
        ('C', 'H', Chem.BondType.SINGLE): 99,
        ('C', 'O', Chem.BondType.SINGLE): 86,
        ('C', 'O', Chem.BondType.DOUBLE): 179,
        ('C', 'N', Chem.BondType.SINGLE): 73,
        ('C', 'N', Chem.BondType.DOUBLE): 147,
        ('C', 'N', Chem.BondType.TRIPLE): 213,
        ('C', 'F', Chem.BondType.SINGLE): 116,
        ('C', 'Cl', Chem.BondType.SINGLE): 81,
        ('C', 'Br', Chem.BondType.SINGLE): 68,
        ('C', 'I', Chem.BondType.SINGLE): 51,
        
        # O-X键
        ('O', 'O', Chem.BondType.SINGLE): 35,  # 过氧键（弱）
        ('O', 'H', Chem.BondType.SINGLE): 111,
        
        # N-X键
        ('N', 'N', Chem.BondType.SINGLE): 38,  # 偶氮键（弱）
        ('N', 'N', Chem.BondType.DOUBLE): 100,
        ('N', 'H', Chem.BondType.SINGLE): 93,
        ('N', 'O', Chem.BondType.SINGLE): 55,
        
        # S-X键
        ('S', 'S', Chem.BondType.SINGLE): 54,  # 二硫键（弱）
        ('S', 'H', Chem.BondType.SINGLE): 81,
    }
    
    # 标准化原子顺序
    key1 = tuple(sorted([atom1, atom2])) + (bond_type,)
    key2 = (atom1, atom2, bond_type)
    key3 = (atom2, atom1, bond_type)
    
    # 查找键能
    if key1 in bond_energies:
        return bond_energies[key1]
    elif key2 in bond_energies:
        return bond_energies[key2]
    elif key3 in bond_energies:
        return bond_energies[key3]
    else:
        # 默认值（基于键类型）
        if bond_type == Chem.BondType.SINGLE:
            return 80
        elif bond_type == Chem.BondType.DOUBLE:
            return 145
        elif bond_type == Chem.BondType.TRIPLE:
            return 200
        else:
            return 80


def calculate_radical_stabilization(atom: Chem.Atom, mol: Chem.Mol, 
                                   breaking_bond: Chem.Bond) -> float:
    """
    计算自由基稳定化能
    
    稳定化因素：
    1. 超共轭（与相邻C-H键）
    2. 共振（烯丙位、苄位）
    3. 碳的级数（叔 > 仲 > 伯）
    
    Args:
        atom: 将形成自由基的原子
        mol: 分子对象
        breaking_bond: 正在断裂的键
    
    Returns:
        稳定化能（kcal/mol，正值表示稳定）
    """
    stabilization = 0.0
    
    # 只考虑碳自由基
    if atom.GetSymbol() != 'C':
        return 0.0
    
    # 1. 检查是否是烯丙位或苄位（共振稳定）
    is_allylic = False
    is_benzylic = False
    
    for neighbor in atom.GetNeighbors():
        if neighbor.GetIdx() == breaking_bond.GetBeginAtomIdx() or \
           neighbor.GetIdx() == breaking_bond.GetEndAtomIdx():
            continue  # 跳过正在断裂的键
        
        # 检查是否连接到sp2碳（双键）
        if neighbor.GetHybridization() == Chem.HybridizationType.SP2:
            bond_to_neighbor = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
            if bond_to_neighbor and bond_to_neighbor.GetBondType() == Chem.BondType.SINGLE:
                # 检查neighbor是否在双键中
                for nn in neighbor.GetNeighbors():
                    if nn.GetIdx() != atom.GetIdx():
                        bond_nn = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nn.GetIdx())
                        if bond_nn and bond_nn.GetBondType() == Chem.BondType.DOUBLE:
                            is_allylic = True
                            break
        
        # 检查是否连接到芳香环
        if neighbor.GetIsAromatic():
            is_benzylic = True
    
    if is_allylic:
        stabilization += 15.0  # 烯丙位自由基稳定化
    if is_benzylic:
        stabilization += 13.0  # 苄位自由基稳定化
    
    # 2. 碳的级数（超共轭效应）
    # 计算连接的碳原子数（不包括正在断裂的键）
    num_carbon_neighbors = 0
    for neighbor in atom.GetNeighbors():
        if neighbor.GetIdx() == breaking_bond.GetBeginAtomIdx() or \
           neighbor.GetIdx() == breaking_bond.GetEndAtomIdx():
            continue
        if neighbor.GetSymbol() == 'C':
            num_carbon_neighbors += 1
    
    # 叔碳自由基 > 仲碳 > 伯碳
    if num_carbon_neighbors == 3:  # 叔碳
        stabilization += 7.0
    elif num_carbon_neighbors == 2:  # 仲碳
        stabilization += 4.0
    elif num_carbon_neighbors == 1:  # 伯碳
        stabilization += 0.0
    
    return stabilization


# ==================== 环张力计算 ====================

def calculate_ring_strain(ring_atoms: List[int], mol: Chem.Mol) -> float:
    """
    计算环张力能量
    
    改进：直接使用经验值，更准确
    
    Args:
        ring_atoms: 环中原子的索引列表
        mol: 分子对象
    
    Returns:
        环张力能量（kcal/mol）
    """
    ring_size = len(ring_atoms)
    
    # 直接使用经验值（更准确）
    return get_empirical_ring_strain(ring_size)


def get_empirical_ring_strain(ring_size: int) -> float:
    """
    获取经验环张力值
    
    基于实验测量的环张力能
    
    Returns:
        环张力（kcal/mol）
    """
    empirical_strains = {
        3: 27.5,  # 环丙烷
        4: 26.3,  # 环丁烷
        5: 6.2,   # 环戊烷
        6: 0.1,   # 环己烷（几乎无张力）
        7: 6.2,   # 环庚烷
        8: 9.7,   # 环辛烷
    }
    
    return empirical_strains.get(ring_size, 0.0)


# ==================== 自由基检测 ====================

def detect_radicals(mol: Chem.Mol) -> List[int]:
    """
    检测分子中的自由基位点
    
    Returns:
        含有未配对电子的原子索引列表
    """
    radical_atoms = []
    
    for atom in mol.GetAtoms():
        if atom.GetNumRadicalElectrons() > 0:
            radical_atoms.append(atom.GetIdx())
    
    return radical_atoms


def is_radical(mol: Chem.Mol) -> bool:
    """检查分子是否是自由基"""
    return len(detect_radicals(mol)) > 0


# ==================== 共轭体系分析 ====================

def find_conjugated_systems(mol: Chem.Mol) -> List[List[int]]:
    """
    找到分子中的共轭体系
    
    共轭体系：交替的单双键或芳香体系
    
    Returns:
        共轭体系列表，每个体系是原子索引的列表
    """
    conjugated_systems = []
    visited = set()
    
    for atom in mol.GetAtoms():
        if atom.GetIdx() in visited:
            continue
        
        # 检查是否是sp2或芳香原子
        if atom.GetHybridization() == Chem.HybridizationType.SP2 or atom.GetIsAromatic():
            # 开始BFS找到整个共轭体系
            system = []
            queue = [atom.GetIdx()]
            
            while queue:
                current_idx = queue.pop(0)
                if current_idx in visited:
                    continue
                
                visited.add(current_idx)
                system.append(current_idx)
                
                current_atom = mol.GetAtomWithIdx(current_idx)
                for neighbor in current_atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited:
                        # 检查是否也是sp2或芳香
                        if neighbor.GetHybridization() == Chem.HybridizationType.SP2 or neighbor.GetIsAromatic():
                            queue.append(neighbor_idx)
            
            if len(system) > 1:
                conjugated_systems.append(system)
    
    return conjugated_systems


# ==================== 导出函数 ====================

__all__ = [
    # 量子化学
    'calculate_huckel_energy',
    'find_pi_system',
    'estimate_homo_lumo_gap',
    
    # 电负性和极性
    'calculate_pauling_electronegativity',
    'calculate_bond_polarity',
    
    # 键能
    'estimate_bond_dissociation_energy',
    'get_base_bond_energy',
    'calculate_radical_stabilization',
    
    # 环张力
    'calculate_ring_strain',
    'get_empirical_ring_strain',
    
    # 自由基
    'detect_radicals',
    'is_radical',
    
    # 共轭体系
    'find_conjugated_systems',
]
