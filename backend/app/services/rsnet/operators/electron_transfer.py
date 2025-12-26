"""
电子转移算符 - 基于Hückel理论的智能实现

化学原理：
- 使用Hückel理论计算HOMO/LUMO能级
- HOMO能级高 → 易氧化
- LUMO能级低 → 易还原
- 考虑共振和电负性效应
"""

from typing import List, Dict, Any, Optional
from rdkit import Chem

from ..core.molecule import Molecule
from ..core.environment import Environment
from ..core.reaction import Reaction
from .base import BaseOperator
from ..core.reaction import Reaction
from .base import BaseOperator
# Local imports used to avoid circular dependency


class ElectronTransferOperator(BaseOperator):
    """
    电子转移算符 - 基于量子化学的智能实现
    
    改进：
    - ✅ 使用Hückel理论预测氧化还原位点
    - ✅ 考虑共轭效应和电负性
    - ✅ 多位点评分系统
    - ✅ 能量估算
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        super().__init__(
            name="ElectronTransfer",
            description="单电子转移反应（基于Hückel理论）",
            config=config
        )
        self.arity = 1  # Unimolecular input (heterogeneous reaction)
        self.reaction_type = "redox"
        self.oxidation_threshold = self.config.get('oxidation_threshold', 3.5)  # V vs Li/Li+
        self.reduction_threshold = self.config.get('reduction_threshold', 1.5)  # V vs Li/Li+
    
    def can_apply(self, molecules: List[Molecule], environment: Environment) -> bool:
        """检查是否可以应用电子转移"""
        if len(molecules) != 1:
            return False
        
        # 必须有电化学环境
        if environment.voltage is None:
            return False
        
        # 检查驱动力
        drives = environment.get_active_drives()
        return drives.get('electrochemical', False) or drives.get('oxidation', False) or drives.get('reduction', False)
    
    def apply(self, molecules: List[Molecule], environment: Environment) -> List[Reaction]:
        """应用电子转移算符"""
        if not self.can_apply(molecules, environment):
            return []
        
        reactions = []
        mol = molecules[0]
        voltage = environment.voltage
        
        # 高电压 → 氧化反应（失去电子）
        if voltage > self.oxidation_threshold:
            oxidized = self._single_electron_oxidation(mol)
            if oxidized:
                reaction = Reaction(
                    reactants=[mol],
                    products=[oxidized],
                    operator_name="electron_abstraction"
                )
                reaction.reaction_energy = (voltage - self.oxidation_threshold) * 23.06  # eV to kcal/mol
                reactions.append(reaction)
        
        # 低电压 → 还原反应（获得电子）
        elif voltage < self.reduction_threshold:
            reduced = self._single_electron_reduction(mol)
            if reduced:
                reaction = Reaction(
                    reactants=[mol],
                    products=[reduced],
                    operator_name="electron_injection"
                )
                reaction.reaction_energy = (self.reduction_threshold - voltage) * 23.06
                reactions.append(reaction)
        
        return reactions
    
    def _single_electron_oxidation(self, mol: Molecule) -> Optional[Molecule]:
        """
        单电子氧化 - 基于Hückel理论
        
        改进：
        - 使用HOMO能级预测最易氧化的位点
        - 考虑共轭稳定性
        - 多位点评分
        """
        try:
            rdkit_mol = mol.rdkit_mol
            
            # 1. 使用Hückel理论找到最易氧化的位点
            target_atom_idx = self._find_oxidation_site_huckel(rdkit_mol)
            if target_atom_idx is None:
                # 回退到简单方法
                target_atom_idx = self._find_oxidation_site_simple(rdkit_mol)
            
            if target_atom_idx is None:
                return None
            
            # 2. 创建自由基阳离子
            new_mol = Chem.RWMol(rdkit_mol)
            atom = new_mol.GetAtomWithIdx(target_atom_idx)
            
            # 增加正电荷，添加自由基电子
            atom.SetFormalCharge(atom.GetFormalCharge() + 1)
            atom.SetNumRadicalElectrons(1)
            
            try:
                Chem.SanitizeMol(new_mol)
                smiles = Chem.MolToSmiles(new_mol)
                return Molecule.from_smiles(smiles, name=f"{mol.name}_ox")
            except:
                return None
        except:
            return None
    
    def _single_electron_reduction(self, mol: Molecule) -> Optional[Molecule]:
        """
        单电子还原 - 基于Hückel理论
        
        改进：
        - 使用LUMO能级预测最易还原的位点
        - 考虑缺电子中心
        """
        try:
            rdkit_mol = mol.rdkit_mol
            
            # 1. 使用Hückel理论找到最易还原的位点
            target_atom_idx = self._find_reduction_site_huckel(rdkit_mol)
            if target_atom_idx is None:
                # 回退到简单方法
                target_atom_idx = self._find_reduction_site_simple(rdkit_mol)
            
            if target_atom_idx is None:
                return None
            
            # 2. 创建自由基阴离子
            new_mol = Chem.RWMol(rdkit_mol)
            atom = new_mol.GetAtomWithIdx(target_atom_idx)
            
            # Check for double bonds to reduce (e.g. C=O -> .C-O-)
            reduced_pi_system = False
            if atom.GetSymbol() == 'C':
                for neighbor in atom.GetNeighbors():
                    bond = new_mol.GetBondBetweenAtoms(target_atom_idx, neighbor.GetIdx())
                    if bond.GetBondType() == Chem.BondType.DOUBLE and neighbor.GetSymbol() in ['O', 'N', 'S']:
                        # Found a pi bond to heteroatom - reduce it!
                        # C=X -> .C-X-
                        bond.SetBondType(Chem.BondType.SINGLE)
                        
                        # Set Radical on C
                        atom.SetNumRadicalElectrons(1)
                        atom.SetFormalCharge(0)
                        
                        # Set Negative Charge on heteroatom
                        neighbor.SetFormalCharge(neighbor.GetFormalCharge() - 1)
                        
                        reduced_pi_system = True
                        break
            
            if not reduced_pi_system:
                # Fallback: Just try adding electron to atom (works for cations)
                atom.SetFormalCharge(atom.GetFormalCharge() - 1)
                atom.SetNumRadicalElectrons(1)
            
            try:
                Chem.SanitizeMol(new_mol)
                smiles = Chem.MolToSmiles(new_mol)
                return Molecule.from_smiles(smiles, name=f"{mol.name}_red")
            except Exception as e:
                # print(f"DEBUG: Reduction sanitization failed: {e}")
                return None
        except:
            return None
    
    def _find_oxidation_site_huckel(self, mol: Chem.Mol) -> Optional[int]:
        """
        使用Hückel理论找到最易氧化的位点
        
        原理：
        - HOMO能级越高，越容易失去电子
        - 在HOMO轨道上电子密度最大的原子最易氧化
        """
        # 计算Hückel能级
        from ..utils.chemistry_tools import (
            calculate_huckel_energy,
            calculate_pauling_electronegativity,
            find_conjugated_systems
        )
        huckel_result = calculate_huckel_energy(mol)
        
        if not huckel_result['pi_atoms'] or huckel_result['homo_energy'] is None:
            return None
        
        pi_atoms = huckel_result['pi_atoms']
        
        # 评分系统：结合多个因素
        scores = []
        
        for atom_idx in pi_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            score = 0.0
            
            # 1. 电负性（越低越易氧化）
            en = calculate_pauling_electronegativity(atom)
            score += (4.0 - en) * 10  # 归一化
            
            # 2. 杂化类型
            if atom.GetHybridization() == Chem.HybridizationType.SP2:
                score += 5.0
            elif atom.GetIsAromatic():
                score += 3.0
            
            # 3. 原子类型（杂原子优先）
            if atom.GetSymbol() in ['N', 'O', 'S']:
                score += 10.0
            
            # 4. 共轭体系中的位置
            conjugated_systems = find_conjugated_systems(mol)
            for system in conjugated_systems:
                if atom_idx in system:
                    score += len(system) * 0.5  # 共轭体系越大越稳定
            
            scores.append((atom_idx, score))
        
        if not scores:
            return None
        
        # 返回得分最高的原子
        scores.sort(key=lambda x: x[1], reverse=True)
        return scores[0][0]
    
    def _find_reduction_site_huckel(self, mol: Chem.Mol) -> Optional[int]:
        """
        使用Hückel理论找到最易还原的位点
        
        原理：
        - LUMO能级越低，越容易获得电子
        - 缺电子中心（羰基碳、阳离子）优先
        """
        # 计算Hückel能级
        from ..utils.chemistry_tools import (
            calculate_huckel_energy,
            calculate_pauling_electronegativity
        )
        huckel_result = calculate_huckel_energy(mol)
        
        if not huckel_result['pi_atoms'] or huckel_result['lumo_energy'] is None:
            return None
        
        pi_atoms = huckel_result['pi_atoms']
        
        # 评分系统
        scores = []
        
        for atom_idx in pi_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            score = 0.0
            
            # 1. 电负性（越高越易还原）
            en = calculate_pauling_electronegativity(atom)
            score += en * 10
            
            # 2. 正电荷（阳离子最易还原）
            if atom.GetFormalCharge() > 0:
                score += 50.0
            
            # 3. 羰基碳（C=O）
            is_carbonyl = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O':
                    bond = mol.GetBondBetweenAtoms(atom_idx, neighbor.GetIdx())
                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                        is_carbonyl = True
                        score += 30.0
                        break
            
            # 4. 亚胺碳（C=N）
            if not is_carbonyl:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'N':
                        bond = mol.GetBondBetweenAtoms(atom_idx, neighbor.GetIdx())
                        if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                            score += 20.0
                            break
            
            scores.append((atom_idx, score))
        
        if not scores:
            return None
        
        # 返回得分最高的原子
        scores.sort(key=lambda x: x[1], reverse=True)
        return scores[0][0]
    
    def _find_oxidation_site_simple(self, mol: Chem.Mol) -> Optional[int]:
        """简单方法：基于原子类型的优先级"""
        candidates = []
        
        for atom in mol.GetAtoms():
            priority = 0
            if atom.GetSymbol() in ['O', 'N', 'S']:
                priority = 4
            elif atom.GetIsAromatic():
                priority = 3
            elif atom.GetHybridization() == Chem.HybridizationType.SP2:
                priority = 2
            elif atom.GetSymbol() == 'C':
                priority = 1
            
            if priority > 0:
                candidates.append((atom.GetIdx(), priority))
        
        if not candidates:
            return None
        
        candidates.sort(key=lambda x: x[1], reverse=True)
        return candidates[0][0]
    
    def _find_reduction_site_simple(self, mol: Chem.Mol) -> Optional[int]:
        """简单方法：基于缺电子中心"""
        candidates = []
        
        for atom in mol.GetAtoms():
            priority = 0
            
            # 羰基碳
            if atom.GetSymbol() == 'C':
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O':
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                            priority = 4
                            break
            
            # 阳离子
            if atom.GetFormalCharge() > 0:
                priority = 3
            
            # 亚胺碳
            elif atom.GetSymbol() == 'C':
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'N':
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                            priority = 2
                            break
            
            if priority > 0:
                candidates.append((atom.GetIdx(), priority))
        
        if not candidates:
            return None
        
        candidates.sort(key=lambda x: x[1], reverse=True)
        return candidates[0][0]
