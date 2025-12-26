"""
自由基反应算符 - 基于BDE的智能实现

化学原理：
- 使用键解离能（BDE）预测H抽取位点
- 自由基稳定性：烯丙位 > 苄位 > 叔碳 > 仲碳 > 伯碳
- 自由基偶联：检测自由基位点并连接
- 考虑位阻和电子效应
"""

from typing import List, Dict, Any, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem

from ..core.molecule import Molecule
from ..core.environment import Environment
from ..core.reaction import Reaction
from .base import BaseOperator
from ..utils.chemistry_tools import (
    estimate_bond_dissociation_energy,
    detect_radicals,
    is_radical,
    calculate_radical_stabilization,
    find_conjugated_systems
)


class RadicalReactionOperator(BaseOperator):
    """
    自由基反应算符 - 基于BDE的智能实现
    
    改进：
    - ✅ 使用BDE计算预测最易断裂的C-H键
    - ✅ 实现自由基偶联机制
    - ✅ 考虑自由基稳定性
    - ✅ 智能位点选择
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        super().__init__(
            name="RadicalReaction",
            description="自由基反应（基于BDE计算）",
            config=config
        )
        self.arity = 2  # Mixed: Unimolecular (abstraction) + Bimolecular (coupling)
        self.reaction_type = "radical"
        self.allow_h_abstraction = self.config.get('allow_h_abstraction', True)
        self.allow_coupling = self.config.get('allow_coupling', True)
        self.allow_chain_growth = self.config.get('allow_chain_growth', True)
        self.bde_threshold = self.config.get('bde_threshold', 95.0)  # kcal/mol
    
    def can_apply(self, molecules: List[Molecule], environment: Environment) -> bool:
        """检查是否有自由基环境"""
        # 需要高温或自由基环境
        if environment.temperature < 400.0:
            drives = environment.get_active_drives()
            if not drives.get('radical_environment', False):
                return False
        
        return True
    
    def apply(self, molecules: List[Molecule], environment: Environment) -> List[Reaction]:
        """应用自由基算符"""
        if not self.can_apply(molecules, environment):
            return []
        
        reactions = []
        
        # 1. 氢抽取反应（单分子）
        if self.allow_h_abstraction and len(molecules) >= 1:
            for mol in molecules:
                abstracted = self._hydrogen_abstraction_bde(mol, environment)
                if abstracted:
                    reaction = Reaction(
                        reactants=[mol],
                        products=[abstracted],
                        operator_name="h_abstraction"
                    )
                    # 能量 = BDE（吸热）
                    reaction.reaction_energy = 10.0  # 简化
                    reactions.append(reaction)
        
        # 2. 自由基偶联（双分子）
        if self.allow_coupling and len(molecules) >= 2:
            # Optimization: Identify all radicals first to avoid checking every pair
            radicals = [m for m in molecules if is_radical(m.rdkit_mol)]
            
            # Limit number of radicals to consider for coupling to avoid N^2 explosion
            # If we have too many radicals, prioritize smaller/lighter ones (likely more mobile)
            MAX_RADICALS_TO_COUPLE = 10
            if len(radicals) > MAX_RADICALS_TO_COUPLE:
                radicals.sort(key=lambda m: m.num_heavy_atoms)
                radicals = radicals[:MAX_RADICALS_TO_COUPLE]
            
            # Check pairwise coupling for selected radicals
            for i, mol1 in enumerate(radicals):
                for mol2 in radicals[i:]:
                    # Check for "irrelevant" radicals?
                    # e.g. very large polymers coupling might be sterically hindered
                    if mol1.num_heavy_atoms + mol2.num_heavy_atoms > 30:
                        continue 

                    coupled = self._radical_coupling(mol1, mol2)
                    if coupled:
                        reactants = [mol1, mol2] if mol1 != mol2 else [mol1]
                        reaction = Reaction(
                            reactants=reactants,
                            products=[coupled],
                            operator_name="radical_coupling"
                        )
                        reaction.reaction_energy = -80.0  # 偶联放热
                        reactions.append(reaction)
        
        return reactions
    
    def _hydrogen_abstraction_bde(self, mol: Molecule, environment: Environment) -> Optional[Molecule]:
        """
        基于BDE的氢抽取
        
        改进：
        - 计算所有C-H键的BDE
        - 选择BDE最低的键
        - 考虑温度效应
        """
        try:
            rdkit_mol = mol.rdkit_mol
            
            # 找到所有C-H键及其BDE
            ch_bonds = []
            for bond in rdkit_mol.GetBonds():
                atom1 = bond.GetBeginAtom()
                atom2 = bond.GetEndAtom()
                
                # 检查是否是C-H键
                if (atom1.GetSymbol() == 'C' and atom2.GetSymbol() == 'H') or \
                   (atom1.GetSymbol() == 'H' and atom2.GetSymbol() == 'C'):
                    
                    # 计算BDE
                    bde = estimate_bond_dissociation_energy(bond, rdkit_mol)
                    
                    # 找到碳原子
                    c_atom = atom1 if atom1.GetSymbol() == 'C' else atom2
                    h_atom = atom2 if atom1.GetSymbol() == 'C' else atom1
                    
                    ch_bonds.append({
                        'bond': bond,
                        'c_idx': c_atom.GetIdx(),
                        'h_idx': h_atom.GetIdx(),
                        'bde': bde
                    })
            
            if not ch_bonds:
                return None
            
            # 选择BDE最低的键（最易断裂）
            ch_bonds.sort(key=lambda x: x['bde'])
            weakest = ch_bonds[0]
            
            # 检查是否在阈值内（温度越高，可以断裂更强的键）
            temperature_bonus = (environment.temperature - 400.0) * 0.05
            effective_threshold = self.bde_threshold + temperature_bonus
            
            if weakest['bde'] > effective_threshold:
                return None
            
            # 创建自由基
            new_mol = Chem.RWMol(rdkit_mol)
            c_atom = new_mol.GetAtomWithIdx(weakest['c_idx'])
            
            # 移除氢原子
            new_mol.RemoveAtom(weakest['h_idx'])
            
            # 添加自由基电子到碳
            # 注意：移除原子后索引会改变，需要重新获取
            # 简化处理：直接设置自由基
            try:
                Chem.SanitizeMol(new_mol)
                # 找到原来的碳原子（现在索引可能变了）
                for atom in new_mol.GetAtoms():
                    if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() < rdkit_mol.GetAtomWithIdx(weakest['c_idx']).GetTotalNumHs():
                        atom.SetNumRadicalElectrons(1)
                        break
                
                smiles = Chem.MolToSmiles(new_mol)
                return Molecule.from_smiles(smiles, name=f"{mol.name}_radical")
            except:
                return None
                
        except:
            return None
    
    def _radical_coupling(self, mol1: Molecule, mol2: Molecule) -> Optional[Molecule]:
        """
        自由基偶联 - 智能实现
        
        改进：
        - 检测自由基位点
        - 连接两个自由基
        - 考虑位阻
        """
        try:
            rdkit_mol1 = mol1.rdkit_mol
            rdkit_mol2 = mol2.rdkit_mol
            
            # 检测自由基位点
            radicals1 = detect_radicals(rdkit_mol1)
            radicals2 = detect_radicals(rdkit_mol2)
            
            if not radicals1 and not radicals2:
                return None  # 没有自由基
            
            # 情况1：两个都是自由基
            if radicals1 and radicals2:
                return self._couple_two_radicals(mol1, mol2, radicals1[0], radicals2[0])
            
            # 情况2：只有一个是自由基（自由基加成）
            elif radicals1:
                return self._radical_addition(mol1, mol2, radicals1[0])
            else:
                return self._radical_addition(mol2, mol1, radicals2[0])
                
        except:
            return None
    
    def _couple_two_radicals(self, mol1: Molecule, mol2: Molecule, 
                            rad1_idx: int, rad2_idx: int) -> Optional[Molecule]:
        """
        两个自由基偶联形成σ键
        
        R1· + R2· → R1-R2
        """
        try:
            # 使用RDKit的组合分子功能
            rdkit_mol1 = mol1.rdkit_mol
            rdkit_mol2 = mol2.rdkit_mol
            
            # 创建组合分子
            combo = Chem.CombineMols(rdkit_mol1, rdkit_mol2)
            editable = Chem.EditableMol(combo)
            
            # mol2的原子索引需要偏移
            offset = rdkit_mol1.GetNumAtoms()
            rad2_idx_in_combo = rad2_idx + offset
            
            # 添加新键连接两个自由基位点
            editable.AddBond(rad1_idx, rad2_idx_in_combo, Chem.BondType.SINGLE)
            
            # 获取新分子
            new_mol = editable.GetMol()
            
            # 移除自由基电子
            new_mol_rw = Chem.RWMol(new_mol)
            new_mol_rw.GetAtomWithIdx(rad1_idx).SetNumRadicalElectrons(0)
            new_mol_rw.GetAtomWithIdx(rad2_idx_in_combo).SetNumRadicalElectrons(0)
            
            try:
                Chem.SanitizeMol(new_mol_rw)
                smiles = Chem.MolToSmiles(new_mol_rw)
                return Molecule.from_smiles(smiles, name=f"{mol1.name}_{mol2.name}_coupled")
            except:
                return None
                
        except:
            return None
    
    def _radical_addition(self, radical_mol: Molecule, target_mol: Molecule, 
                         rad_idx: int) -> Optional[Molecule]:
        """
        自由基加成到双键
        
        R· + CH2=CH2 → R-CH2-CH2·
        """
        try:
            rdkit_radical = radical_mol.rdkit_mol
            rdkit_target = target_mol.rdkit_mol
            
            # 找到目标分子中的双键
            double_bonds = []
            for bond in rdkit_target.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    atom1 = bond.GetBeginAtom()
                    atom2 = bond.GetEndAtom()
                    # 只考虑C=C双键
                    if atom1.GetSymbol() == 'C' and atom2.GetSymbol() == 'C':
                        double_bonds.append({
                            'bond': bond,
                            'atom1_idx': atom1.GetIdx(),
                            'atom2_idx': atom2.GetIdx()
                        })
            
            if not double_bonds:
                return None  # 没有双键可以加成
            
            # 选择第一个双键（可以改进为选择最活泼的）
            db = double_bonds[0]
            
            # 组合分子
            combo = Chem.CombineMols(rdkit_radical, rdkit_target)
            editable = Chem.EditableMol(combo)
            
            # 计算偏移
            offset = rdkit_radical.GetNumAtoms()
            
            # 自由基加成到双键的一端
            # R· + C=C → R-C-C·
            editable.AddBond(rad_idx, db['atom1_idx'] + offset, Chem.BondType.SINGLE)
            
            # 将双键变为单键
            editable.RemoveBond(db['atom1_idx'] + offset, db['atom2_idx'] + offset)
            editable.AddBond(db['atom1_idx'] + offset, db['atom2_idx'] + offset, Chem.BondType.SINGLE)
            
            new_mol = editable.GetMol()
            new_mol_rw = Chem.RWMol(new_mol)
            
            # 移除原自由基电子，添加新自由基电子到另一端
            new_mol_rw.GetAtomWithIdx(rad_idx).SetNumRadicalElectrons(0)
            new_mol_rw.GetAtomWithIdx(db['atom2_idx'] + offset).SetNumRadicalElectrons(1)
            
            try:
                Chem.SanitizeMol(new_mol_rw)
                smiles = Chem.MolToSmiles(new_mol_rw)
                return Molecule.from_smiles(smiles, name=f"{radical_mol.name}_add_{target_mol.name}")
            except:
                return None
                
        except:
            return None
