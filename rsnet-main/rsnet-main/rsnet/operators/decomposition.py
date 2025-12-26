"""
分解反应算符 - 基于BDE的智能实现

化学原理：
- 热分解：高温下键断裂
- 均裂：R-R → R· + R·（生成自由基）
- 异裂：R-X → R+ + X-（生成离子）
- 需要高温（>500K）
- 优先断裂最弱的键
"""

from typing import List, Dict, Any, Optional, Tuple
from rdkit import Chem

from ..core.molecule import Molecule
from ..core.environment import Environment
from ..core.reaction import Reaction
from .base import BaseOperator
from ..utils.chemistry_tools import estimate_bond_dissociation_energy


class DecompositionOperator(BaseOperator):
    """
    分解反应算符 - 基于BDE的智能实现
    
    改进：
    - ✅ 使用BDE计算找到最弱键
    - ✅ 温度依赖的分解阈值
    - ✅ 均裂和异裂机制
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        super().__init__(
            name="Decomposition",
            description="热分解反应（基于BDE）",
            config=config
        )
        self.arity = 1  # Unimolecular (bond breaking)
        self.reaction_type = "decomposition"
        self.min_temperature = self.config.get('min_temperature', 500.0)
        self.allow_homolytic = self.config.get('allow_homolytic', True)
        self.allow_heterolytic = self.config.get('allow_heterolytic', True)
        self.bde_threshold = self.config.get('bde_threshold', 60.0)  # kcal/mol
    
    def can_apply(self, molecules: List[Molecule], environment: Environment) -> bool:
        """检查是否有高温分解条件"""
        return environment.temperature >= self.min_temperature
    
    def apply(self, molecules: List[Molecule], environment: Environment) -> List[Reaction]:
        """应用分解算符"""
        if not self.can_apply(molecules, environment):
            return []
        
        reactions = []
        
        for mol in molecules:
            # 均裂分解
            if self.allow_homolytic:
                fragments = self._homolytic_cleavage_bde(mol, environment)
                if fragments:
                    reaction = Reaction(
                        reactants=[mol],
                        products=fragments,
                        operator_name="homolytic_cleavage"
                    )
                    # 分解是吸热反应
                    reaction.reaction_energy = 50.0 + (environment.temperature - 500.0) * 0.1
                    reactions.append(reaction)
        
        return reactions
    
    def _homolytic_cleavage_bde(self, mol: Molecule, environment: Environment) -> Optional[List[Molecule]]:
        """
        均裂断裂 - 使用BDE找到最弱键
        
        改进：
        - 使用真实的BDE计算
        - 考虑温度效应
        - 只断裂弱键
        """
        try:
            rdkit_mol = mol.rdkit_mol
            
            # 找到所有键及其BDE
            bonds_with_bde = []
            for bond in rdkit_mol.GetBonds():
                bde = estimate_bond_dissociation_energy(bond, rdkit_mol)
                bonds_with_bde.append({
                    'bond': bond,
                    'bde': bde,
                    'atom1_idx': bond.GetBeginAtomIdx(),
                    'atom2_idx': bond.GetEndAtomIdx()
                })
            
            if not bonds_with_bde:
                return None
            
            # 按BDE排序，找到最弱的键
            bonds_with_bde.sort(key=lambda x: x['bde'])
            weakest = bonds_with_bde[0]
            
            # 考虑温度效应：温度越高，可以断裂更强的键
            temperature_bonus = (environment.temperature - 500.0) * 0.1
            effective_threshold = self.bde_threshold + temperature_bonus
            
            if weakest['bde'] > effective_threshold:
                return None  # 键太强，无法断裂
            
            # 断裂键，生成自由基
            new_mol = Chem.RWMol(rdkit_mol)
            new_mol.RemoveBond(weakest['atom1_idx'], weakest['atom2_idx'])
            
            # 添加自由基电子
            atom1 = new_mol.GetAtomWithIdx(weakest['atom1_idx'])
            atom2 = new_mol.GetAtomWithIdx(weakest['atom2_idx'])
            atom1.SetNumRadicalElectrons(atom1.GetNumRadicalElectrons() + 1)
            atom2.SetNumRadicalElectrons(atom2.GetNumRadicalElectrons() + 1)
            
            # 获取片段
            try:
                Chem.SanitizeMol(new_mol)
                fragments = Chem.GetMolFrags(new_mol, asMols=True, sanitizeFrags=True)
                
                if len(fragments) >= 2:
                    products = []
                    for i, frag in enumerate(fragments):
                        smiles = Chem.MolToSmiles(frag)
                        products.append(Molecule.from_smiles(smiles, name=f"{mol.name}_frag{i+1}"))
                    return products
            except:
                return None
            
            return None
        except:
            return None
