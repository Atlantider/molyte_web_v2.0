"""
Operator registry for RSNet.

This module provides a centralized registry for managing reaction operators
and their activation based on driving forces and molecular features.
"""

from typing import List, Dict, Any, Optional, Set
from ..core.molecule import Molecule
from ..core.environment import Environment
from .base import BaseOperator
from .hydrogen_transfer import HydrogenTransferOperator
from .bond_breaking import BondBreakingOperator
from .cyclization import CyclizationOperator
from .addition import AdditionOperator
from .rearrangement import RearrangementOperator
from .redox import RedoxOperator
# 扩展算符 - 每个算符都在独立的文件中
from .electron_transfer import ElectronTransferOperator
from .ring_opening import RingOpeningOperator
from .elimination import EliminationOperator
from .radical_reaction import RadicalReactionOperator
from .coordination import CoordinationOperator
from .polymerization import PolymerizationOperator
from .decomposition import DecompositionOperator
from .clustering import ClusteringOperator
from .substitution import SubstitutionOperator


class OperatorRegistry:
    """
    Registry for managing reaction operators.
    
    This class provides centralized management of all available operators
    and methods to select appropriate operators based on conditions.
    """
    
    def __init__(self):
        """Initialize the operator registry."""
        self._operators = {}
        self._operator_instances = {}
        self._setup_default_operators()
    
    def _setup_default_operators(self):
        """Set up default operators with their configurations."""
        
        # 1. Hydrogen transfer operator
        self.register_operator(
            "hydrogen_transfer",
            HydrogenTransferOperator,
            {
                'max_distance': 4.0,
                'min_distance': 1.5,
                'require_3d': False
            }
        )
        
        # 2. Bond breaking operator
        self.register_operator(
            "bond_breaking",
            BondBreakingOperator,
            {
                'min_bond_strength': 50.0,  # kcal/mol
                'allow_ring_breaking': True,
                'prefer_weak_bonds': True
            }
        )
        
        # 3. Cyclization operator
        self.register_operator(
            "cyclization",
            CyclizationOperator,
            {
                'min_ring_size': 3,
                'max_ring_size': 8,
                'allow_spiro': True,
                'allow_fused': True
            }
        )
        
        # 4. Addition operator
        self.register_operator(
            "addition",
            AdditionOperator,
            {
                'allow_michael_addition': True,
                'allow_cycloaddition': True,
                'require_activation': False
            }
        )
        
        # 5. Rearrangement operator
        self.register_operator(
            "rearrangement",
            RearrangementOperator,
            {
                'min_ring_size': 3,
                'max_ring_size': 8,
                'allow_skeletal': True
            }
        )
        
        # 6. Redox operator
        self.register_operator(
            "redox",
            RedoxOperator,
            {
                'max_electron_transfer': 2,
                'allow_radical_formation': True
            }
        )
        
        # 7. Electron transfer operator (扩展)
        self.register_operator(
            "electron_transfer",
            ElectronTransferOperator,
            {
                'oxidation_threshold': 3.0,  # V
                'reduction_threshold': 1.0   # V
            }
        )
        
        # 8. Ring opening operator (扩展)
        self.register_operator(
            "ring_opening",
            RingOpeningOperator,
            {
                'min_ring_size': 3,
                'max_ring_size': 5,
                'prefer_strained': True
            }
        )
        
        # 9. Elimination operator (扩展)
        self.register_operator(
            "elimination",
            EliminationOperator,
            {
                'min_temperature': 350.0,
                'allow_beta_elimination': True,
                'allow_dehydration': True
            }
        )
        
        # 10. Radical reaction operator (扩展)
        self.register_operator(
            "radical_reaction",
            RadicalReactionOperator,
            {
                'allow_h_abstraction': True,
                'allow_coupling': True,
                'allow_chain_growth': True
            }
        )
        
        # 11. Coordination operator (扩展)
        self.register_operator(
            "coordination",
            CoordinationOperator,
            {
                'allow_exchange': True,
                'allow_bridging': True,
                'max_coordination_number': 6
            }
        )
        
        # 12. Polymerization operator (扩展)
        self.register_operator(
            "polymerization",
            PolymerizationOperator,
            {
                'max_chain_length': 10,
                'allow_branching': True,
                'allow_crosslinking': False
            }
        )
        
        # 13. Decomposition operator (扩展)
        self.register_operator(
            "decomposition",
            DecompositionOperator,
            {
                'min_temperature': 500.0,
                'allow_homolytic': True,
                'allow_heterolytic': True
            }
        )

        # 14. Clustering operator (扩展 - 多分子)
        self.register_operator(
            "clustering",
            ClusteringOperator,
            {
                'max_coordination': 4
            }
        )
        
        # 15. Substitution operator (Sn2/Transesterification)
        self.register_operator(
            "substitution",
            SubstitutionOperator,
            {
                'allow_sn2': True,
                'allow_transesterification': True
            }
        )
    
        # 16. Salt Dissociation operator
        from .salt_dissociation import SaltDissociationOperator
        self.register_operator(
            "salt_dissociation",
            SaltDissociationOperator,
            {}
        )
    
    def register_operator(
        self, 
        name: str, 
        operator_class: type, 
        config: Optional[Dict[str, Any]] = None
    ):
        """
        Register a new operator.
        
        Args:
            name: Operator name
            operator_class: Operator class
            config: Configuration parameters
        """
        self._operators[name] = {
            'class': operator_class,
            'config': config or {}
        }
    
    def get_operator(self, name: str) -> Optional[BaseOperator]:
        """
        Get an operator instance by name.
        
        Args:
            name: Operator name
            
        Returns:
            Operator instance or None if not found
        """
        if name not in self._operators:
            return None
        
        # Create instance if not already created
        if name not in self._operator_instances:
            operator_info = self._operators[name]
            self._operator_instances[name] = operator_info['class'](
                config=operator_info['config']
            )
        
        return self._operator_instances[name]
    
    def get_all_operators(self) -> List[BaseOperator]:
        """
        Get all registered operators.
        
        Returns:
            List of all operator instances
        """
        operators = []
        for name in self._operators.keys():
            operator = self.get_operator(name)
            if operator:
                operators.append(operator)
        return operators
    
    def get_active_operators(
        self, 
        molecules: List[Molecule], 
        environment: Environment
    ) -> List[BaseOperator]:
        """
        Get operators that can be applied to the given molecules and environment.
        
        Args:
            molecules: List of molecules
            environment: Reaction environment
            
        Returns:
            List of applicable operators
        """
        active_operators = []
        
        for name in self._operators.keys():
            operator = self.get_operator(name)
            if operator and operator.can_apply(molecules, environment):
                active_operators.append(operator)
        
        return active_operators
    
    def get_operators_by_driving_force(
        self, 
        driving_forces: Dict[str, float],
        threshold: float = 0.1
    ) -> List[BaseOperator]:
        """
        Get operators based on active driving forces.
        
        Args:
            driving_forces: Dictionary of driving force strengths
            threshold: Minimum strength threshold
            
        Returns:
            List of relevant operators
        """
        relevant_operators = []
        
        # Define operator-driving force relationships
        operator_drives = {
            'hydrogen_transfer': ['thermal', 'radical_environment'],
            'bond_breaking': ['thermal', 'high_temperature', 'radical_environment', 'ring_strain'],
            'cyclization': ['thermal', 'ring_strain', 'pi_system_reaction'],
            'addition': ['electrochemical', 'solution_phase', 'functional_group_reaction'],
            'rearrangement': ['thermal', 'high_temperature', 'radical_environment'],
            'redox': ['oxidation', 'reduction', 'electrochemical', 'radical_environment'],
            'substitution': ['thermal', 'solution_phase']
        }
        
        for operator_name, required_drives in operator_drives.items():
            # Check if any required driving force is above threshold
            if any(driving_forces.get(drive, 0.0) >= threshold for drive in required_drives):
                operator = self.get_operator(operator_name)
                if operator:
                    relevant_operators.append(operator)
        
        return relevant_operators
    
    def get_operator_priority(
        self, 
        operator: BaseOperator, 
        molecules: List[Molecule], 
        environment: Environment
    ) -> float:
        """
        Calculate priority score for an operator.
        
        Args:
            operator: Operator to evaluate
            molecules: List of molecules
            environment: Reaction environment
            
        Returns:
            Priority score (0.0 to 1.0)
        """
        if not operator.can_apply(molecules, environment):
            return 0.0
        
        # Get driving forces
        drives = environment.get_active_drives()
        
        # Base priority based on operator type
        base_priorities = {
            'BondBreaking': 0.8,  # High priority for decomposition
            'HydrogenTransfer': 0.7,  # Common reaction
            'Redox': 0.6,  # Important in electrochemical systems
            'Addition': 0.5,  # Moderate priority
            'Cyclization': 0.4,  # Lower priority
            'Rearrangement': 0.3,  # Lowest priority (usually slower)
            'Substitution': 0.6
        }
        
        base_priority = base_priorities.get(operator.name, 0.5)
        
        # Adjust based on driving forces
        drive_bonus = 0.0
        if operator.name == 'Redox' and (drives.get('oxidation') or drives.get('reduction')):
            drive_bonus += 0.2
        elif operator.name == 'BondBreaking' and drives.get('high_temperature'):
            drive_bonus += 0.2
        elif operator.name == 'HydrogenTransfer' and drives.get('radical_environment'):
            drive_bonus += 0.1
        
        return min(1.0, base_priority + drive_bonus)
    
    def get_recommended_operators(
        self, 
        molecules: List[Molecule], 
        environment: Environment,
        max_operators: int = 5
    ) -> List[BaseOperator]:
        """
        Get recommended operators sorted by priority.
        
        Args:
            molecules: List of molecules
            environment: Reaction environment
            max_operators: Maximum number of operators to return
            
        Returns:
            List of recommended operators sorted by priority
        """
        active_operators = self.get_active_operators(molecules, environment)
        
        # Calculate priorities and sort
        operator_priorities = []
        for operator in active_operators:
            priority = self.get_operator_priority(operator, molecules, environment)
            operator_priorities.append((operator, priority))
        
        # Sort by priority (descending)
        operator_priorities.sort(key=lambda x: x[1], reverse=True)
        
        # Return top operators
        return [op for op, _ in operator_priorities[:max_operators]]
    
    def list_operators(self) -> Dict[str, Dict[str, Any]]:
        """
        List all registered operators with their information.
        
        Returns:
            Dictionary of operator information
        """
        operator_info = {}
        
        for name, info in self._operators.items():
            operator = self.get_operator(name)
            operator_info[name] = {
                'class': info['class'].__name__,
                'config': info['config'],
                'description': operator.description if operator else "No description"
            }
        
        return operator_info
    
    def get_active_operators_smart(
        self,
        molecules: List[Molecule],
        environment: Environment,
        min_score: float = 0.1
    ) -> List[tuple]:
        """
        智能获取激活的算符及其优先级分数
        
        基于化学和电化学原理的智能激活系统：
        1. 环境驱动力匹配
        2. 分子特征匹配
        3. 化学可行性验证
        
        Args:
            molecules: 分子列表
            environment: 反应环境
            min_score: 最小激活分数阈值
            
        Returns:
            List[(operator, score)] 按分数降序排列
        """
        from .activation_rules import get_activation_rule, FEATURE_DETECTORS
        from ..features.driving_forces import get_driving_forces
        from ..utils.structure_analysis import get_structure_tags
        from ..utils.chemistry_tools import detect_radicals, is_radical
        
        # 1. 获取环境驱动力
        env_drives = environment.get_active_drives()
        drive_strengths = get_driving_forces(molecules, environment)
        
        # 2. 获取分子特征
        mol_features = []
        for mol in molecules:
            tags = get_structure_tags(mol)
            # 添加额外的自由基检测
            tags['is_radical'] = is_radical(mol.rdkit_mol)
            mol_features.append(tags)
        
        # 3. 评估每个算符
        operator_scores = []
        for op_name in self._operators.keys():
            operator = self.get_operator(op_name)
            if not operator:
                continue
            
            # 获取激活规则
            rule = get_activation_rule(op_name)
            if not rule:
                continue
            
            # 计算激活分数
            score = self._calculate_activation_score(
                operator, op_name, rule, molecules, environment,
                env_drives, drive_strengths, mol_features
            )
            
            if score >= min_score:
                operator_scores.append((operator, score))
        
        # 4. 按优先级排序
        operator_scores.sort(key=lambda x: x[1], reverse=True)
        
        return operator_scores
    
    def _calculate_activation_score(
        self,
        operator: BaseOperator,
        operator_name: str,
        rule: Dict[str, Any],
        molecules: List[Molecule],
        environment: Environment,
        env_drives: Dict[str, bool],
        drive_strengths: Dict[str, float],
        mol_features: List[Dict]
    ) -> float:
        """
        计算算符激活分数
        
        分数 = 环境匹配度 × 分子匹配度 × 化学可行性 × 权重
        
        基于化学原理：
        - 电化学反应需要电化学驱动力
        - 热反应需要足够的温度
        - 分子结构必须满足反应要求
        """
        # 1. 环境匹配度
        env_score = self._evaluate_environment_match(
            rule, env_drives, drive_strengths, environment
        )
        
        if env_score == 0.0:
            return 0.0  # 环境不匹配，直接返回
        
        # 2. 分子匹配度
        mol_score = self._evaluate_molecule_match(
            rule, molecules, mol_features
        )
        
        # 3. 化学可行性（调用operator.can_apply）
        feasibility = 1.0 if operator.can_apply(molecules, environment) else 0.0
        
        # 4. 算符权重
        weight = rule.get('weight', 0.5)
        
        # 综合分数
        total_score = env_score * 0.4 + mol_score * 0.3 + feasibility * 0.2 + weight * 0.1
        
        return total_score
    
    def _evaluate_environment_match(
        self,
        rule: Dict[str, Any],
        env_drives: Dict[str, bool],
        drive_strengths: Dict[str, float],
        environment: Environment
    ) -> float:
        """
        评估环境匹配度
        
        化学原理：
        - 必需驱动力必须存在
        - 增强驱动力提高匹配度
        - 驱动力强度影响分数
        """
        from .activation_rules import get_drive_weight
        
        required_drives = rule.get('required_drives', [])
        enhancing_drives = rule.get('enhancing_drives', [])
        
        # 检查必需驱动力
        if required_drives:
            # 至少一个必需驱动力必须激活
            has_required = any(env_drives.get(drive, False) for drive in required_drives)
            if not has_required:
                return 0.0
        
        # 计算驱动力分数
        total_strength = 0.0
        total_weight = 0.0
        
        # 必需驱动力（高权重）
        for drive in required_drives:
            if env_drives.get(drive, False):
                strength = drive_strengths.get(drive, 0.5)
                weight = get_drive_weight(drive)
                total_strength += strength * weight * 1.5  # 必需驱动力加权
                total_weight += weight * 1.5
        
        # 增强驱动力（正常权重）
        for drive in enhancing_drives:
            if env_drives.get(drive, False):
                strength = drive_strengths.get(drive, 0.5)
                weight = get_drive_weight(drive)
                total_strength += strength * weight
                total_weight += weight
        
        if total_weight == 0.0:
            return 0.3  # 没有特定驱动力，给予基础分数
        
        return min(1.0, total_strength / total_weight)
    
    def _evaluate_molecule_match(
        self,
        rule: Dict[str, Any],
        molecules: List[Molecule],
        mol_features: List[Dict]
    ) -> float:
        """
        评估分子匹配度
        
        化学原理：
        - 分子必须具有反应所需的结构特征
        - 多个分子取最高匹配度
        """
        from .activation_rules import FEATURE_DETECTORS
        
        molecular_checks = rule.get('molecular_checks', {})
        
        if not molecular_checks:
            return 0.7  # 没有特定要求，给予中等分数
        
        # 对每个分子评估
        max_score = 0.0
        
        for tags in mol_features:
            mol_score = 0.0
            total_weight = 0.0
            
            for check_name, check_weight in molecular_checks.items():
                # 使用特征检测器
                detector = FEATURE_DETECTORS.get(check_name)
                if detector and detector(tags):
                    mol_score += check_weight
                total_weight += check_weight
            
            if total_weight > 0.0:
                normalized_score = mol_score / total_weight
                max_score = max(max_score, normalized_score)
        
        return max_score



# Global operator registry instance
OPERATOR_REGISTRY = OperatorRegistry()


def get_active_operators(
    molecules: List[Molecule], 
    environment: Environment
) -> List[BaseOperator]:
    """
    Convenience function to get active operators.
    
    Args:
        molecules: List of molecules
        environment: Reaction environment
        
    Returns:
        List of applicable operators
    """
    return OPERATOR_REGISTRY.get_active_operators(molecules, environment)


def get_recommended_operators(
    molecules: List[Molecule], 
    environment: Environment,
    max_operators: int = 5
) -> List[BaseOperator]:
    """
    Convenience function to get recommended operators.
    
    Args:
        molecules: List of molecules
        environment: Reaction environment
        max_operators: Maximum number of operators to return
        
    Returns:
        List of recommended operators
    """
    return OPERATOR_REGISTRY.get_recommended_operators(
        molecules, environment, max_operators
    )


def get_active_operators_smart(
    molecules: List[Molecule],
    environment: Environment,
    min_score: float = 0.1
) -> List[tuple]:
    """
    智能获取激活的算符及其优先级分数（便捷函数）
    
    Args:
        molecules: 分子列表
        environment: 反应环境
        min_score: 最小激活分数阈值
        
    Returns:
        List[(operator, score)] 按分数降序排列
    """
    return OPERATOR_REGISTRY.get_active_operators_smart(
        molecules, environment, min_score
    )

