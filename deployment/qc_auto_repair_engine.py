"""
QC 自动修复规则引擎

这是一个完整的规则表驱动的自动修复系统，用于处理 Gaussian 计算中的常见错误。
规则按优先级从高到低排列，每个错误类型最多重试 2 次，超过总重试次数则标记为人工检查。

规则表结构：
- error_patterns: 错误特征关键词（支持多条 AND 组合）
- error_category: 错误分类
- priority: 优先级（1-10，数字越小优先级越高）
- auto_fix_strategies: 按顺序尝试的修复策略列表
- max_retries: 该错误类型最多重试次数
- requires_manual_check: 是否需要人工检查
"""

import re
import logging
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field
from enum import Enum

logger = logging.getLogger(__name__)


class ErrorCategory(Enum):
    """错误分类"""
    SCF_CONVERGENCE = "scf_convergence"
    GEOMETRY_OPTIMIZATION = "geometry_optimization"
    FREQUENCY = "frequency"
    INTEGRAL_PRECISION = "integral_precision"
    MEMORY_DISK = "memory_disk"
    SOLVATION = "solvation"
    INTERNAL_CRASH = "internal_crash"
    UNKNOWN = "unknown"


@dataclass
class RepairStrategy:
    """修复策略"""
    name: str  # 策略名称
    keywords_to_add: str  # 要添加的关键词
    keywords_to_remove: Optional[List[str]] = None  # 要移除的关键词
    restart_from_chk: bool = False  # 是否从 chk 文件重启
    use_cartesian: bool = False  # 是否强制使用笛卡尔坐标
    description: str = ""  # 策略描述
    # Slurm 参数修改（可选）
    slurm_memory_gb: Optional[int] = None  # 修改内存大小（GB），None 表示不修改
    slurm_time_seconds: Optional[int] = None  # 修改时间限制（秒），None 表示不修改
    slurm_cpus: Optional[int] = None  # 修改 CPU 数量，None 表示不修改


@dataclass
class ErrorRule:
    """错误规则"""
    error_patterns: List[str]  # 错误特征关键词列表（AND 组合）
    error_category: ErrorCategory
    priority: int  # 1-10，数字越小优先级越高
    description: str
    auto_fix_strategies: List[RepairStrategy]  # 按顺序尝试的修复策略
    max_retries: int = 2  # 该错误类型最多重试次数
    requires_manual_check: bool = False  # 是否需要人工检查


# ==================== 规则表定义 ====================

ERROR_RULES = [
    # ========== 优先级 1: 内部坐标系统崩溃（最严重） ==========
    ErrorRule(
        error_patterns=[
            r'Tors failed|Bend failed',
            r'FormBX had a problem',
            r'segmentation violation|segmentation fault'
        ],
        error_category=ErrorCategory.GEOMETRY_OPTIMIZATION,
        priority=1,
        description="内部坐标系统崩溃（Tors failed + FormBX + 段错误）",
        auto_fix_strategies=[
            RepairStrategy(
                name="笛卡尔坐标重启（第1级）",
                keywords_to_add="Opt=(Cartesian,MaxStep=5) SCF=XQC Integral=UltraFine NoSymm Geom=AllCheck Guess=Read",
                restart_from_chk=True,
                use_cartesian=True,
                description="从 chk 文件重启，强制使用笛卡尔坐标，减小步长"
            ),
            RepairStrategy(
                name="笛卡尔坐标重启（第2级）",
                keywords_to_add="Opt=(Cartesian,CalcFC,MaxStep=5) SCF=XQC Integral=UltraFine NoSymm Geom=AllCheck Guess=Read",
                restart_from_chk=True,
                use_cartesian=True,
                description="计算初始 Hessian，进一步稳定优化"
            )
        ],
        max_retries=2,
        requires_manual_check=True  # 失败后需要人工检查几何
    ),

    # ========== 优先级 2: 内部坐标在坐标转换/Hessian更新时崩溃 ==========
    ErrorRule(
        error_patterns=[
            r'Cartesian Forces.*Max.*RMS',
            r'Berny optimization',
            r'FormGI|FormBX|segmentation violation|segmentation fault'
        ],
        error_category=ErrorCategory.GEOMETRY_OPTIMIZATION,
        priority=2,
        description="内部坐标在坐标转换/Hessian更新时崩溃",
        auto_fix_strategies=[
            RepairStrategy(
                name="笛卡尔坐标重启（第1级）",
                keywords_to_add="Opt=(Cartesian,MaxStep=5) SCF=XQC Integral=UltraFine NoSymm Geom=AllCheck Guess=Read",
                restart_from_chk=True,
                use_cartesian=True,
                description="从 chk 文件重启，强制使用笛卡尔坐标"
            ),
            RepairStrategy(
                name="笛卡尔坐标重启（第2级）",
                keywords_to_add="Opt=(Cartesian,CalcFC,MaxStep=5) SCF=XQC Integral=UltraFine NoSymm Geom=AllCheck Guess=Read",
                restart_from_chk=True,
                use_cartesian=True,
                description="计算初始 Hessian"
            )
        ],
        max_retries=2,
        requires_manual_check=True
    ),

    # ========== 优先级 3: 一般的 Tors failed / FormBX 错误 ==========
    ErrorRule(
        error_patterns=[r'Tors failed|FormBX had a problem|Bend failed'],
        error_category=ErrorCategory.GEOMETRY_OPTIMIZATION,
        priority=3,
        description="几何优化中的内部坐标错误（Tors failed/FormBX）",
        auto_fix_strategies=[
            RepairStrategy(
                name="对称分子优化策略（第1级）",
                keywords_to_add="Opt=(NoMicro,CalcFC,MaxCycles=300) Integral=UltraFine",
                description="禁用微步优化，计算初始 Hessian，超精细积分"
            ),
            RepairStrategy(
                name="笛卡尔坐标重启（第2级）",
                keywords_to_add="Opt=(Cartesian,MaxStep=5) Integral=UltraFine NoSymm Geom=AllCheck Guess=Read",
                restart_from_chk=True,
                use_cartesian=True,
                description="切换到笛卡尔坐标"
            )
        ],
        max_retries=2,
        requires_manual_check=True
    ),

    # ========== 优先级 4: SCF 不收敛 ==========
    ErrorRule(
        error_patterns=[r'No convergence in SCF|SCF has not converged|Convergence failure in SCF'],
        error_category=ErrorCategory.SCF_CONVERGENCE,
        priority=4,
        description="SCF 迭代未收敛",
        auto_fix_strategies=[
            RepairStrategy(
                name="SCF 收敛策略（第1级）",
                keywords_to_add="SCF=(XQC,MaxCycle=512)",
                description="使用二次收敛算法，增加迭代次数"
            ),
            RepairStrategy(
                name="SCF 收敛策略（第2级）",
                keywords_to_add="SCF=(XQC,MaxCycle=512) Guess=Mix Integral=UltraFine",
                description="添加混合初始猜测，超精细积分"
            )
        ],
        max_retries=2,
        requires_manual_check=True
    ),

    # ========== 优先级 5: SCF 密度矩阵振荡 ==========
    ErrorRule(
        error_patterns=[r'Density matrix convergence failure|SCF is confused'],
        error_category=ErrorCategory.SCF_CONVERGENCE,
        priority=5,
        description="SCF 密度矩阵振荡",
        auto_fix_strategies=[
            RepairStrategy(
                name="阻尼策略（第1级）",
                keywords_to_add="SCF=(Damp,MaxCycle=512)",
                description="使用阻尼算法稳定收敛"
            ),
            RepairStrategy(
                name="阻尼策略（第2级）",
                keywords_to_add="SCF=(Damp,MaxCycle=512) NoSymm",
                description="关闭对称性，进一步稳定"
            )
        ],
        max_retries=2,
        requires_manual_check=True
    ),

    # ========== 优先级 6: 优化步数用完 ==========
    ErrorRule(
        error_patterns=[r'Optimization stopped.*Number of steps exceeded|Max optimization cycles reached'],
        error_category=ErrorCategory.GEOMETRY_OPTIMIZATION,
        priority=6,
        description="几何优化步数用完",
        auto_fix_strategies=[
            RepairStrategy(
                name="增加优化步数（第1级）",
                keywords_to_add="Opt=(MaxCycles=400,CalcFC)",
                description="增加最大优化步数，计算初始 Hessian"
            ),
            RepairStrategy(
                name="增加优化步数（第2级）",
                keywords_to_add="Opt=(MaxCycles=600,CalcFC,VeryTight)",
                description="进一步增加步数，使用更严格的收敛标准"
            )
        ],
        max_retries=2,
        requires_manual_check=True
    ),

    # ========== 优先级 7: 线性相关 ==========
    ErrorRule(
        error_patterns=[r'Linear dependencies in overlap matrix|Small eigenvalue of overlap matrix'],
        error_category=ErrorCategory.INTEGRAL_PRECISION,
        priority=7,
        description="基组线性相关",
        auto_fix_strategies=[
            RepairStrategy(
                name="线性相关处理（第1级）",
                keywords_to_add="SCF=NoVarAcc Pop=Full",
                description="禁用可变精度，输出完整种群分析"
            ),
            RepairStrategy(
                name="线性相关处理（第2级）",
                keywords_to_add="SCF=NoVarAcc Pop=Full Integral=UltraFine",
                description="使用超精细积分"
            )
        ],
        max_retries=2,
        requires_manual_check=True  # 可能需要降级基组
    ),

    # ========== 优先级 8: 积分精度不足 ==========
    ErrorRule(
        error_patterns=[r'Integral accuracy insufficient|Tried to raise integral accuracy'],
        error_category=ErrorCategory.INTEGRAL_PRECISION,
        priority=8,
        description="积分精度不足",
        auto_fix_strategies=[
            RepairStrategy(
                name="提高积分精度（第1级）",
                keywords_to_add="Integral=UltraFine",
                description="使用超精细积分网格"
            ),
            RepairStrategy(
                name="提高积分精度（第2级）",
                keywords_to_add="Integral=SuperFineGrid",
                description="使用超级精细积分网格"
            )
        ],
        max_retries=2,
        requires_manual_check=False
    ),

    # ========== 优先级 9: 内存不足 ==========
    ErrorRule(
        error_patterns=[r'galloc: could not allocate|Not enough memory available'],
        error_category=ErrorCategory.MEMORY_DISK,
        priority=9,
        description="内存不足",
        auto_fix_strategies=[
            RepairStrategy(
                name="调整内存（第1级）",
                keywords_to_add="",  # 需要在外部处理 %Mem
                description="自动降低 %Mem 到节点内存的 50%"
            )
        ],
        max_retries=1,
        requires_manual_check=True  # 需要管理员介入
    ),

    # ========== 优先级 10: 通用崩溃 ==========
    ErrorRule(
        error_patterns=[r'Error termination via Lnk1e|l9999.exe|segmentation violation'],
        error_category=ErrorCategory.INTERNAL_CRASH,
        priority=10,
        description="高斯内部崩溃",
        auto_fix_strategies=[
            RepairStrategy(
                name="通用急救策略（第1级）",
                keywords_to_add="NoSymm Integral=UltraFine Geom=AllCheck Guess=Read",
                restart_from_chk=True,
                description="关闭对称性，超精细积分，从 chk 重启"
            ),
            RepairStrategy(
                name="通用急救策略（第2级）",
                keywords_to_add="Opt=(Cartesian,MaxStep=5,CalcFC) NoSymm Integral=UltraFine Geom=AllCheck Guess=Read",
                restart_from_chk=True,
                use_cartesian=True,
                description="强制笛卡尔坐标，计算初始 Hessian"
            )
        ],
        max_retries=2,
        requires_manual_check=True
    )
]


class QCAutoRepairEngine:
    """QC 自动修复规则引擎"""

    def __init__(self):
        """初始化引擎"""
        # 按优先级排序规则
        self.rules = sorted(ERROR_RULES, key=lambda r: r.priority)
        self.logger = logging.getLogger(__name__)

    def analyze_error(self, log_content: str) -> Optional[Tuple[ErrorRule, int]]:
        """
        分析错误并返回匹配的规则和策略索引

        Args:
            log_content: Gaussian 日志文件内容

        Returns:
            (匹配的规则, 推荐的策略索引) 或 None
        """
        for rule in self.rules:
            # 检查所有错误特征是否都匹配（AND 组合）
            if all(re.search(pattern, log_content, re.IGNORECASE) 
                   for pattern in rule.error_patterns):
                self.logger.info(f"匹配规则: {rule.description}")
                return (rule, 0)  # 返回第一个策略索引

        return None

    def get_repair_strategy(self, rule: ErrorRule, strategy_index: int) -> Optional[RepairStrategy]:
        """获取修复策略"""
        if strategy_index < len(rule.auto_fix_strategies):
            return rule.auto_fix_strategies[strategy_index]
        return None

    def should_continue_retry(self, error_count: int, total_retries: int, max_total: int = 3) -> bool:
        """
        判断是否应该继续重试

        Args:
            error_count: 该错误类型的重试次数
            total_retries: 总重试次数
            max_total: 最大总重试次数

        Returns:
            是否应该继续重试
        """
        return error_count < 2 and total_retries < max_total

