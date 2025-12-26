"""
Quota checking and estimation utilities

提供配额检查和核时预估功能
"""
from typing import Dict, Any
from app.models.user import User
from sqlalchemy.orm import Session
import logging

logger = logging.getLogger(__name__)


def check_submission_quota(
    user: User,
    db: Session,
    estimated_hours: float,
    min_balance_required: float = 1.0
) -> Dict[str, Any]:
    """
    检查用户是否有足够配额提交任务
    
    Args:
        user: 用户对象
        db: 数据库会话
        estimated_hours: 预估任务核时消耗
        min_balance_required: 最小余额要求(默认1核时)
    
    Returns:
        {
            'can_submit': bool,
            'reason': str,  # 不能提交的原因
            'available_quota': float,
            'required_quota': float,
            'debt': float  # 欠费金额(如果有)
        }
    """
    from app.services.quota_service import QuotaService
    
    # 获取可用配额
    available_quota = QuotaService.get_available_quota(user, db)
    
    # 检查1: 是否有欠费
    if user.balance_cpu_hours < 0:
        debt = abs(user.balance_cpu_hours)
        return {
            'can_submit': False,
            'reason': f'账户欠费{debt:.2f}核时，请先充值',
            'available_quota': available_quota,
            'required_quota': estimated_hours,
            'debt': debt
        }
    
    # 检查2: 可用配额是否足够(预估核时 + 最小余额)
    required_quota = estimated_hours + min_balance_required
    if available_quota < required_quota:
        return {
            'can_submit': False,
            'reason': f'可用配额不足。需要{required_quota:.2f}核时，当前可用{available_quota:.2f}核时',
            'available_quota': available_quota,
            'required_quota': required_quota,
            'debt': 0.0
        }
    
    # 检查通过
    return {
        'can_submit': True,
        'reason': '',
        'available_quota': available_quota,
        'required_quota': required_quota,
        'debt': 0.0
    }


def estimate_qc_job_hours(
    basis_set: str,
    functional: str,
    atom_count: int = 10,
    calculation_type: str = 'OPT'  # SP, OPT, FREQ
) -> float:
    """
    预估QC任务核时消耗
    
    基于经验公式估算
    
    Args:
        basis_set: 基组
        functional: 泛函
        atom_count: 原子数
        calculation_type: 计算类型
    
    Returns:
        预估核时数(小时)
    """
    # 基础时间(单点能量,小时)
    base_time = {
        'STO-3G': 0.01,
        '3-21G': 0.02,
        '6-31G': 0.05,
        '6-31G(d)': 0.1,
        '6-31+G(d)': 0.15,
        '6-31++G(d,p)': 0.2,
        '6-311G(d,p)': 0.3,
        '6-311++G(d,p)': 0.4,
        '6-311++G(2d,2p)': 0.5,
    }.get(basis_set, 0.2)
    
    # 泛函系数
    functional_factor = {
        'HF': 0.5,
        'B3LYP': 1.0,
        'wB97XD': 1.5,
        'M06-2X': 1.3,
        'PBE': 0.9,
    }.get(functional, 1.0)
    
    # 原子数系数(非线性,N^2复杂度)
    atom_factor = (atom_count / 10) ** 2
    
    # 计算类型系数
    calc_type_factor = {
        'SP': 1.0,      # 单点能量
        'OPT': 10.0,    # 几何优化
        'FREQ': 15.0,   # 频率计算
        'RESP': 5.0,    # RESP电荷
    }.get(calculation_type, 1.0)
    
    # CPU数(假设16核)
    cpu_count = 16
    
    # 估算核时 = 基础时间 × 泛函系数 × 原子系数 × 计算类型系数 × CPU数
    estimated_hours = base_time * functional_factor * atom_factor * calc_type_factor * cpu_count
    
    # 添加20%安全余量
    estimated_hours *= 1.2
    
    # 限制范围: 最小0.5核时,最大200核时
    return max(0.5, min(200.0, estimated_hours))


def estimate_md_job_hours(
    steps: int,
    system_size: int = 1000,
    ensemble: str = 'NPT'
) -> float:
    """
    预估MD任务核时消耗
    
    Args:
        steps: MD步数
        system_size: 系统大小(原子数)
        ensemble: 系综类型
    
    Returns:
        预估核时数(小时)
    """
    # MD步数系数(每1000步的单核时间)
    time_per_1000steps = 0.05  # 1000步约0.05小时(单核)
    
    # 系统大小系数(非线性)
    size_factor = (system_size / 1000) ** 1.5
    
    # 系综类型系数
    ensemble_factor = {
        'NVE': 0.8,
        'NVT': 1.0,
        'NPT': 1.2,
    }.get(ensemble, 1.0)
    
    # CPU数(MD通常用更多CPU)
    cpu_count = 32
    
    # 估算核时
    estimated_hours = (steps / 1000) * time_per_1000steps * size_factor * ensemble_factor * cpu_count
    
    # 安全余量30%
    estimated_hours *= 1.3
    
    # 限制范围
    return max(1.0, min(500.0, estimated_hours))
