"""
XTB 预优化模块

在生成 Gaussian 输入文件前，使用 XTB 对 cluster minus 结构进行快速优化。
这可以：
1. 修复 cluster 减去 ligand 导致的结构失衡
2. 减少 Gaussian 的优化步数
3. 加快 Gaussian 收敛
4. 降低计算成本
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any

logger = logging.getLogger(__name__)


def preoptimize_cluster_minus_with_xtb(
    xyz_content: str,
    work_dir: Path,
    charge: int = 1,
    multiplicity: int = 1,
    enable_xtb: bool = True,
    xtb_path: str = "xtb"
) -> Dict[str, Any]:
    """
    使用 XTB 预优化 cluster minus 结构
    
    Args:
        xyz_content: XYZ 格式的 cluster minus 结构
        work_dir: 工作目录
        charge: 分子电荷
        multiplicity: 自旋多重度
        enable_xtb: 是否启用 XTB 优化
        xtb_path: XTB 可执行文件路径
    
    Returns:
        字典，包含：
        - xyz_content: 优化后的 XYZ 内容（如果优化成功）或原始内容
        - optimized: 是否进行了优化
        - energy: XTB 能量（如果优化成功）
        - iterations: 优化迭代次数
        - rmsd: 原子位移 RMSD
        - error: 错误信息（如果有）
    """
    result = {
        'xyz_content': xyz_content,
        'optimized': False,
        'energy': None,
        'iterations': 0,
        'rmsd': None,
        'error': None
    }
    
    if not enable_xtb:
        logger.info("XTB 预优化已禁用，使用原始结构")
        return result
    
    try:
        from app.services.xtb_optimizer import XTBOptimizer
        
        logger.info("开始使用 XTB 预优化 cluster minus 结构...")
        
        optimizer = XTBOptimizer(xtb_path=xtb_path, charge=charge, multiplicity=multiplicity)
        
        xtb_result = optimizer.optimize_xyz(
            xyz_content,
            work_dir,
            max_iterations=500,
            convergence="normal"
        )
        
        if xtb_result:
            result['xyz_content'] = xtb_result['optimized_xyz']
            result['optimized'] = True
            result['energy'] = xtb_result['energy']
            result['iterations'] = xtb_result['iterations']
            result['rmsd'] = xtb_result['rmsd']
            
            logger.info(
                f"✅ XTB 预优化成功: "
                f"能量={xtb_result['energy']:.6f} Hartree, "
                f"迭代={xtb_result['iterations']}, "
                f"RMSD={xtb_result['rmsd']:.6f} Å"
            )
        else:
            result['error'] = "XTB 优化失败"
            logger.warning("❌ XTB 优化失败，使用原始结构")
    
    except ImportError:
        result['error'] = "无法导入 XTBOptimizer"
        logger.warning("❌ 无法导入 XTBOptimizer，使用原始结构")
    except Exception as e:
        result['error'] = str(e)
        logger.error(f"❌ XTB 预优化异常: {e}", exc_info=True)
    
    return result


def should_use_xtb_preoptimization(config: Dict[str, Any]) -> bool:
    """
    根据配置判断是否应该使用 XTB 预优化
    
    Args:
        config: 任务配置
    
    Returns:
        是否应该使用 XTB 预优化
    """
    # 从配置中读取设置，默认启用
    return config.get('use_xtb_preoptimization', True)


def get_xtb_config(config: Dict[str, Any]) -> Dict[str, Any]:
    """
    从任务配置中获取 XTB 相关设置
    
    Args:
        config: 任务配置
    
    Returns:
        XTB 配置字典
    """
    xtb_config = config.get('xtb_config', {})
    
    return {
        'enable_xtb': xtb_config.get('enable', True),
        'xtb_path': xtb_config.get('path', 'xtb'),
        'max_iterations': xtb_config.get('max_iterations', 500),
        'convergence': xtb_config.get('convergence', 'normal'),
    }

