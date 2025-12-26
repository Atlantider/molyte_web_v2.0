"""
部分结果上传模块

支持上传部分成功的计算结果（容错机制核心）
"""
import logging
from datetime import datetime
from typing import Dict, Any, List, Optional, TYPE_CHECKING
from pathlib import Path

if TYPE_CHECKING:
    from worker.core.config import WorkerConfig


logger = logging.getLogger(__name__)


class PartialUploader:
    """部分结果上传器
    
    核心功能：
    1. 即使部分计算失败，也保存成功的结果
    2. 记录失败详情供后续分析
    3. 标注结果为partial供前端显示
    """
    
    def __init__(self, config: 'WorkerConfig'):
        """初始化部分结果上传器"""
        self.config = config
        self.logger = logging.getLogger(self.__class__.__name__)
        
        from worker.uploaders.cos_client import COSClient
        self.cos = COSClient(config)
    
    def upload_desolvation_results(
        self,
        job_id: int,
        results: Dict[str, Any],
        is_partial: bool = False
    ):
        """
        上传去溶剂化能结果（支持部分结果）
        
        Args:
            job_id: 后处理任务 ID
            results: 计算结果
            is_partial: 是否为部分结果
        """
        cos_prefix = f"postprocess/desolvation/{job_id}/"
        
        # 添加元数据
        result_data = {
            **results,
            'is_partial': is_partial,
            'upload_time': datetime.now().isoformat(),
            'partial_metadata': {
                'total_calculations': results.get('total_ligands', 0),
                'successful_calculations': results.get('calculated_ligands', 0),
                'failed_calculations': len(results.get('failed_details', [])),
                'completion_rate': self._calculate_completion_rate(results)
            } if is_partial else None
        }
        
        # 上传结果 JSON
        success = self.cos.upload_json(
            f"{cos_prefix}results.json",
            result_data
        )
        
        if success:
            self.logger.info(
                f"{'部分' if is_partial else '完整'}结果已上传: {cos_prefix}"
            )
        else:
            self.logger.error(f"结果上传失败: {cos_prefix}")
        
        # 如果有失败详情，单独保存
        if is_partial and results.get('failed_details'):
            self.cos.upload_json(
                f"{cos_prefix}failed_details.json",
                {
                    'job_id': job_id,
                    'failed_calculations': results['failed_details'],
                    'recorded_at': datetime.now().isoformat()
                }
            )
        
        # 上传可视化文件（如果有）
        if results.get('plot_files'):
            for plot in results['plot_files']:
                if isinstance(plot, dict):
                    local_path = Path(plot.get('local_path', ''))
                    if local_path.exists():
                        self.cos.upload_file(
                            local_path,
                            f"{cos_prefix}{plot.get('filename', local_path.name)}"
                        )
    
    def upload_redox_results(
        self,
        job_id: int,
        results: Dict[str, Any],
        is_partial: bool = False
    ):
        """上传氧化还原电位结果"""
        cos_prefix = f"postprocess/redox/{job_id}/"
        
        result_data = {
            **results,
            'is_partial': is_partial,
            'upload_time': datetime.now().isoformat()
        }
        
        self.cos.upload_json(f"{cos_prefix}results.json", result_data)
        
        self.logger.info(
            f"{'部分' if is_partial else '完整'} Redox 结果已上传: {cos_prefix}"
        )
    
    def upload_binding_results(
        self,
        job_id: int,
        results: Dict[str, Any]
    ):
        """上传 Binding 结果"""
        cos_prefix = f"postprocess/binding/{job_id}/"
        
        self.cos.upload_json(
            f"{cos_prefix}results.json",
            {
                **results,
                'upload_time': datetime.now().isoformat()
            }
        )
        
        self.logger.info(f"Binding 结果已上传: {cos_prefix}")
    
    def _calculate_completion_rate(self, results: Dict) -> float:
        """计算完成率"""
        total = results.get('total_ligands', 0)
        calculated = results.get('calculated_ligands', 0)
        
        if total == 0:
            return 0.0
        return round(calculated / total * 100, 2)
