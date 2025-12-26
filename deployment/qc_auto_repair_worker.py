#!/usr/bin/env python3
"""
QC 自动修复 Worker 模块

这个模块实现了完整的自动化重试流程：
1. 检测失败的 QC 任务
2. 分析错误日志
3. 应用修复策略
4. 生成新的 GJF 文件
5. 自动提交到 Slurm（使用原始配置）
6. 监控重试任务
7. 自动后处理结果
8. 返回结果到前端

集成到 Polling Worker 中使用。
"""

import re
import logging
import subprocess
import shutil
import time
from pathlib import Path
from typing import Dict, Optional, Any, Callable
from dataclasses import dataclass

from qc_auto_repair_engine import QCAutoRepairEngine


@dataclass
class AutoRepairResult:
    """自动修复结果"""
    success: bool
    error_message: str = ""
    matched_rule_description: str = ""
    applied_strategy_name: str = ""
    new_gjf_path: Optional[Path] = None
    slurm_job_id: Optional[str] = None
    retry_count: int = 0
    should_continue_retry: bool = False


class QCAutoRepairWorker:
    """QC 自动修复 Worker"""

    def __init__(self, logger=None, post_process_callback: Optional[Callable] = None):
        """
        初始化

        Args:
            logger: 日志记录器
            post_process_callback: 后处理回调函数，签名为 (job_id, work_dir) -> bool
                                  返回 True 表示后处理成功，False 表示失败
        """
        self.engine = QCAutoRepairEngine()
        self.logger = logger or logging.getLogger(__name__)
        self.post_process_callback = post_process_callback
    
    def auto_repair_and_resubmit(self, job_id: int, work_dir: Path, 
                                 retry_count: int = 0) -> AutoRepairResult:
        """
        自动修复失败的 QC 任务并重新提交到 Slurm
        
        Args:
            job_id: 任务 ID
            work_dir: 工作目录
            retry_count: 当前重试次数
        
        Returns:
            AutoRepairResult 修复结果
        """
        try:
            # 1. 查找日志文件
            log_files = list(work_dir.glob("*_out.log")) + list(work_dir.glob("*.log"))
            log_files = [f for f in log_files if f.name not in ['qc_out.log', 'qc_err.log']]
            
            if not log_files:
                return AutoRepairResult(
                    success=False,
                    error_message="未找到 Gaussian 日志文件"
                )
            
            log_file = log_files[0]
            log_content = log_file.read_text(errors='ignore')
            
            # 2. 分析错误
            result = self.engine.analyze_error(log_content)
            if not result:
                return AutoRepairResult(
                    success=False,
                    error_message="未匹配到任何修复规则"
                )
            
            rule, strategy_idx = result
            strategy = rule.auto_fix_strategies[strategy_idx]
            
            self.logger.info(f"任务 {job_id}: 匹配规则 '{rule.description}'")
            self.logger.info(f"任务 {job_id}: 应用策略 '{strategy.name}'")
            
            # 3. 检查是否应该继续重试
            should_continue = self.engine.should_continue_retry(retry_count, rule.max_retries)
            if not should_continue:
                return AutoRepairResult(
                    success=False,
                    error_message=f"已达到最大重试次数 ({rule.max_retries})",
                    matched_rule_description=rule.description,
                    applied_strategy_name=strategy.name,
                    retry_count=retry_count,
                    should_continue_retry=False
                )
            
            # 4. 查找 GJF 文件
            gjf_files = list(work_dir.glob("*.gjf"))
            if not gjf_files:
                return AutoRepairResult(
                    success=False,
                    error_message="未找到 GJF 文件"
                )
            
            gjf_file = gjf_files[0]
            
            # 5. 生成修复后的 GJF
            new_gjf_content = self._generate_repaired_gjf(gjf_file, strategy)
            
            # 6. 备份原始 GJF
            backup_gjf = work_dir / f"{gjf_file.stem}_retry{retry_count}.gjf"
            if not backup_gjf.exists():
                shutil.copy(gjf_file, backup_gjf)
            
            # 7. 写入新的 GJF
            gjf_file.write_text(new_gjf_content)
            self.logger.info(f"任务 {job_id}: 已生成修复后的 GJF 文件")
            
            # 8. 提交到 Slurm
            job_script = work_dir / "job.sh"
            if not job_script.exists():
                return AutoRepairResult(
                    success=False,
                    error_message="未找到 job.sh 脚本",
                    matched_rule_description=rule.description,
                    applied_strategy_name=strategy.name,
                    new_gjf_path=gjf_file,
                    retry_count=retry_count + 1,
                    should_continue_retry=True
                )
            
            slurm_result = self._submit_to_slurm(work_dir, job_script, strategy)
            if not slurm_result['success']:
                return AutoRepairResult(
                    success=False,
                    error_message=f"Slurm 提交失败: {slurm_result.get('error')}",
                    matched_rule_description=rule.description,
                    applied_strategy_name=strategy.name,
                    new_gjf_path=gjf_file,
                    retry_count=retry_count + 1,
                    should_continue_retry=True
                )
            
            slurm_job_id = slurm_result['slurm_job_id']
            self.logger.info(f"任务 {job_id}: 已提交到 Slurm (Job ID: {slurm_job_id})")
            
            return AutoRepairResult(
                success=True,
                matched_rule_description=rule.description,
                applied_strategy_name=strategy.name,
                new_gjf_path=gjf_file,
                slurm_job_id=slurm_job_id,
                retry_count=retry_count + 1,
                should_continue_retry=True
            )
        
        except Exception as e:
            self.logger.error(f"任务 {job_id}: 自动修复失败 - {e}", exc_info=True)
            return AutoRepairResult(
                success=False,
                error_message=str(e)
            )
    
    def _generate_repaired_gjf(self, gjf_path: Path, strategy) -> str:
        """生成修复后的 GJF 文件内容"""
        content = gjf_path.read_text()
        lines = content.split('\n')
        
        # 找到 route 行
        route_idx = -1
        for i, line in enumerate(lines):
            if line.strip().startswith('#'):
                route_idx = i
                break
        
        if route_idx == -1:
            raise ValueError("未找到 route 行")
        
        original_route = lines[route_idx].strip()
        new_route = self._modify_route_line(original_route, strategy)
        lines[route_idx] = new_route
        
        return '\n'.join(lines)
    
    def _modify_route_line(self, original_route: str, strategy) -> str:
        """修改 route 行"""
        modified = original_route
        
        # 移除可能冲突的关键词
        modified = re.sub(r'\s*Opt=\([^)]*\)', '', modified)
        modified = re.sub(r'\s+NoSymm\b', '', modified)
        modified = re.sub(r'\s+Integral=\w+', '', modified)
        modified = re.sub(r'\s+Geom=\w+', '', modified)
        modified = re.sub(r'\s+Guess=\w+', '', modified)
        modified = re.sub(r'\s+SCF=\w+', '', modified)
        
        # 清理多余空格
        modified = re.sub(r'\s+', ' ', modified).strip()
        
        # 添加新的关键词
        new_route = modified.rstrip() + ' ' + strategy.keywords_to_add
        new_route = re.sub(r'\s+', ' ', new_route).strip()
        
        return new_route
    
    def _submit_to_slurm(self, work_dir: Path, job_script: Path, strategy=None) -> Dict[str, Any]:
        """
        提交任务到 Slurm

        如果修复策略中指定了 Slurm 参数修改，则修改 job.sh 中的相应参数。
        否则使用原始的 job.sh 配置。
        """
        try:
            # 如果策略中指定了 Slurm 参数修改，则修改 job.sh
            if strategy and (strategy.slurm_memory_gb or strategy.slurm_time_seconds or strategy.slurm_cpus):
                job_script = self._modify_slurm_parameters(work_dir, job_script, strategy)

            result = subprocess.run(
                ['sbatch', str(job_script)],
                cwd=str(work_dir),
                capture_output=True,
                text=True,
                timeout=60
            )

            if result.returncode != 0:
                return {'success': False, 'error': f'sbatch failed: {result.stderr}'}

            # 解析 Slurm job ID
            match = re.search(r'Submitted batch job (\d+)', result.stdout)
            if match:
                return {'success': True, 'slurm_job_id': match.group(1)}
            else:
                return {'success': False, 'error': f'Could not parse Slurm job ID'}

        except Exception as e:
            return {'success': False, 'error': str(e)}

    def _modify_slurm_parameters(self, work_dir: Path, job_script: Path, strategy) -> Path:
        """
        修改 job.sh 中的 Slurm 参数

        如果策略中指定了参数修改，则创建一个临时的修改后的 job.sh。
        否则返回原始的 job.sh。
        """
        job_script_content = job_script.read_text()
        modified = False

        # 修改内存
        if strategy.slurm_memory_gb:
            old_mem = re.search(r'#SBATCH\s+--mem=\d+G', job_script_content)
            if old_mem:
                job_script_content = job_script_content.replace(
                    old_mem.group(0),
                    f'#SBATCH --mem={strategy.slurm_memory_gb}G'
                )
                self.logger.info(f"修改内存: {old_mem.group(0)} -> #SBATCH --mem={strategy.slurm_memory_gb}G")
                modified = True

        # 修改时间限制
        if strategy.slurm_time_seconds:
            old_time = re.search(r'#SBATCH\s+--time=\d+', job_script_content)
            if old_time:
                job_script_content = job_script_content.replace(
                    old_time.group(0),
                    f'#SBATCH --time={strategy.slurm_time_seconds}'
                )
                self.logger.info(f"修改时间: {old_time.group(0)} -> #SBATCH --time={strategy.slurm_time_seconds}")
                modified = True

        # 修改 CPU 数量
        if strategy.slurm_cpus:
            old_cpus = re.search(r'#SBATCH\s+--cpus-per-task=\d+', job_script_content)
            if old_cpus:
                job_script_content = job_script_content.replace(
                    old_cpus.group(0),
                    f'#SBATCH --cpus-per-task={strategy.slurm_cpus}'
                )
                self.logger.info(f"修改 CPU: {old_cpus.group(0)} -> #SBATCH --cpus-per-task={strategy.slurm_cpus}")
                modified = True

        # 如果有修改，写入临时文件
        if modified:
            temp_job_script = work_dir / "job_modified.sh"
            temp_job_script.write_text(job_script_content)
            self.logger.info(f"已创建修改后的 job 脚本: {temp_job_script}")
            return temp_job_script

        # 没有修改，返回原始脚本
        return job_script

    def monitor_and_postprocess_retry(self, job_id: int, work_dir: Path,
                                      slurm_job_id: str, max_wait_seconds: int = 3600) -> Dict[str, Any]:
        """
        监控重试任务并自动后处理结果

        这个方法会：
        1. 监控 Slurm 任务状态
        2. 等待任务完成
        3. 检查计算是否成功
        4. 自动调用后处理回调函数
        5. 返回最终结果

        Args:
            job_id: 任务 ID
            work_dir: 工作目录
            slurm_job_id: Slurm Job ID
            max_wait_seconds: 最大等待时间（秒）

        Returns:
            {
                'success': bool,
                'status': 'COMPLETED' | 'FAILED' | 'TIMEOUT',
                'error_message': str,
                'postprocess_success': bool,
                'postprocess_error': str
            }
        """
        try:
            self.logger.info(f"任务 {job_id}: 开始监控重试任务 (Slurm Job: {slurm_job_id})")

            start_time = time.time()
            poll_interval = 30  # 每 30 秒检查一次

            while time.time() - start_time < max_wait_seconds:
                # 检查 Slurm 任务状态
                status = self._check_slurm_job_status(slurm_job_id)

                if status == 'COMPLETED':
                    self.logger.info(f"任务 {job_id}: Slurm 任务已完成，检查计算结果")

                    # 检查计算是否成功
                    if self._check_gaussian_success(work_dir):
                        self.logger.info(f"任务 {job_id}: Gaussian 计算成功")

                        # 调用后处理回调函数
                        postprocess_success = False
                        postprocess_error = ""

                        if self.post_process_callback:
                            try:
                                self.logger.info(f"任务 {job_id}: 开始后处理结果")
                                postprocess_success = self.post_process_callback(job_id, work_dir)

                                if postprocess_success:
                                    self.logger.info(f"任务 {job_id}: 后处理成功，结果已返回前端")
                                else:
                                    postprocess_error = "后处理回调返回 False"
                                    self.logger.error(f"任务 {job_id}: 后处理失败 - {postprocess_error}")
                            except Exception as e:
                                postprocess_error = str(e)
                                self.logger.error(f"任务 {job_id}: 后处理异常 - {e}", exc_info=True)
                        else:
                            self.logger.warning(f"任务 {job_id}: 未设置后处理回调函数")

                        return {
                            'success': postprocess_success,
                            'status': 'COMPLETED',
                            'error_message': "",
                            'postprocess_success': postprocess_success,
                            'postprocess_error': postprocess_error
                        }
                    else:
                        error_msg = "Gaussian 计算失败或未正常完成"
                        self.logger.error(f"任务 {job_id}: {error_msg}")
                        return {
                            'success': False,
                            'status': 'FAILED',
                            'error_message': error_msg,
                            'postprocess_success': False,
                            'postprocess_error': error_msg
                        }

                elif status == 'FAILED':
                    error_msg = "Slurm 任务失败"
                    self.logger.error(f"任务 {job_id}: {error_msg}")
                    return {
                        'success': False,
                        'status': 'FAILED',
                        'error_message': error_msg,
                        'postprocess_success': False,
                        'postprocess_error': error_msg
                    }

                elif status == 'RUNNING':
                    elapsed = time.time() - start_time
                    self.logger.info(f"任务 {job_id}: 任务运行中... ({elapsed:.0f}s/{max_wait_seconds}s)")

                # 等待后再检查
                time.sleep(poll_interval)

            # 超时
            error_msg = f"任务监控超时 ({max_wait_seconds}s)"
            self.logger.error(f"任务 {job_id}: {error_msg}")
            return {
                'success': False,
                'status': 'TIMEOUT',
                'error_message': error_msg,
                'postprocess_success': False,
                'postprocess_error': error_msg
            }

        except Exception as e:
            error_msg = f"监控任务时发生异常: {str(e)}"
            self.logger.error(f"任务 {job_id}: {error_msg}", exc_info=True)
            return {
                'success': False,
                'status': 'FAILED',
                'error_message': error_msg,
                'postprocess_success': False,
                'postprocess_error': error_msg
            }

    def _check_slurm_job_status(self, slurm_job_id: str) -> str:
        """检查 Slurm 任务状态"""
        try:
            result = subprocess.run(
                ['squeue', '-j', slurm_job_id, '-h', '-o', '%T'],
                capture_output=True,
                text=True,
                timeout=10
            )

            if result.returncode == 0 and result.stdout.strip():
                status = result.stdout.strip().upper()
                # 映射 Slurm 状态到标准状态
                if status in ['RUNNING', 'PENDING', 'CONFIGURING']:
                    return 'RUNNING'
                elif status in ['COMPLETED', 'COMPLETING']:
                    return 'COMPLETED'
                elif status in ['FAILED', 'CANCELLED', 'TIMEOUT', 'NODE_FAIL']:
                    return 'FAILED'
                else:
                    return status
            else:
                # 任务不在队列中，可能已完成
                return 'COMPLETED'

        except Exception as e:
            self.logger.warning(f"检查 Slurm 任务状态失败: {e}")
            return 'UNKNOWN'

    def _check_gaussian_success(self, work_dir: Path) -> bool:
        """检查 Gaussian 计算是否成功"""
        try:
            # 查找日志文件
            log_files = list(work_dir.glob("*_out.log")) + list(work_dir.glob("*.log"))
            log_files = [f for f in log_files if f.name not in ['qc_out.log', 'qc_err.log']]

            if not log_files:
                self.logger.warning(f"未找到 Gaussian 日志文件")
                return False

            log_content = log_files[0].read_text(errors='ignore')

            # 检查是否有成功标记
            if 'Normal termination' in log_content:
                return True

            # 检查是否有失败标记
            if 'Error termination' in log_content or 'segmentation violation' in log_content:
                return False

            # 无法确定
            self.logger.warning(f"无法确定计算状态")
            return False

        except Exception as e:
            self.logger.error(f"检查 Gaussian 成功状态失败: {e}")
            return False

