"""
Slurm 工具模块

提供 Slurm 作业管理功能
"""
import logging
import subprocess
import re
from typing import Dict, Any, List, Optional
from pathlib import Path


logger = logging.getLogger(__name__)


class SlurmManager:
    """Slurm 作业管理器"""
    
    def __init__(self):
        """初始化 Slurm 管理器"""
        self.logger = logging.getLogger(self.__class__.__name__)
    
    def submit_job(self, work_dir: Path) -> Dict[str, Any]:
        """
        提交作业到 Slurm
        
        Args:
            work_dir: 工作目录，应包含 job.sh
            
        Returns:
            {'success': bool, 'slurm_job_id': str, 'error': str}
        """
        job_script = work_dir / 'job.sh'
        
        if not job_script.exists():
            return {'success': False, 'error': f'Job script not found: {job_script}'}
        
        try:
            result = subprocess.run(
                ['sbatch', str(job_script)],
                cwd=str(work_dir),
                capture_output=True,
                text=True,
                timeout=30
            )
            
            if result.returncode == 0:
                # 解析 "Submitted batch job 12345"
                match = re.search(r'Submitted batch job (\d+)', result.stdout)
                if match:
                    slurm_job_id = match.group(1)
                    self.logger.info(f"作业提交成功: {slurm_job_id}")
                    return {'success': True, 'slurm_job_id': slurm_job_id}
                else:
                    return {'success': False, 'error': f'Cannot parse job ID: {result.stdout}'}
            else:
                return {'success': False, 'error': result.stderr}
                
        except subprocess.TimeoutExpired:
            return {'success': False, 'error': 'sbatch timeout'}
        except Exception as e:
            return {'success': False, 'error': str(e)}
    
    def get_job_status(self, slurm_job_id: str) -> str:
        """
        获取作业状态
        
        Returns:
            状态字符串: PENDING, RUNNING, COMPLETED, FAILED, CANCELLED, TIMEOUT
        """
        try:
            result = subprocess.run(
                ['sacct', '-j', slurm_job_id, '--format=State', '--noheader', '-P'],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode == 0:
                # 取第一行状态
                lines = result.stdout.strip().split('\n')
                if lines:
                    status = lines[0].strip()
                    # 标准化状态
                    if 'COMPLETED' in status:
                        return 'COMPLETED'
                    elif 'RUNNING' in status:
                        return 'RUNNING'
                    elif 'PENDING' in status:
                        return 'PENDING'
                    elif 'FAILED' in status:
                        return 'FAILED'
                    elif 'CANCELLED' in status:
                        return 'CANCELLED'
                    elif 'TIMEOUT' in status:
                        return 'TIMEOUT'
                    else:
                        return status
            
            return 'UNKNOWN'
            
        except Exception as e:
            self.logger.error(f"获取作业状态失败: {e}")
            return 'UNKNOWN'
    
    def get_job_cpu_hours(self, slurm_job_id: str) -> float:
        """
        获取作业消耗的 CPU 核时
        
        Returns:
            CPU 核时数
        """
        try:
            result = subprocess.run(
                [
                    'sacct', '-j', slurm_job_id,
                    '--format=CPUTimeRAW,AllocCPUS',
                    '--noheader', '-P'
                ],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode == 0:
                lines = result.stdout.strip().split('\n')
                if lines:
                    parts = lines[0].split('|')
                    if len(parts) >= 2:
                        cpu_seconds = int(parts[0])
                        # 转换为核时
                        return cpu_seconds / 3600.0
            
            return 0.0
            
        except Exception as e:
            self.logger.error(f"获取 CPU 核时失败: {e}")
            return 0.0
    
    def cancel_job(self, slurm_job_id: str) -> bool:
        """取消作业"""
        try:
            result = subprocess.run(
                ['scancel', slurm_job_id],
                capture_output=True,
                text=True,
                timeout=10
            )
            return result.returncode == 0
        except Exception as e:
            self.logger.error(f"取消作业失败: {e}")
            return False
    
    def get_partitions(self) -> List[Dict[str, Any]]:
        """获取分区信息"""
        try:
            result = subprocess.run(
                [
                    'sinfo',
                    '--format=%P|%a|%D|%C',
                    '--noheader'
                ],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            partitions = []
            if result.returncode == 0:
                for line in result.stdout.strip().split('\n'):
                    parts = line.split('|')
                    if len(parts) >= 4:
                        name = parts[0].rstrip('*')
                        state = parts[1]
                        nodes = int(parts[2])
                        # CPUS: A/I/O/T (Allocated/Idle/Other/Total)
                        cpu_parts = parts[3].split('/')
                        
                        partitions.append({
                            'name': name,
                            'state': state,
                            'total_nodes': nodes,
                            'available_nodes': nodes,  # 简化
                            'total_cpus': int(cpu_parts[-1]) if len(cpu_parts) >= 4 else 0,
                            'available_cpus': int(cpu_parts[1]) if len(cpu_parts) >= 2 else 0
                        })
            
            return partitions
            
        except Exception as e:
            self.logger.error(f"获取分区信息失败: {e}")
            return []
    
    def get_running_jobs(self) -> Dict[str, str]:
        """
        获取当前运行的作业
        
        Returns:
            {slurm_job_id: job_name}
        """
        try:
            result = subprocess.run(
                [
                    'squeue', '-u', '$USER',
                    '--format=%i|%j',
                    '--noheader'
                ],
                capture_output=True,
                text=True,
                timeout=10,
                shell=True
            )
            
            jobs = {}
            if result.returncode == 0:
                for line in result.stdout.strip().split('\n'):
                    if line:
                        parts = line.split('|')
                        if len(parts) >= 2:
                            jobs[parts[0]] = parts[1]
            
            return jobs
            
        except Exception as e:
            self.logger.error(f"获取运行中作业失败: {e}")
            return {}
    
    def generate_qc_job_script(
        self,
        output_path: Path,
        job_name: str,
        partition: str,
        cpus: int,
        time_limit: int,
        work_dir: Path
    ):
        """生成 QC 作业脚本"""
        gjf_file = list(work_dir.glob('*.gjf'))[0].name if list(work_dir.glob('*.gjf')) else 'input.gjf'
        log_file = gjf_file.replace('.gjf', '.log')
        
        script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --time={time_limit // 60}:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# 加载 Gaussian 环境
source /opt/gaussian/g16/bsd/g16.profile

# 设置工作目录
cd {work_dir}

# 运行 Gaussian
g16 < {gjf_file} > {log_file}

echo "Gaussian calculation completed at $(date)"
"""
        
        with open(output_path, 'w') as f:
            f.write(script)
        
        self.logger.info(f"生成作业脚本: {output_path}")
    
    def generate_md_job_script(
        self,
        output_path: Path,
        job_name: str,
        partition: str,
        cpus: int,
        time_limit: int,
        work_dir: Path,
        input_file: str = 'in.lammps'
    ):
        """生成 MD 作业脚本"""
        script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
#SBATCH --nodes=1
#SBATCH --ntasks={cpus}
#SBATCH --time={time_limit // 60}:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# 加载 LAMMPS 环境
module load lammps/2023

# 设置工作目录
cd {work_dir}

# 运行 LAMMPS
mpirun -np {cpus} lmp -in {input_file}

echo "LAMMPS simulation completed at $(date)"
"""
        
        with open(output_path, 'w') as f:
            f.write(script)
        
        self.logger.info(f"生成 MD 作业脚本: {output_path}")
