"""
Reaction Network Worker Integration
Add these methods to the PollingWorker class in polling_worker.py

Insert after the other _process_*_job methods (around line 6300)
"""

def _process_reaction_network_job(self, job: Dict):
    """处理反应网络任务"""
    job_id = job['id']
    self.logger.info(f"开始处理反应网络任务 {job_id}")
    
    try:
        # 0. 检查是否已在处理
        if job_id in self.running_jobs:
            self.logger.info(f"反应网络任务 {job_id} 已在处理中，跳过")
            return
        
        # 1. 更新状态为QUEUED
        self._update_job_status(job_id, 'QUEUED', 'reaction_network')
        
        # 2. 获取任务配置
        config = job.get('config', {})
        job_name = config.get('job_name', f'RN_{job_id}')
        initial_smiles = config.get('initial_smiles', [])
        
        if not initial_smiles:
            raise Exception("初始分子SMILES列表不能为空")
        
        # 环境参数
        temperature = config.get('temperature', 300.0)
        electrode_type = config.get('electrode_type', 'anode')
        voltage = config.get('voltage', 0.1)
        
        # 网络生成参数
        max_generations = config.get('max_generations', 3)
        max_species = config.get('max_species', 50)
        energy_cutoff = config.get('energy_cutoff', 80.0)
        
        # Slurm资源配置
        slurm_partition = config.get('slurm_partition', 'cpu')
        slurm_cpus = config.get('slurm_cpus', 16)
        slurm_time = config.get('slurm_time', 7200)
        
        # 3. 创建工作目录
        rn_work_base = Path(self.config['local'].get('rn_work_base_path', 
                                                      self.config['local']['work_base_path']))
        work_dir = rn_work_base / f"RN-{job_id}-{job_name}"
        work_dir.mkdir(parents=True, exist_ok=True)
        
        # 4. 生成RSNet Python脚本
        script_path = work_dir / "run_rsnet.py"
        self._generate_rsnet_script(
            script_path,
            initial_smiles=initial_smiles,
            temperature=temperature,
            electrode_type=electrode_type,
            voltage=voltage,
            max_generations=max_generations,
            max_species=max_species,
            energy_cutoff=energy_cutoff,
            output_dir=str(work_dir)
        )
        self.logger.info(f"生成RSNet脚本: {script_path}")
        
        # 5. 生成Slurm作业脚本
        job_script = work_dir / "job.sh"
        self._generate_rsnet_slurm_script(
            job_script,
            work_dir=work_dir,
            cpus=slurm_cpus,
            time_limit=slurm_time,
            partition=slurm_partition
        )
        self.logger.info(f"生成Slurm脚本: {job_script}")
        
        # 6. 提交到Slurm
        slurm_result = self._submit_to_slurm(work_dir)
        
        if slurm_result['success']:
            slurm_job_id = slurm_result['slurm_job_id']
            self.logger.info(f"反应网络任务 {job_id} 已提交到Slurm: {slurm_job_id}")
            
            # 7. 更新数据库
            self._update_job_status(
                job_id, 'RUNNING', 'reaction_network',
                slurm_job_id=slurm_job_id,
                work_dir=str(work_dir)
            )
            
            # 8. 添加到运行队列
            self.running_jobs[job_id] = {
                'type': 'reaction_network',
                'slurm_job_id': slurm_job_id,
                'work_dir': str(work_dir),
                'start_time': time.time(),
                'resp_cpu_hours': 0.0
            }
        else:
            raise Exception(f"Slurm提交失败: {slurm_result['error']}")
            
    except Exception as e:
        self.logger.error(f"处理反应网络任务 {job_id} 失败: {e}", exc_info=True)
        self._update_job_status(job_id, 'FAILED', 'reaction_network', error_message=str(e))


def _generate_rsnet_script(self, script_path: Path, **kwargs):
    """生成RSNet计算Python脚本"""
    
    # 处理SMILES列表
    smiles_str = ',\n        '.join([f'"{s}"' for s in kwargs['initial_smiles']])
    
    script_content = f'''#!/usr/bin/env python3
"""
RSNet反应网络生成脚本
自动生成，请勿手动修改
"""

import sys
import json
import os
from pathlib import Path

# 添加rsnet到路径
rsnet_path = Path('/opt/molyte_web_v2.0/rsnet-main/rsnet-main')
if rsnet_path.exists():
    sys.path.insert(0, str(rsnet_path))
else:
    print("Error: RSNet path not found!")
    sys.exit(1)

try:
    from rsnet_simple_api import generate_reaction_network
    
    print("=== 开始反应网络生成 ===")
    print(f"初始分子数: {len(kwargs['initial_smiles'])}")
    print(f"温度: {kwargs['temperature']} K")
    print(f"电极类型: {kwargs['electrode_type']}")
    print(f"最大代数: {kwargs['max_generations']}")
    
    result = generate_reaction_network(
        smiles_list=[
            {smiles_str}
        ],
        temperature={kwargs['temperature']},
        electrode_type="{kwargs['electrode_type']}",
        voltage={kwargs['voltage']},
        max_generations={kwargs['max_generations']},
        max_species={kwargs['max_species']},
        energy_cutoff={kwargs['energy_cutoff']},
        visualize=True,
        save_results=True,
        output_dir="{kwargs['output_dir']}"
    )
    
    print(f"\\n=== 反应网络生成完成 ===")
    print(f"生成分子数: {{len(result['molecules'])}}")
    print(f"发现反应数: {{len(result['reactions'])}}")
    print(f"最大代数: {{result['statistics']['max_generation']}}")
    
    # 保存详细结果
    output_file = Path("{kwargs['output_dir']}") / "network_result.json"
    with open(output_file, "w", encoding='utf-8') as f:
        json.dump({{
            'statistics': result['statistics'],
            'molecules': [{{
                'name': m.name,
                'smiles': m.smiles,
                'generation': getattr(m, 'generation', 0),
                'energy': getattr(m, 'energy', None)
            }} for m in result['molecules']],
            'reactions': [{{
                'reactants': [r.smiles for r in rxn.reactants],
                'products': [p.smiles for p in rxn.products],
                'operator': getattr(rxn, 'operator_name', 'unknown'),
                'energy': getattr(rxn, 'energy', None)
            }} for rxn in result['reactions']],
            'visualization_path': result.get('visualization_path', ''),
            'generation_time': result.get('generation_time', 0)
        }}, f, indent=2, ensure_ascii=False)
    
    print(f"结果已保存到: {{output_file}}")
    print("反应网络生成任务完成！")
    
except Exception as e:
    print(f"Error: {{e}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''
    
    with open(script_path, 'w', encoding='utf-8') as f:
        f.write(script_content)
    
    import os
    os.chmod(script_path, 0o755)


def _generate_rsnet_slurm_script(self, script_path: Path, work_dir: Path, 
                                 cpus: int, time_limit: int, partition: str):
    """生成反应网络的Slurm作业脚本"""
    
    script_content = f'''#!/bin/bash
#SBATCH --job-name=RSNet
#SBATCH --output=rsnet_out.log
#SBATCH --error=rsnet_err.log
#SBATCH --partition={partition}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --time={time_limit}
#SBATCH --mem=32G

echo "=== RSNet反应网络生成任务 ==="
echo "工作目录: {work_dir}"
echo "CPU核心数: {cpus}"
echo "分区: {partition}"
echo "开始时间: $(date)"

# 进入工作目录
cd "{work_dir}"

# 激活Python环境
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate molyte

# 设置环境变量
export PYTHONPATH="/opt/molyte_web_v2.0/rsnet-main/rsnet-main:$PYTHONPATH"
export OMP_NUM_THREADS={cpus}

# 运行RSNet脚本
echo "\\n=== 开始运行RSNet ==="
python run_rsnet.py

exit_code=$?

if [ $exit_code -eq 0 ]; then
    echo "\\n=== RSNet任务成功完成 ==="
else
    echo "\\n=== RSNet任务失败 (退出码: $exit_code) ==="
fi

echo "结束时间: $(date)"
exit $exit_code
'''
    
    with open(script_path, 'w', encoding='utf-8') as f:
        f.write(script_content)
    
    import os
    os.chmod(script_path, 0o755)


def _handle_reaction_network_completion(self, job_id: int, job_info: Dict):
    """处理反应网络任务完成"""
    
    try:
        work_dir = Path(job_info['work_dir'])
        result_file = work_dir / "network_result.json"
        
        if not result_file.exists():
            raise Exception(f"结果文件不存在: {result_file}")
        
        # 读取结果
        with open(result_file, 'r', encoding='utf-8') as f:
            result = json.load(f)
        
        # 保存分子到数据库
        self._save_reaction_network_molecules(job_id, result['molecules'])
        
        # 保存反应到数据库
        self._save_reaction_network_reactions(job_id, result['reactions'])
        
        # 更新任务统计信息
        num_molecules = len(result['molecules'])
        num_reactions = len(result['reactions'])
        max_generation = result['statistics'].get('max_generation', 0)
        
        # 获取CPU核时
        slurm_job_id = job_info.get('slurm_job_id')
        cpu_hours = self._get_job_cpu_hours(slurm_job_id) if slurm_job_id else 0.0
        
        # 更新任务状态
        self._update_reaction_network_job(
            job_id,
            status='COMPLETED',
            num_molecules=num_molecules,
            num_reactions=num_reactions,
            max_generation_reached=max_generation,
            network_json_path=str(result_file),
            visualization_png_path=result.get('visualization_path', ''),
            cpu_hours=cpu_hours
        )
        
        self.logger.info(f"反应网络任务 {job_id} 完成: {num_molecules} 分子, {num_reactions} 反应")
        
    except Exception as e:
        self.logger.error(f"处理反应网络任务 {job_id} 完成时出错: {e}", exc_info=True)
        self._update_job_status(job_id, 'FAILED', 'reaction_network', 
                               error_message=f"结果处理失败: {str(e)}")


def _save_reaction_network_molecules(self, job_id: int, molecules: List[Dict]):
    """保存分子到数据库"""
    # 这里需要通过API调用保存到数据库
    # 简化实现：调用后端API
    try:
        data = {
            'job_id': job_id,
            'molecules': molecules
        }
        response = requests.post(
            f"{self.api_base_url}/reaction-network/jobs/{job_id}/molecules/batch",
            headers=self.api_headers,
            json=data,
            timeout=30
        )
        if response.status_code != 200:
            self.logger.warning(f"保存分子失败: {response.text}")
    except Exception as e:
        self.logger.error(f"保存分子时出错: {e}")


def _save_reaction_network_reactions(self, job_id: int, reactions: List[Dict]):
    """保存反应到数据库"""
    try:
        data = {
            'job_id': job_id,
            'reactions': reactions
        }
        response = requests.post(
            f"{self.api_base_url}/reaction-network/jobs/{job_id}/reactions/batch",
            headers=self.api_headers,
            json=data,
            timeout=30
        )
        if response.status_code != 200:
            self.logger.warning(f"保存反应失败: {response.text}")
    except Exception as e:
        self.logger.error(f"保存反应时出错: {e}")


def _update_reaction_network_job(self, job_id: int, **kwargs):
    """更新反应网络任务"""
    try:
        response = requests.patch(
            f"{self.api_base_url}/reaction-network/jobs/{job_id}",
            headers=self.api_headers,
            json=kwargs,
            timeout=10
        )
        if response.status_code != 200:
            self.logger.warning(f"更新任务 {job_id} 失败: {response.text}")
    except Exception as e:
        self.logger.error(f"更新任务时出错: {e}")
