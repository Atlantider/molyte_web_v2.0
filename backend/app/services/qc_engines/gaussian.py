"""
Gaussian QC引擎实现

封装现有的Gaussian计算逻辑
"""
from pathlib import Path
from typing import Tuple
import os

from app.services.qc_engines import QCEngine, QCCalculationInput, QCCalculationResult
from app.tasks import qc_submission, qc_postprocess


class GaussianEngine(QCEngine):
    """Gaussian QC引擎"""
    
    @property
    def name(self) -> str:
        return "gaussian"
    
    def generate_input(
        self, 
        params: QCCalculationInput, 
        work_dir: Path
    ) -> Tuple[Path, str]:
        """
        生成Gaussian输入文件
        
        使用现有的generate_gaussian_input函数
        """
        # 创建临时QCJob对象用于复用现有代码
        # 这是过渡方案,后续可以重构
        from app.models.qc import QCJob
        
        temp_job = QCJob()
        temp_job.molecule_name = params.molecule_name or "molecule"
        temp_job.smiles = params.smiles
        temp_job.charge = params.charge
        temp_job.spin_multiplicity = params.multiplicity
        temp_job.functional = params.functional
        temp_job.basis_set = params.basis_set
        
        # 构建config
        temp_job.config = {
            "solvent_config": {
                "model": params.solvent_model,
                "solvent_name": params.solvent_name,
            },
            "accuracy_level": "custom"
        }
        
        if params.custom_solvent:
            temp_job.config["solvent_config"].update(params.custom_solvent)
        
        if params.coordinates:
            temp_job.config["xyz_content"] = params.coordinates
        
        # 使用现有函数生成输入文件
        gjf_path, safe_name = qc_submission.generate_gaussian_input(
            temp_job, work_dir, None
        )
        
        return Path(gjf_path), safe_name
    
    def parse_output(self, output_file: Path) -> QCCalculationResult:
        """
        解析Gaussian输出文件
        
        使用现有的extract_gaussian_results函数
        """
        # 查找.log文件
        log_file = output_file
        if not log_file.exists():
            # 尝试查找同目录下的.log文件
            work_dir = output_file.parent
            log_files = list(work_dir.glob("*.log"))
            if log_files:
                log_file = log_files[0]
        
        if not log_file.exists():
            return QCCalculationResult(
                energy_au=0.0,
                success=False,
                error_message=f"Gaussian output file not found: {log_file}"
            )
        
        # 使用现有函数解析结果
        results = qc_postprocess.extract_gaussian_results(str(log_file))
        
        if results.get("energy_au") is None:
            return QCCalculationResult(
                energy_au=0.0,
                success=False,
                error_message="Failed to extract energy from Gaussian output"
            )
        
        return QCCalculationResult(
            energy_au=results["energy_au"],
            homo=results.get("homo"),
            lumo=results.get("lumo"),
            success=True
        )
    
    def get_slurm_script(
        self,
        input_file: Path,
        job_name: str,
        partition: str,
        nodes: int,
        cpus: int,
        time_limit: int
    ) -> str:
        """生成Gaussian Slurm脚本"""
        safe_name = input_file.stem
        
        script_content = f"""#!/bin/bash
#SBATCH --job-name={job_name[:64]}
#SBATCH --output=qc_out.log
#SBATCH --error=qc_err.log
#SBATCH --partition={partition}
#SBATCH --nodes={nodes}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --time={time_limit}

# 进入工作目录
cd $SLURM_SUBMIT_DIR

# 设置Gaussian环境
export g16root=/public/software
export GAUSS_SCRDIR=/public/software/g16/scratch
source /public/software/g16/bsd/g16.profile

# 运行Gaussian
ulimit -s unlimited
g16 < "{safe_name}.gjf" > "{safe_name}_out.log" 2>&1

# 转换checkpoint文件
if [ -f "{safe_name}.chk" ]; then
    formchk "{safe_name}.chk" "{safe_name}.fchk"
fi

echo "QC calculation completed"
"""
        
        return script_content
