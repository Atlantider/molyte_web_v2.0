"""
QC计算引擎抽象基类

提供统一的接口支持多种QC计算引擎(Gaussian, PySCF等)
"""
from abc import ABC, abstractmethod
from typing import Dict, Any, Tuple, List, Optional
from pathlib import Path
from pydantic import BaseModel


class QCCalculationInput(BaseModel):
    """QC计算输入参数"""
    smiles: str
    charge: int
    multiplicity: int
    functional: str
    basis_set: str
    solvent_model: str  # gas, pcm, smd, custom
    solvent_name: Optional[str] = None
    custom_solvent: Optional[Dict] = None
    calc_type: str = "opt freq"  # opt, freq, opt freq, sp
    coordinates: Optional[str] = None  # XYZ坐标(可选)
    molecule_name: Optional[str] = None


class QCCalculationResult(BaseModel):
    """QC计算结果"""
    energy_au: float
    homo: Optional[float] = None
    lumo: Optional[float] = None
    dipole: Optional[float] = None
    optimized_coords: Optional[str] = None
    frequencies: Optional[List[float]] = None
    success: bool = True
    error_message: Optional[str] = None


class QCEngine(ABC):
    """QC计算引擎抽象基类"""
    
    @property
    @abstractmethod
    def name(self) -> str:
        """引擎名称"""
        pass
    
    @abstractmethod
    def generate_input(
        self, 
        params: QCCalculationInput, 
        work_dir: Path
    ) -> Tuple[Path, str]:
        """
        生成输入文件
        
        Args:
            params: 计算参数
            work_dir: 工作目录
            
        Returns:
            (input_file_path, job_name)
        """
        pass
    
    @abstractmethod
    def parse_output(
        self,
        output_file: Path
    ) -> QCCalculationResult:
        """
        解析输出文件
        
        Args:
            output_file: 输出文件路径
            
        Returns:
            计算结果
        """
        pass
    
    @abstractmethod
    def get_slurm_script(
        self,
        input_file: Path,
        job_name: str,
        partition: str,
        nodes: int,
        cpus: int,
        time_limit: int
    ) -> str:
        """
        生成Slurm提交脚本
        
        Args:
            input_file: 输入文件路径
            job_name: 任务名称
            partition: 分区名称
            nodes: 节点数
            cpus: CPU核数
            time_limit: 时间限制(分钟)
            
        Returns:
            Slurm脚本内容
        """
        pass
