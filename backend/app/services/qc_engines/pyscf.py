"""
PySCF QC引擎实现

开源量子化学计算引擎

增强功能:
- Gaussian参数兼容性转换
- XTB结构预优化集成
- 智能SCF收敛参数
- 结果验证和错误处理
"""
from pathlib import Path
from typing import Tuple, Optional, Dict, Any
import json
import logging

from app.services.qc_engines import QCEngine, QCCalculationInput, QCCalculationResult
from app.utils.qc_parameter_converter import (
    convert_basis_set,
    convert_functional,
    get_solvent_dielectric,
    get_scf_convergence_params,
    should_use_xtb_preopt
)

logger = logging.getLogger(__name__)


class PySCFEngine(QCEngine):
    """PySCF QC引擎"""
    
    @property
    def name(self) -> str:
        return "pyscf"
    
    def generate_input(
        self, 
        params: QCCalculationInput, 
        work_dir: Path
    ) -> Tuple[Path, str]:
        """生成PySCF Python脚本"""
        from app.tasks.qc_submission import sanitize_filename
        
        safe_name = sanitize_filename(params.molecule_name or "molecule")
        script_path = work_dir / f"{safe_name}_pyscf.py"
        
        # 生成PySCF计算脚本
        script_content = self._generate_pyscf_script(params, safe_name)
        script_path.write_text(script_content)
        
        return script_path, safe_name
    
    def _generate_pyscf_script(self, params: QCCalculationInput, safe_name: str) -> str:
        """生成PySCF计算脚本内容(增强版)"""
        
        # 参数转换(Gaussian → PySCF)
        pyscf_basis = convert_basis_set(params.basis_set)
        pyscf_functional = convert_functional(params.functional)
        
        # 获取分子坐标(含XTB预优化)
        coords_block = self._get_optimized_coordinates(params)
        
        # SCF收敛参数
        scf_params = self._get_scf_params(params)
        
        # 溶剂模型配置
        solvent_setup = self._get_solvent_setup(params)
        
        return f"""#!/usr/bin/env python
# PySCF QC Calculation Script
# Molecule: {params.molecule_name}
# SMILES: {params.smiles}

import json
import sys
from pyscf import gto, dft, scf
import numpy as np

try:
    # 分子定义
    mol = gto.M(
        atom='''{coords_block}''',
        basis='{pyscf_basis}',
        charge={params.charge},
        spin={params.multiplicity - 1}
    )
    
    # DFT计算
    if {params.multiplicity} == 1:
        mf = dft.RKS(mol)
    else:
        mf = dft.UKS(mol)
    mf.xc = '{pyscf_functional}'
    
    # SCF收敛参数
    {scf_params}
    
    {solvent_setup}
    
    # 能量计算
    energy = mf.kernel()
    
    if not mf.converged:
        raise RuntimeError("SCF not converged")
    
    # 轨道能量
    mo_energy = mf.mo_energy
    if isinstance(mo_energy, tuple):  # UKS
        homo_idx = mol.nelectron[0] - 1
        homo = float(mo_energy[0][homo_idx] * 27.2114) if homo_idx >= 0 else None
        lumo = float(mo_energy[0][homo_idx + 1] * 27.2114) if homo_idx + 1 < len(mo_energy[0]) else None
    else:  # RKS
        homo_idx = mol.nelectron // 2 - 1
        homo = float(mo_energy[homo_idx] * 27.2114) if homo_idx >= 0 else None
        lumo = float(mo_energy[homo_idx + 1] * 27.2114) if homo_idx + 1 < len(mo_energy) else None
    
    # 偶极矩
    dip = mf.dip_moment()
    dipole = float(np.linalg.norm(dip))
    
    # 输出结果(JSON格式)
    results = {{
        'energy_au': float(energy),
        'homo': homo,
        'lumo': lumo,
        'dipole': dipole,
        'success': True,
        'converged': mf.converged
    }}
    
    with open('pyscf_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print("PySCF calculation completed successfully")
    print(f"Energy: {{energy}} Hartree")
    if homo is not None:
        print(f"HOMO: {{homo}} eV")
    if lumo is not None:
        print(f"LUMO: {{lumo}} eV")
    
    sys.exit(0)

except Exception as e:
    error_results = {{
        'energy_au': 0.0,
        'success': False,
        'error_message': str(e)
    }}
    with open('pyscf_results.json', 'w') as f:
        json.dump(error_results, f, indent=2)
    print(f"PySCF calculation failed: {{e}}")
    sys.exit(1)
"""
    
    
    def _get_optimized_coordinates(self, params: QCCalculationInput) -> str:
        """
        获取或生成分子坐标(含XTB预优化)
        
        工作流:
        1. 如果有XYZ坐标 → XTB优化(如适用)
        2. 如果只有SMILES → RDKit生成 → XTB优化
        3. 失败抛出异常
        """
        xyz_content = None
        
        # 1. 获取初始坐标
        if params.coordinates:
            xyz_content = params.coordinates
            logger.debug("Using provided XYZ coordinates")
        else:
            # 从SMILES生成
            xyz_content = self._smiles_to_xyz(params.smiles, params)
            logger.debug("Generated XYZ from SMILES")
        
        # 2. XTB预优化(默认启用)
        use_xtb = getattr(params, 'use_xtb_preopt', True)
        if use_xtb and self._should_use_xtb(params):
            try:
                optimized = self._xtb_optimize_coords(xyz_content, params)
                if optimized:
                    logger.info("XTB pre-optimization successful")
                    return optimized
                else:
                    logger.warning("XTB optimization failed, using original coords")
            except Exception as e:
                logger.warning(f"XTB optimization error: {e}, using original coords")
        
        # 3. 转换为PySCF格式
        return self._xyz_to_pyscf_format(xyz_content)
    
    def _should_use_xtb(self, params) -> bool:
        """判断是否应该使用XTB预优化"""
        # 估算原子数
        if params.coordinates:
            num_atoms = len([l for l in params.coordinates.split('\n') if l.strip() and not l.startswith('#')])
        else:
            # 从SMILES估算
            num_atoms = len(params.smiles) // 2  # 粗略估计
        
        # 大分子或带电体系建议使用XTB
        return num_atoms > 10 or abs(params.charge) > 0
    
    def _xtb_optimize_coords(self, xyz_content: str, params) -> Optional[str]:
        """使用XTB优化坐标"""
        try:
            from app.services.xtb_optimizer import XTBOptimizer
            from pathlib import Path
            import tempfile
            
            with tempfile.TemporaryDirectory() as tmpdir:
                optimizer = XTBOptimizer(
                    charge=params.charge,
                    multiplicity=params.multiplicity
                )
                
                result = optimizer.optimize_xyz(
                    xyz_content,
                    Path(tmpdir),
                    convergence='normal'
                )
                
                if result and result.get('converged'):
                    logger.info(f"XTB optimization converged: E={result['energy']:.6f} Ha")
                    return result['optimized_xyz']
        except ImportError:
            logger.debug("XTB optimizer not available")
        except Exception as e:
            logger.warning(f"XTB optimization failed: {e}")
        
        return None
    
    def _smiles_to_xyz(self, smiles: str, params) -> str:
        """从SMILES生成XYZ坐标 (使用渐进式生成)"""
        try:
            from app.utils.coordinate_generator import generate_3d_coordinates
            
            result = generate_3d_coordinates(
                smiles=smiles,
                molecule_name=params.molecule_name,
                charge=params.charge,
                multiplicity=params.multiplicity,
                enable_xtb=True
            )
            
            if result.source == 'random':
                logger.warning(
                    f"⚠ PySCF using random coordinates for {smiles}\n"
                    f"   Quality: {result.quality}, Min dist: {result.min_distance:.2f} Å"
                )
            else:
                logger.info(f"✓ PySCF coords generated via {result.source}")
            
            return result.xyz_content
            
        except Exception as e:
            raise ValueError(f"Failed to generate XYZ from SMILES {smiles}: {e}")
    
    def _xyz_to_pyscf_format(self, xyz_content: str) -> str:
        """将XYZ格式转换为PySCF atom字符串格式"""
        lines = xyz_content.strip().split('\n')
        
        # 跳过前两行(原子数和注释)
        coord_lines = []
        for line in lines[2:]:
            parts = line.split()
            if len(parts) >= 4:
                # PySCF格式: "Symbol x y z"
                coord_lines.append(f"{parts[0]} {parts[1]} {parts[2]} {parts[3]}")
        
        # 用分号连接(PySCF格式)
        return '; '.join(coord_lines)
    
    def _get_scf_params(self, params: QCCalculationInput) -> str:
        """生成SCF收敛参数代码"""
        # 估算分子大小
        if params.coordinates:
            num_atoms = len([l for l in params.coordinates.split('\n') if l.strip() and not l.startswith('#')])
        else:
            num_atoms = len(params.smiles) // 2
        
        molecule_info = {
            'num_atoms': num_atoms,
            'charge': params.charge,
            'multiplicity': params.multiplicity
        }
        
        scf_opts = get_scf_convergence_params(molecule_info)
        
        # 生成代码
        code_lines = []
        
        if 'max_cycle' in scf_opts:
            code_lines.append(f"mf.max_cycle = {scf_opts['max_cycle']}")
        
        if 'conv_tol' in scf_opts:
            code_lines.append(f"mf.conv_tol = {scf_opts['conv_tol']}")
        
        if 'level_shift' in scf_opts:
            code_lines.append(f"mf.level_shift = {scf_opts['level_shift']}")
        
        if 'diis' in scf_opts and scf_opts['diis']:
            code_lines.append("mf.diis = True")
        
        return '\n    '.join(code_lines) if code_lines else "# Using default SCF parameters"
    
    def _get_solvent_setup(self, params: QCCalculationInput) -> str:
        """生成溶剂模型设置代码"""
        if params.solvent_model == 'gas':
            return "# Gas phase calculation"
        
        if params.solvent_model in ['pcm', 'smd']:
            # 获取介电常数
            eps = get_solvent_dielectric(params.solvent_name)
            return f"""    # Solvent model: {params.solvent_model}
    from pyscf import solvent
    mf = solvent.PCM(mf)
    mf.with_solvent.eps = {eps}"""
        
        return "# No solvent model"
    
    def parse_output(self, output_file: Path) -> QCCalculationResult:
        """解析PySCF JSON输出(增强版)"""
        # 智能文件查找
        results_file = self._find_results_file(output_file.parent)
        
        if not results_file:
            return QCCalculationResult(
                energy_au=0.0,
                success=False,
                error_message="PySCF results file not found"
            )
        
        try:
            with open(results_file) as f:
                data = json.load(f)
            
            # 验证结果
            is_valid, error_msg = self._validate_results(data)
            if not is_valid:
                logger.warning(f"Invalid PySCF results: {error_msg}")
                data['success'] = False
                data['error_message'] = error_msg
            
            return QCCalculationResult(**data)
        except Exception as e:
            return QCCalculationResult(
                energy_au=0.0,
                success=False,
                error_message=f"Failed to parse PySCF results: {str(e)}"
            )
    
    def _find_results_file(self, work_dir: Path) -> Optional[Path]:
        """智能查找PySCF结果文件"""
        # 1. 标准名称
        standard = work_dir / "pyscf_results.json"
        if standard.exists():
            return standard
        
        # 2. 任何*_results.json文件
        results_files = list(work_dir.glob("*_results.json"))
        if results_files:
            logger.info(f"Found results file: {results_files[0].name}")
            return results_files[0]
        
        # 3. 任何.json文件(验证内容)
        json_files = list(work_dir.glob("*.json"))
        for json_file in json_files:
            try:
                with open(json_file) as f:
                    data = json.load(f)
                if 'energy_au' in data or 'energy' in data:
                    logger.info(f"Found valid results in: {json_file.name}")
                    return json_file
            except:
                continue
        
        logger.error(f"No results file found in {work_dir}")
        return None
    
    def _validate_results(self, results: Dict[str, Any]) -> Tuple[bool, str]:
        """验证PySCF结果的合理性"""
        # 必需字段
        if 'energy_au' not in results:
            return False, "Missing energy_au field"
        
        energy = results['energy_au']
        
        # 合理性检查
        if energy == 0.0:
            return False, "Energy is zero (suspicious)"
        
        if abs(energy) > 100000:
            return False, f"Energy too large: {energy} Ha"
        
        # 收敛检查
        if 'converged' in results and not results['converged']:
            return False, "SCF not converged"
        
        # HOMO-LUMO gap检查
        if 'homo' in results and 'lumo' in results:
            if results['homo'] is not None and results['lumo'] is not None:
                gap = results['lumo'] - results['homo']
                if gap < 0:
                    return False, f"Negative HOMO-LUMO gap: {gap} eV"
                if gap < 0.05:
                    logger.warning(f"Very small HOMO-LUMO gap: {gap} eV")
        
        return True, "OK"

    
    def get_slurm_script(
        self,
        input_file: Path,
        job_name: str,
        partition: str,
        nodes: int,
        cpus: int,
        time_limit: int
    ) -> str:
        """生成PySCF Slurm脚本"""
        return f"""#!/bin/bash
#SBATCH --job-name={job_name[:64]}
#SBATCH --output=pyscf_out.log
#SBATCH --error=pyscf_err.log
#SBATCH --partition={partition}
#SBATCH --nodes={nodes}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --time={time_limit}

# 进入工作目录
cd $SLURM_SUBMIT_DIR

# 激活Python环境(如果有conda)
# source activate molyte

# 设置OpenMP线程数
export OMP_NUM_THREADS={cpus}

# 运行PySCF
python {input_file.name} > pyscf_calculation.log 2>&1

echo "PySCF calculation completed"
"""
