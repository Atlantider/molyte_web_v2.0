"""
Gaussian 工具模块

提供 Gaussian 输入文件生成和输出解析功能
"""
import logging
import re
from typing import Dict, Any, Optional
from pathlib import Path


logger = logging.getLogger(__name__)


class GaussianUtils:
    """Gaussian 工具类"""
    
    # Hartree 到 kcal/mol 的转换因子
    HARTREE_TO_KCAL = 627.509
    
    def __init__(self):
        self.logger = logging.getLogger(self.__class__.__name__)
    
    def generate_input_file(
        self,
        output_path: Path,
        molecule_name: str,
        xyz_content: str = '',
        smiles: str = '',
        charge: int = 0,
        spin_multiplicity: int = 1,
        functional: str = 'B3LYP',
        basis_set: str = '6-31++g(d,p)',
        solvent_model: str = 'gas',
        solvent_name: str = '',
        nprocs: int = 16,
        mem: str = '32GB',
        job_type: str = 'opt freq'
    ):
        """
        生成 Gaussian 输入文件
        
        Args:
            output_path: 输出文件路径
            molecule_name: 分子名称
            xyz_content: XYZ 坐标内容
            smiles: SMILES 字符串（如果没有 xyz_content）
            charge: 电荷
            spin_multiplicity: 自旋多重度
            functional: 泛函
            basis_set: 基组
            solvent_model: 溶剂模型 (gas, pcm, smd)
            solvent_name: 溶剂名称
            nprocs: CPU 核心数
            mem: 内存
            job_type: 计算类型
        """
        # 构建方法行
        method = f"#{functional}/{basis_set} {job_type}"
        
        # 添加溶剂模型
        if solvent_model.lower() == 'pcm' and solvent_name:
            method += f" SCRF=(PCM,Solvent={solvent_name})"
        elif solvent_model.lower() == 'smd' and solvent_name:
            method += f" SCRF=(SMD,Solvent={solvent_name})"
        
        # 获取坐标
        if xyz_content:
            coords = self._parse_xyz_content(xyz_content)
        elif smiles:
            coords = self._smiles_to_coords(smiles)
        else:
            raise ValueError("Must provide either xyz_content or smiles")
        
        # 构建输入文件
        content = f"""%NProcShared={nprocs}
%Mem={mem}
%Chk={molecule_name}.chk
{method}

{molecule_name}

{charge} {spin_multiplicity}
{coords}

"""
        
        with open(output_path, 'w') as f:
            f.write(content)
        
        self.logger.info(f"生成 Gaussian 输入文件: {output_path}")
    
    def parse_output(self, work_dir: Path) -> Dict[str, Any]:
        """
        解析 Gaussian 输出
        
        Args:
            work_dir: 工作目录
            
        Returns:
            解析结果
        """
        log_files = list(work_dir.glob('*.log'))
        if not log_files:
            return {'success': False, 'error': 'No log file found'}
        
        log_file = log_files[0]
        
        try:
            with open(log_file, 'r') as f:
                content = f.read()
            
            result = {
                'success': False,
                'log_file': str(log_file)
            }
            
            # 检查是否正常终止
            if 'Normal termination of Gaussian' in content:
                result['success'] = True
            else:
                # 检查错误类型
                if 'Convergence failure' in content or 'SCF Done' not in content:
                    result['error'] = 'SCF convergence failure'
                elif 'Error' in content:
                    error_match = re.search(r'Error[^\n]+', content)
                    result['error'] = error_match.group(0) if error_match else 'Unknown error'
                else:
                    result['error'] = 'Calculation did not complete normally'
                return result
            
            # 提取能量
            energy_match = re.search(r'SCF Done:.*?=\s*([-\d.]+)', content)
            if energy_match:
                result['energy_au'] = float(energy_match.group(1))
                result['energy_kcal'] = result['energy_au'] * self.HARTREE_TO_KCAL
            
            # 提取 HOMO/LUMO
            homo_match = re.search(r'Alpha\s+occ\.\s+eigenvalues.*?([-\d.]+)\s*$', content, re.MULTILINE)
            lumo_match = re.search(r'Alpha\s+virt\.\s+eigenvalues.*?([-\d.]+)', content)
            
            if homo_match:
                result['homo_au'] = float(homo_match.group(1))
            if lumo_match:
                result['lumo_au'] = float(lumo_match.group(1))
            
            # 提取优化后的坐标
            coords = self._extract_optimized_coords(content)
            if coords:
                result['optimized_xyz'] = coords
            
            # 提取偶极矩
            dipole_match = re.search(r'Tot=\s*([\d.]+)', content)
            if dipole_match:
                result['dipole_moment'] = float(dipole_match.group(1))
            
            # 提取零点能
            zpe_match = re.search(r'Zero-point correction=\s*([\d.]+)', content)
            if zpe_match:
                result['zero_point_energy'] = float(zpe_match.group(1))
            
            return result
            
        except Exception as e:
            self.logger.error(f"解析输出失败: {e}", exc_info=True)
            return {'success': False, 'error': str(e)}
    
    def _parse_xyz_content(self, xyz_content: str) -> str:
        """解析 XYZ 内容，返回 Gaussian 格式坐标"""
        lines = xyz_content.strip().split('\n')
        
        # 跳过 XYZ 头部（原子数和注释行）
        start_idx = 0
        try:
            int(lines[0].strip())
            start_idx = 2
        except ValueError:
            start_idx = 0
        
        coords = []
        for line in lines[start_idx:]:
            parts = line.split()
            if len(parts) >= 4:
                atom = parts[0]
                x, y, z = parts[1], parts[2], parts[3]
                coords.append(f" {atom}    {x}    {y}    {z}")
        
        return '\n'.join(coords)
    
    def _smiles_to_coords(self, smiles: str) -> str:
        """将 SMILES 转换为 3D 坐标"""
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Cannot parse SMILES: {smiles}")
            
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
            
            conf = mol.GetConformer()
            coords = []
            for i, atom in enumerate(mol.GetAtoms()):
                pos = conf.GetAtomPosition(i)
                coords.append(f" {atom.GetSymbol()}    {pos.x:.6f}    {pos.y:.6f}    {pos.z:.6f}")
            
            return '\n'.join(coords)
            
        except ImportError:
            raise ImportError("RDKit is required for SMILES to coords conversion")
    
    def _extract_optimized_coords(self, content: str) -> Optional[str]:
        """从输出中提取优化后的坐标"""
        # 查找 "Standard orientation" 或 "Input orientation"
        pattern = r'Standard orientation:.*?-{50,}\n(.*?)-{50,}'
        matches = re.findall(pattern, content, re.DOTALL)
        
        if matches:
            # 取最后一个（最终优化的结构）
            table = matches[-1]
            coords = []
            
            for line in table.strip().split('\n'):
                parts = line.split()
                if len(parts) == 6:
                    try:
                        atomic_num = int(parts[1])
                        x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                        
                        # 原子序数到元素符号
                        elements = {
                            1: 'H', 6: 'C', 7: 'N', 8: 'O', 9: 'F',
                            15: 'P', 16: 'S', 17: 'Cl', 3: 'Li'
                        }
                        symbol = elements.get(atomic_num, 'X')
                        coords.append(f"{symbol}  {x:.6f}  {y:.6f}  {z:.6f}")
                    except ValueError:
                        continue
            
            if coords:
                return '\n'.join([str(len(coords)), ''] + coords)
        
        return None
