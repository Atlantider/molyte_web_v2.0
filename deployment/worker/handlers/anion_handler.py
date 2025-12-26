"""
阴离子生成任务处理器

处理自定义阴离子力场参数生成
"""
import logging
from typing import Dict, Any

from worker.handlers.base import BaseHandler


logger = logging.getLogger(__name__)


class AnionHandler(BaseHandler):
    """阴离子生成任务处理器"""
    
    JOB_TYPE = "ANION_GENERATION"
    
    def process(self, job: Dict[str, Any]) -> bool:
        """处理阴离子生成任务"""
        job_id = job['id']
        config = job.get('config', {})
        status = config.get('status', 'PENDING')
        
        self.logger.info(f"开始处理阴离子生成任务 {job_id} (状态: {status})")
        
        try:
            if status == 'PENDING':
                # Phase 1: 解析输入，创建 QC 任务
                return self._process_pending(job_id, config)
            elif status == 'QC_PENDING':
                # Phase 2: 处理 QC 结果，生成力场
                return self._process_qc_pending(job_id, config)
            else:
                self.logger.warning(f"未知状态: {status}")
                return False
                
        except Exception as e:
            self.handle_error(job_id, e)
            return False
    
    def handle_completion(self, job_id: int, job_info: Dict, slurm_status: str):
        """处理 QC 完成"""
        pass
    
    def _process_pending(self, job_id: int, config: Dict) -> bool:
        """Phase 1: 解析输入，创建 QC 任务"""
        try:
            anion_name = config.get('anion_name')
            identifier_type = config.get('identifier_type')
            identifier_value = config.get('identifier_value')
            charge = config.get('charge', -1)
            
            # 解析分子结构
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            if identifier_type == 'smiles':
                mol = Chem.MolFromSmiles(identifier_value)
            elif identifier_type == 'name':
                # 通过名称查找
                mol = self._get_mol_from_name(identifier_value)
            else:
                raise ValueError(f"Unknown identifier type: {identifier_type}")
            
            if mol is None:
                raise ValueError(f"Cannot parse molecule: {identifier_value}")
            
            # 生成 3D 坐标
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
            
            # 生成 XYZ
            xyz_content = self._mol_to_xyz(mol)
            
            # 创建 QC 任务
            qc_result = self.client.post(
                "/api/v1/qc/jobs",
                data={
                    'molecule_name': anion_name,
                    'xyz_content': xyz_content,
                    'charge': charge,
                    'spin_multiplicity': 1,
                    'basis_set': '6-31++g(d,p)',
                    'functional': 'B3LYP',
                    'job_type': 'opt_freq'
                }
            )
            
            if qc_result:
                # 更新状态为 QC_PENDING
                self.client.update_job_status(
                    job_id=job_id,
                    job_type='ANION_GENERATION',
                    status='QC_PENDING',
                    qc_job_id=qc_result['id']
                )
                return True
            else:
                raise ValueError("Failed to create QC job")
                
        except Exception as e:
            self.handle_error(job_id, e)
            return False
    
    def _process_qc_pending(self, job_id: int, config: Dict) -> bool:
        """Phase 2: 处理 QC 结果，生成力场"""
        try:
            qc_job_id = config.get('qc_job_id')
            anion_name = config.get('anion_name')
            
            # 获取 QC 结果
            qc_result = self.client.get(f"/api/v1/qc/{qc_job_id}")
            
            if not qc_result or qc_result.get('status') != 'COMPLETED':
                self.logger.info(f"QC 任务 {qc_job_id} 尚未完成")
                return False
            
            # 生成力场文件
            from worker.utils.file_utils import FileUtils
            
            file_utils = FileUtils(self.config)
            ff_result = file_utils.generate_anion_forcefield(
                anion_name=anion_name,
                qc_result=qc_result
            )
            
            if ff_result.get('success'):
                # 注册阴离子
                register_result = self.client.post(
                    f"/api/v1/workers/anion/{job_id}/register",
                    data=ff_result
                )
                
                if register_result:
                    self.update_status(job_id, 'COMPLETED')
                    return True
            
            self.update_status(
                job_id, 'FAILED',
                error_message=ff_result.get('error', 'Force field generation failed')
            )
            return False
            
        except Exception as e:
            self.handle_error(job_id, e)
            return False
    
    def _get_mol_from_name(self, name: str):
        """通过名称获取分子"""
        # 简单的名称到SMILES映射
        name_to_smiles = {
            'PF6': 'F[P-](F)(F)(F)(F)F',
            'TFSI': 'FC(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F',
            'BF4': 'F[B-](F)(F)F',
        }
        
        smiles = name_to_smiles.get(name.upper())
        if smiles:
            from rdkit import Chem
            return Chem.MolFromSmiles(smiles)
        return None
    
    def _mol_to_xyz(self, mol) -> str:
        """将RDKit分子转为XYZ格式"""
        from rdkit import Chem
        
        conf = mol.GetConformer()
        lines = [str(mol.GetNumAtoms()), ""]
        
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            lines.append(f"{atom.GetSymbol()}  {pos.x:.6f}  {pos.y:.6f}  {pos.z:.6f}")
        
        return "\n".join(lines)
