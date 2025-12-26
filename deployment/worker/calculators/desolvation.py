"""
去溶剂化能计算模块

支持部分结果计算
"""
import logging
from typing import Dict, Any, List, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from worker.core.client import APIClient


logger = logging.getLogger(__name__)

# 物理常数
HARTREE_TO_KCAL = 627.509


class DesolvationCalculator:
    """去溶剂化能计算器"""
    
    def __init__(self, client: 'APIClient'):
        """初始化计算器"""
        self.client = client
        self.logger = logging.getLogger(self.__class__.__name__)
    
    def calculate_stepwise_partial(
        self,
        config: Dict[str, Any],
        successful_qc_ids: List[int],
        failed_details: List[Dict]
    ) -> Dict[str, Any]:
        """
        计算Stepwise模式去溶剂化能（支持部分结果）
        
        Args:
            config: 任务配置
            successful_qc_ids: 成功的 QC 任务 ID 列表
            failed_details: 失败详情列表
            
        Returns:
            计算结果
        """
        cluster_qc_id = config.get('cluster_qc_job_id')
        ligand_qc_jobs = config.get('ligand_qc_jobs', {})
        cluster_minus_ids = config.get('cluster_minus_job_ids', [])
        cluster_minus_type_mapping = config.get('cluster_minus_type_mapping', {})
        cluster_data = config.get('cluster_data', {})
        
        # 获取 E_cluster
        e_cluster = self._get_qc_energy(cluster_qc_id)
        
        if e_cluster is None:
            raise ValueError("Cannot get cluster energy")
        
        per_ligand_results = []
        skipped_ligands = []
        
        # 遍历每个配体
        ligands = cluster_data.get('ligands', [])
        for i, ligand_info in enumerate(ligands):
            ligand_type = ligand_info.get('ligand_type')
            ligand_label = ligand_info.get('ligand_label')
            ligand_charge = ligand_info.get('charge', 0)
            
            try:
                # 获取 E_ligand
                ligand_key = f"{ligand_type}_{ligand_charge}"
                ligand_qc_id = ligand_qc_jobs.get(ligand_key)
                
                if not ligand_qc_id or ligand_qc_id not in successful_qc_ids:
                    raise ValueError(f"Ligand QC {ligand_qc_id} not successful")
                
                e_ligand = self._get_qc_energy(ligand_qc_id)
                if e_ligand is None:
                    raise ValueError(f"Cannot get ligand energy for {ligand_type}")
                
                # 获取 E_cluster_minus
                cluster_minus_qc_id = cluster_minus_type_mapping.get(ligand_type)
                if not cluster_minus_qc_id and i < len(cluster_minus_ids):
                    cluster_minus_qc_id = cluster_minus_ids[i]
                
                if not cluster_minus_qc_id or cluster_minus_qc_id not in successful_qc_ids:
                    raise ValueError(f"Cluster minus QC not successful")
                
                e_cluster_minus = self._get_qc_energy(cluster_minus_qc_id)
                if e_cluster_minus is None:
                    raise ValueError("Cannot get cluster minus energy")
                
                # 计算去溶剂化能
                delta_e_au = e_cluster - (e_cluster_minus + e_ligand)
                delta_e_kcal = delta_e_au * HARTREE_TO_KCAL
                
                per_ligand_results.append({
                    'ligand_id': i,
                    'ligand_type': ligand_type,
                    'ligand_label': ligand_label,
                    'e_ligand': e_ligand,
                    'e_cluster_minus': e_cluster_minus,
                    'delta_e': delta_e_kcal,
                    'status': 'calculated'
                })
                
            except Exception as e:
                self.logger.warning(f"Skipping ligand {ligand_label}: {e}")
                skipped_ligands.append({
                    'ligand_id': i,
                    'ligand_type': ligand_type,
                    'ligand_label': ligand_label,
                    'status': 'skipped',
                    'reason': str(e)
                })
        
        return {
            'e_cluster': e_cluster,
            'per_ligand_results': per_ligand_results,
            'skipped_ligands': skipped_ligands,
            'total_ligands': len(ligands),
            'calculated_ligands': len(per_ligand_results),
            'failed_details': failed_details
        }
    
    def calculate_full_partial(
        self,
        config: Dict[str, Any],
        successful_qc_ids: List[int],
        failed_details: List[Dict]
    ) -> Dict[str, Any]:
        """
        计算Full模式去溶剂化能（支持部分结果）
        
        改进：即使部分配体失败，也计算可用部分的结果
        """
        cluster_qc_id = config.get('cluster_qc_job_id')
        center_ion_job_id = config.get('center_ion_job_id')
        ligand_qc_jobs = config.get('ligand_qc_jobs', {})
        cluster_data = config.get('cluster_data', {})
        
        # 获取 E_cluster
        e_cluster = self._get_qc_energy(cluster_qc_id)
        if e_cluster is None:
            raise ValueError("Cannot get cluster energy")
        
        # 获取 E_ion
        e_ion = self._get_qc_energy(center_ion_job_id)
        if e_ion is None:
            raise ValueError("Cannot get center ion energy")
        
        # 计算配体能量
        total_ligand_energy = 0.0
        ligand_energies = {}
        failed_ligands = []
        ligands = cluster_data.get('ligands', [])
        
        for ligand_info in ligands:
            ligand_type = ligand_info.get('ligand_type')
            ligand_label = ligand_info.get('ligand_label')
            ligand_charge = ligand_info.get('charge', 0)
            ligand_key = f"{ligand_type}_{ligand_charge}"
            
            if ligand_key not in ligand_energies:
                ligand_qc_id = ligand_qc_jobs.get(ligand_key)
                
                if ligand_qc_id and ligand_qc_id in successful_qc_ids:
                    e_ligand = self._get_qc_energy(ligand_qc_id)
                    if e_ligand is not None:
                        ligand_energies[ligand_key] = e_ligand
                    else:
                        failed_ligands.append({
                            'ligand_type': ligand_type,
                            'reason': 'Cannot get energy'
                        })
                else:
                    failed_ligands.append({
                        'ligand_type': ligand_type,
                        'reason': 'QC job not successful'
                    })
            
            # 累加能量
            if ligand_key in ligand_energies:
                total_ligand_energy += ligand_energies[ligand_key]
        
        # 即使有失败，也计算可用部分
        is_partial = len(failed_ligands) > 0
        total_desolvation_energy = None
        
        if not is_partial or len(ligand_energies) > 0:
            # 计算去溶剂化能（注意：部分结果是不完整的）
            delta_e_au = e_cluster - (e_ion + total_ligand_energy)
            total_desolvation_energy = delta_e_au * HARTREE_TO_KCAL
        
        return {
            'e_cluster': e_cluster,
            'e_ion': e_ion,
            'total_ligand_energy': total_ligand_energy,
            'ligand_energies': ligand_energies,
            'total_desolvation_energy': total_desolvation_energy,
            'is_partial': is_partial,
            'calculated_ligand_count': len(ligand_energies),
            'total_ligand_count': len(set(l.get('ligand_type') for l in ligands)),
            'failed_ligands': failed_ligands,
            'failed_details': failed_details,
            'warning': (
                'This result is partial - some ligand energies are missing'
                if is_partial else None
            )
        }
    
    def _get_qc_energy(self, qc_job_id: int) -> Optional[float]:
        """获取 QC 任务的能量"""
        if not qc_job_id:
            return None
        
        try:
            result = self.client.get(f"/api/v1/qc/{qc_job_id}/result")
            if result and result.get('energy_au') is not None:
                return result['energy_au']
            return None
        except Exception as e:
            self.logger.error(f"获取 QC {qc_job_id} 能量失败: {e}")
            return None
