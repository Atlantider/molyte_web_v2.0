#!/usr/bin/env python3
"""
修复 cluster 和 cluster_minus 任务的电荷错误

问题：cluster 和 cluster_minus 任务的电荷被错误地设置为中心离子的电荷（例如 +1），
而应该是中心离子加上配体的总电荷。

例如：
- cluster: Li+ + 3EFA + 2TFSI → charge = 1 + 0 + (-1)*2 = -1（但数据库中是 1，错误）
- cluster_minus: Li+ + EC + 3PF6- → charge = 1 + 0 + (-1)*3 = -2（但数据库中是 1，错误）

此脚本会修复所有 cluster 和 cluster_minus 任务的电荷。
"""

import sys
sys.path.insert(0, '/opt/molyte_web_v1.0/backend')

from app.database import SessionLocal
from app.models.qc import QCJob
from app.models.result import SolvationStructure
from app.api.v1.cluster_analysis import MOLECULE_INFO_MAP
import re

# 定义配体的电荷
LIGAND_CHARGES = {
    'EC': 0,      # 中性
    'DEC': 0,     # 中性
    'DME': 0,     # 中性
    'EMC': 0,     # 中性
    'MP': 0,      # 中性
    'MPN': 0,     # 中性
    'EFA': 0,     # 中性
    'VC': 0,      # 中性
    'FS': 0,      # 中性
    'PF6': -1,    # 阴离子
    'BF4': -1,    # 阴离子
    'TFSI': -1,   # 阴离子
    'FSI': -1,    # 阴离子
    'DFOB': -1,   # 阴离子
    'ClO4': -1,   # 阴离子
    'NO3': -1,    # 阴离子
}

# 中心离子的电荷
CENTER_ION_CHARGES = {
    'Li': 1,
    'Na': 1,
    'K': 1,
    'Mg': 2,
    'Ca': 2,
    'Zn': 2,
}

def calculate_cluster_minus_charge(center_ion: str, original_composition: dict, removed_composition: dict) -> int:
    """计算 cluster_minus 的正确电荷"""
    center_charge = CENTER_ION_CHARGES.get(center_ion, 1)
    total_charge = center_charge
    
    for mol_name, count in original_composition.items():
        if mol_name not in removed_composition:
            # 这个分子没有被移除
            mol_charge = LIGAND_CHARGES.get(mol_name, 0)
            total_charge += mol_charge * count
    
    return total_charge

def main():
    db = SessionLocal()

    try:
        # 查找所有 cluster 和 cluster_minus 任务
        cluster_jobs = db.query(QCJob).filter(
            (QCJob.task_type == 'cluster') | (QCJob.task_type.like('cluster_minus_%'))
        ).all()

        print(f"Found {len(cluster_jobs)} cluster/cluster_minus tasks")
        print("-" * 100)

        fixed_count = 0
        errors = []

        for job in cluster_jobs:
            task_type = job.task_type

            # 获取中心离子和原始组成
            center_ion = None
            original_composition = None
            if job.solvation_structure_id:
                struct = db.query(SolvationStructure).filter(
                    SolvationStructure.id == job.solvation_structure_id
                ).first()
                if struct:
                    center_ion = struct.center_ion
                    original_composition = struct.composition or {}

            if not center_ion or not original_composition:
                errors.append(f"Job {job.id}: Missing center_ion or composition")
                continue

            # 计算正确的电荷
            center_charge = CENTER_ION_CHARGES.get(center_ion, 1)
            correct_charge = center_charge

            if task_type == 'cluster':
                # 对于 cluster 任务，计算所有配体的电荷
                for mol_name, count in original_composition.items():
                    if mol_name not in ["Li", "Na", "K", "Mg", "Ca", "Zn"]:
                        mol_charge = LIGAND_CHARGES.get(mol_name, 0)
                        correct_charge += mol_charge * count
            elif task_type.startswith('cluster_minus_'):
                # 对于 cluster_minus 任务，计算剩余配体的电荷
                parts = task_type.replace('cluster_minus_', '').split('_')
                removed_composition = {}
                i = 0
                while i < len(parts):
                    if i + 1 < len(parts) and parts[i + 1].isdigit():
                        mol_name = parts[i]
                        count = int(parts[i + 1])
                        removed_composition[mol_name] = count
                        i += 2
                    else:
                        i += 1

                for mol_name, count in original_composition.items():
                    if mol_name not in removed_composition:
                        mol_charge = LIGAND_CHARGES.get(mol_name, 0)
                        correct_charge += mol_charge * count

            # 如果电荷不匹配，更新
            if job.charge != correct_charge:
                print(f"Fixing Job {job.id}: {job.molecule_name}")
                print(f"  Task Type: {task_type}")
                print(f"  Center Ion: {center_ion}")
                print(f"  Composition: {original_composition}")
                print(f"  Old Charge: {job.charge} → New Charge: {correct_charge}")

                job.charge = correct_charge
                fixed_count += 1
        
        # 提交更改
        if fixed_count > 0:
            db.commit()
            print(f"\n✅ Fixed {fixed_count} tasks")
        else:
            print(f"\n✅ All tasks already have correct charges")
        
        if errors:
            print(f"\n⚠️  Errors encountered:")
            for error in errors:
                print(f"  {error}")
    
    finally:
        db.close()

if __name__ == "__main__":
    main()

