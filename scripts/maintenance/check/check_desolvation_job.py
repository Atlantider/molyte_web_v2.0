#!/usr/bin/env python3
"""
检查去溶剂化任务的状态和配置

使用方法：
python check_desolvation_job.py <job_id>
"""

import sys
import os

# 添加 backend 到路径
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'backend'))

from app.database import SessionLocal
from app.models.job import PostprocessJob
from app.models.qc import QCJob
import json

def check_job(job_id: int):
    db = SessionLocal()
    try:
        # 查询任务
        job = db.query(PostprocessJob).filter(PostprocessJob.id == job_id).first()
        
        if not job:
            print(f"❌ 任务 {job_id} 不存在")
            return
        
        print(f"✅ 任务 {job_id} 信息：")
        print(f"  - 类型: {job.job_type}")
        print(f"  - 状态: {job.status}")
        print(f"  - MD Job ID: {job.md_job_id}")
        print(f"  - 创建时间: {job.created_at}")
        print(f"  - 更新时间: {job.updated_at}")
        print(f"  - 开始时间: {job.started_at}")
        print(f"  - 完成时间: {job.finished_at}")
        print(f"  - 进度: {job.progress}%")
        
        if job.config:
            print(f"\n📋 配置信息：")
            print(json.dumps(job.config, indent=2, ensure_ascii=False))
            
            # 检查 QC 任务
            qc_job_ids = job.config.get('qc_job_ids', [])
            reused_qc_job_ids = job.config.get('reused_qc_job_ids', [])
            
            if qc_job_ids or reused_qc_job_ids:
                print(f"\n🔬 QC 任务信息：")
                print(f"  - 新创建的 QC 任务: {len(qc_job_ids)} 个")
                print(f"  - 复用的 QC 任务: {len(reused_qc_job_ids)} 个")
                
                all_qc_ids = qc_job_ids + reused_qc_job_ids
                for qc_id in all_qc_ids:
                    qc_job = db.query(QCJob).filter(QCJob.id == qc_id).first()
                    if qc_job:
                        reused_mark = " (复用)" if qc_id in reused_qc_job_ids else ""
                        print(f"    - QC Job {qc_id}{reused_mark}: {qc_job.molecule_name} - {qc_job.status}")
                    else:
                        print(f"    - QC Job {qc_id}: ❌ 不存在")
            else:
                print(f"\n⚠️  没有找到 QC 任务 ID")
        else:
            print(f"\n⚠️  没有配置信息")
        
        if job.error_message:
            print(f"\n❌ 错误信息: {job.error_message}")
    
    finally:
        db.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("使用方法: python check_desolvation_job.py <job_id>")
        sys.exit(1)
    
    job_id = int(sys.argv[1])
    check_job(job_id)

