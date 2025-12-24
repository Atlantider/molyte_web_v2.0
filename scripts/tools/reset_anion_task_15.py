#!/usr/bin/env python3
"""
重置阴离子生成任务 15 的状态，用于测试修复
"""
import sys
sys.path.insert(0, '/public/home/xiaoji/molyte_web/backend')

from app.database import SessionLocal
from app.models.forcefield import AnionGenerationJob, AnionGenerationStatus
from datetime import datetime

def reset_task_15():
    db = SessionLocal()
    try:
        # 查找任务 15
        job = db.query(AnionGenerationJob).filter(AnionGenerationJob.id == 15).first()
        
        if not job:
            print("❌ 任务 15 不存在")
            return False
            
        print(f"📋 当前任务状态:")
        print(f"   ID: {job.id}")
        print(f"   Job ID: {job.job_id}")
        print(f"   阴离子: {job.anion_name}")
        print(f"   状态: {job.status}")
        print(f"   消息: {job.message}")
        print(f"   标识符类型: {job.identifier_type}")
        print(f"   标识符值: {job.identifier_value}")
        
        # 重置状态
        job.status = AnionGenerationStatus.PENDING
        job.message = "Task reset for testing PubChem API fix"
        job.started_at = None
        job.finished_at = None
        job.qc_job_id = None
        job.work_dir = None
        job.lt_path = None
        job.pdb_path = None
        
        db.commit()
        
        print(f"✅ 任务 15 状态已重置为 PENDING")
        print(f"   Worker 将在下次轮询时重新处理此任务")
        
        return True
        
    except Exception as e:
        print(f"❌ 重置失败: {e}")
        db.rollback()
        return False
    finally:
        db.close()

if __name__ == "__main__":
    reset_task_15()
