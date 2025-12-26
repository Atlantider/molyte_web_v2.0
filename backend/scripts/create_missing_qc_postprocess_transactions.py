#!/usr/bin/env python3
"""
为QC和后处理任务创建缺失的消费交易记录

这个脚本用于为那些已完成但没有消费交易的QC和后处理任务创建交易记录。
"""

import sys
import os
sys.path.append('/opt/molyte_web_v1.0/backend')
os.environ['DATABASE_URL'] = 'postgresql://molyte:molyte2025@127.0.0.1:5432/molyte_db'

from app.database import SessionLocal
from app.models.qc import QCJob, QCJobStatus
from app.models.job import PostprocessJob, JobStatus, MDJob
from app.models.billing import QuotaTransaction
from app.models.user import User
from datetime import datetime

def create_qc_transactions():
    """为QC任务创建消费交易"""
    db = SessionLocal()
    
    try:
        # 查询所有已完成且有核时数据的QC任务
        completed_jobs = db.query(QCJob).filter(
            QCJob.status == QCJobStatus.COMPLETED,
            QCJob.actual_cpu_hours > 0
        ).all()
        
        print(f"Found {len(completed_jobs)} completed QC jobs with CPU hours\n")
        
        # 检查哪些任务没有消费交易
        jobs_without_transactions = []
        for job in completed_jobs:
            # 查询该任务是否有消费交易
            transaction = db.query(QuotaTransaction).filter(
                QuotaTransaction.user_id == job.user_id,
                QuotaTransaction.type == 'consume',
                QuotaTransaction.reference_type == 'qc_job',
                QuotaTransaction.reference_id == job.id
            ).first()
            
            if not transaction:
                jobs_without_transactions.append(job)
        
        print(f"Found {len(jobs_without_transactions)} QC jobs without consume transactions\n")
        
        if not jobs_without_transactions:
            print("All QC jobs have transactions")
            return 0.0
        
        total_created = 0
        total_cpu_hours = 0.0
        
        for job in jobs_without_transactions:
            # 获取用户
            user = db.query(User).filter(User.id == job.user_id).first()
            if not user:
                print(f"⚠️  User {job.user_id} not found for QC job {job.id}")
                continue
            
            # 获取当前用户余额
            balance_before = user.balance_cpu_hours
            balance_after = balance_before - job.actual_cpu_hours
            
            # 创建消费交易
            transaction = QuotaTransaction(
                user_id=job.user_id,
                type='consume',
                amount=-job.actual_cpu_hours,
                balance_before=balance_before,
                balance_after=balance_after,
                reference_id=job.id,
                reference_type='qc_job',
                description=f'QC Job {job.id} completed ({job.actual_cpu_hours:.2f}h)'
            )
            db.add(transaction)
            
            # 更新用户余额
            user.balance_cpu_hours = balance_after
            
            total_created += 1
            total_cpu_hours += job.actual_cpu_hours
        
        # 提交更改
        db.commit()
        
        print(f"✅ Created {total_created} QC transactions")
        print(f"✅ Total CPU hours consumed: {total_cpu_hours:.2f}h\n")
        
        return total_cpu_hours
        
    except Exception as e:
        print(f"❌ Error creating QC transactions: {e}")
        import traceback
        traceback.print_exc()
        db.rollback()
        return 0.0
    finally:
        db.close()


def create_postprocess_transactions():
    """为后处理任务创建消费交易"""
    db = SessionLocal()
    
    try:
        # 查询所有已完成且有核时数据的后处理任务
        completed_jobs = db.query(PostprocessJob).filter(
            PostprocessJob.status == JobStatus.COMPLETED,
            PostprocessJob.actual_cpu_hours > 0
        ).all()
        
        print(f"Found {len(completed_jobs)} completed postprocess jobs with CPU hours\n")
        
        # 检查哪些任务没有消费交易
        jobs_without_transactions = []
        for job in completed_jobs:
            # 查询该任务是否有消费交易
            transaction = db.query(QuotaTransaction).filter(
                QuotaTransaction.reference_type == 'postprocess_job',
                QuotaTransaction.reference_id == job.id
            ).first()
            
            if not transaction:
                jobs_without_transactions.append(job)
        
        print(f"Found {len(jobs_without_transactions)} postprocess jobs without consume transactions\n")
        
        if not jobs_without_transactions:
            print("All postprocess jobs have transactions")
            return 0.0
        
        total_created = 0
        total_cpu_hours = 0.0
        
        for job in jobs_without_transactions:
            # 获取关联的MD任务以获取用户ID
            md_job = db.query(MDJob).filter(MDJob.id == job.md_job_id).first()
            if not md_job:
                print(f"⚠️  MD Job {job.md_job_id} not found for postprocess job {job.id}")
                continue
            
            # 获取用户
            user = db.query(User).filter(User.id == md_job.user_id).first()
            if not user:
                print(f"⚠️  User {md_job.user_id} not found for postprocess job {job.id}")
                continue
            
            # 获取当前用户余额
            balance_before = user.balance_cpu_hours
            balance_after = balance_before - job.actual_cpu_hours
            
            # 创建消费交易
            transaction = QuotaTransaction(
                user_id=md_job.user_id,
                type='consume',
                amount=-job.actual_cpu_hours,
                balance_before=balance_before,
                balance_after=balance_after,
                reference_id=job.id,
                reference_type='postprocess_job',
                description=f'Postprocess Job {job.id} completed ({job.actual_cpu_hours:.2f}h)'
            )
            db.add(transaction)
            
            # 更新用户余额
            user.balance_cpu_hours = balance_after
            
            total_created += 1
            total_cpu_hours += job.actual_cpu_hours
        
        # 提交更改
        db.commit()
        
        print(f"✅ Created {total_created} postprocess transactions")
        print(f"✅ Total CPU hours consumed: {total_cpu_hours:.2f}h\n")
        
        return total_cpu_hours
        
    except Exception as e:
        print(f"❌ Error creating postprocess transactions: {e}")
        import traceback
        traceback.print_exc()
        db.rollback()
        return 0.0
    finally:
        db.close()


if __name__ == '__main__':
    print("=" * 60)
    print("Creating Missing QC and Postprocess Transactions")
    print("=" * 60)
    print()
    
    qc_hours = create_qc_transactions()
    postprocess_hours = create_postprocess_transactions()
    
    print("=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"QC transactions CPU hours: {qc_hours:.2f}h")
    print(f"Postprocess transactions CPU hours: {postprocess_hours:.2f}h")
    print(f"Total CPU hours: {qc_hours + postprocess_hours:.2f}h")

