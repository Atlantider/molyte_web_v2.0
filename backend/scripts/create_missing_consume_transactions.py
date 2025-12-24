#!/usr/bin/env python3
"""
为历史任务创建缺失的消费交易记录

这个脚本用于为那些已完成但没有消费交易的任务创建交易记录。
"""

import sys
import os
sys.path.append('/opt/molyte_web_v1.0/backend')
os.environ['DATABASE_URL'] = 'postgresql://molyte:molyte2025@127.0.0.1:5432/molyte_db'

from app.database import SessionLocal
from app.models.job import MDJob, JobStatus
from app.models.billing import QuotaTransaction
from app.models.user import User
from datetime import datetime

def create_missing_transactions():
    """为历史任务创建缺失的消费交易"""
    db = SessionLocal()
    
    try:
        # 查询所有已完成且有核时数据的任务
        completed_jobs = db.query(MDJob).filter(
            MDJob.status == JobStatus.COMPLETED,
            MDJob.actual_cpu_hours > 0
        ).all()
        
        print(f"Found {len(completed_jobs)} completed jobs with CPU hours\n")
        
        # 检查哪些任务没有消费交易
        jobs_without_transactions = []
        for job in completed_jobs:
            # 查询该任务是否有消费交易
            transaction = db.query(QuotaTransaction).filter(
                QuotaTransaction.user_id == job.user_id,
                QuotaTransaction.type == 'consume',
                QuotaTransaction.description.contains(f'Job {job.id}') if QuotaTransaction.description else False
            ).first()
            
            if not transaction:
                jobs_without_transactions.append(job)
        
        print(f"Found {len(jobs_without_transactions)} jobs without consume transactions\n")
        
        if not jobs_without_transactions:
            print("All jobs have transactions")
            return
        
        total_created = 0
        total_cpu_hours = 0.0
        
        for job in jobs_without_transactions:
            # 计算总核时
            total_cpu_hours_job = job.actual_cpu_hours + (job.resp_cpu_hours or 0)

            # 获取当前用户余额
            user = db.query(User).filter(User.id == job.user_id).first()
            balance_before = user.balance_cpu_hours
            balance_after = balance_before - total_cpu_hours_job

            # 创建消费交易
            transaction = QuotaTransaction(
                user_id=job.user_id,
                type='consume',
                amount=-total_cpu_hours_job,  # 负数表示消费
                balance_before=balance_before,
                balance_after=balance_after,
                reference_id=job.id,
                reference_type='job',
                description=f'MD Job {job.id} completed (MD: {job.actual_cpu_hours:.2f}h + RESP: {job.resp_cpu_hours or 0:.2f}h)'
            )
            db.add(transaction)

            # 更新用户余额
            user.balance_cpu_hours = balance_after
            
            print(f"Job {job.id} (User {job.user_id}):")
            print(f"  CPU hours: {total_cpu_hours_job:.2f}h")
            print(f"  Transaction created: -{total_cpu_hours_job:.2f}h")
            
            total_created += 1
            total_cpu_hours += total_cpu_hours_job
            print()
        
        # 提交更改
        db.commit()
        
        print(f"✅ Created {total_created} transactions")
        print(f"✅ Total CPU hours consumed: {total_cpu_hours:.2f}h")
        
        # 显示受影响的用户
        print("\nAffected users:")
        affected_user_ids = set(job.user_id for job in jobs_without_transactions)
        for user_id in affected_user_ids:
            user = db.query(User).filter(User.id == user_id).first()
            user_jobs = [j for j in jobs_without_transactions if j.user_id == user_id]
            user_cpu_hours = sum(j.actual_cpu_hours + (j.resp_cpu_hours or 0) for j in user_jobs)
            
            # 重新计算用户余额
            from sqlalchemy import func
            total_recharged = db.query(func.sum(QuotaTransaction.amount)).filter(
                QuotaTransaction.user_id == user_id,
                QuotaTransaction.type.in_(['recharge', 'points_exchange', 'admin_adjust', 'bonus'])
            ).scalar() or 0.0
            
            total_consumed = db.query(func.sum(QuotaTransaction.amount)).filter(
                QuotaTransaction.user_id == user_id,
                QuotaTransaction.type == 'consume'
            ).scalar() or 0.0
            total_consumed = abs(total_consumed)
            
            new_balance = user.free_cpu_hours_granted + total_recharged - total_consumed
            
            print(f"  {user.username}:")
            print(f"    Jobs: {len(user_jobs)}")
            print(f"    CPU hours: {user_cpu_hours:.2f}h")
            print(f"    New balance: {new_balance:.2f}h (was {user.balance_cpu_hours:.2f}h)")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        db.rollback()
    finally:
        db.close()

if __name__ == '__main__':
    create_missing_transactions()

