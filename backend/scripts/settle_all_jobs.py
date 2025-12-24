#!/usr/bin/env python3
"""
批量结算所有未结算的已完成任务

用于修复核时显示问题：确保所有已完成的任务都被正确结算
"""

import sys
import os
from pathlib import Path

# 添加项目根目录到 Python 路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from app.database import SessionLocal
from app.models.job import MDJob, JobStatus
from app.models.user import User
from app.services.billing import BillingService
from app.models.billing import QuotaTransaction

def main():
    """批量结算所有未结算的已完成任务"""
    db = SessionLocal()
    
    try:
        print("=== 批量结算未结算的已完成任务 ===")
        
        # 获取所有未结算的已完成任务
        jobs = db.query(MDJob).filter(
            MDJob.status == JobStatus.COMPLETED,
            MDJob.billed == False
        ).all()
        
        print(f"找到 {len(jobs)} 个未结算的已完成任务")
        
        if not jobs:
            print("没有需要结算的任务")
            return
        
        success_count = 0
        failed_count = 0
        admin_count = 0
        
        for job in jobs:
            user = db.query(User).filter(User.id == job.user_id).first()
            if not user:
                print(f"❌ 任务 {job.id}: 用户不存在")
                failed_count += 1
                continue
                
            print(f"\n结算任务 {job.id}:")
            print(f"  用户: {user.username} (ID: {job.user_id})")
            print(f"  实际核时: {job.actual_cpu_hours}")
            print(f"  RESP核时: {job.resp_cpu_hours}")
            print(f"  用户余额: {user.balance_cpu_hours}")
            
            success, message = BillingService.settle_job(db, job)
            print(f"  结算结果: {message}")
            
            if success:
                if "管理员任务免费" in message:
                    admin_count += 1
                    print("  ✅ 管理员任务免费")
                else:
                    success_count += 1
                    print("  ✅ 结算成功")
                    
                    # 刷新用户数据查看余额变化
                    db.refresh(user)
                    print(f"  结算后余额: {user.balance_cpu_hours}")
            else:
                failed_count += 1
                print("  ❌ 结算失败")
        
        print(f"\n=== 结算统计 ===")
        print(f"总任务数: {len(jobs)}")
        print(f"成功结算: {success_count}")
        print(f"管理员免费: {admin_count}")
        print(f"失败: {failed_count}")
        
        # 检查结算后的消费记录
        print(f"\n=== 检查消费记录 ===")
        consume_transactions = db.query(QuotaTransaction).filter(
            QuotaTransaction.type == 'consume'
        ).all()
        
        total_consumption = sum(abs(t.amount) for t in consume_transactions)
        print(f"总消费交易数: {len(consume_transactions)}")
        print(f"总消费核时: {total_consumption:.2f}")
        
        # 按用户统计消费
        user_consumption = {}
        for t in consume_transactions:
            if t.user_id not in user_consumption:
                user_consumption[t.user_id] = 0
            user_consumption[t.user_id] += abs(t.amount)
        
        print(f"\n=== 用户消费统计 ===")
        for user_id, consumption in user_consumption.items():
            user = db.query(User).filter(User.id == user_id).first()
            username = user.username if user else f"用户{user_id}"
            print(f"{username}: {consumption:.2f} 核时")
            
    except Exception as e:
        print(f"❌ 批量结算失败: {e}")
        import traceback
        traceback.print_exc()
        
    finally:
        db.close()

if __name__ == "__main__":
    main()
