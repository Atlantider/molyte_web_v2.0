#!/usr/bin/env python3
"""
批量结算所有未结算的已完成 Cluster Analysis 任务

修复核时显示为0的问题：
- 查找所有已完成但 cpu_hours_used = 0 的任务
- 重新计算并更新 cpu_hours_used 字段
- 为非管理员用户生成消费记录
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.database import SessionLocal
from app.models.job import AdvancedClusterJob
from app.models.qc import QCJob
from app.models.user import User
from app.services.billing import BillingService
from sqlalchemy import func

def main():
    db = SessionLocal()
    
    try:
        print("=== 批量结算 Cluster Analysis 任务 ===\n")
        
        # 查找所有已完成但未结算的任务
        completed_jobs = db.query(AdvancedClusterJob).filter(
            AdvancedClusterJob.status == 'COMPLETED',
            AdvancedClusterJob.cpu_hours_used == 0.0
        ).all()
        
        print(f"找到 {len(completed_jobs)} 个已完成但未结算的任务")
        
        if not completed_jobs:
            print("没有需要结算的任务")
            return
        
        settled_count = 0
        total_hours = 0.0
        
        for job in completed_jobs:
            print(f"\n处理任务 {job.id}:")
            
            # 获取用户信息
            user = db.query(User).filter(User.id == job.user_id).first()
            if not user:
                print(f"  ❌ 用户 {job.user_id} 不存在，跳过")
                continue
            
            print(f"  用户: {user.username} ({user.role.value})")
            
            # 计算关联 QC 任务的核时
            qc_jobs = db.query(QCJob).filter(
                QCJob.cluster_analysis_job_id == job.id
            ).all()
            
            qc_hours = sum(qc.actual_cpu_hours for qc in qc_jobs)
            print(f"  关联QC任务: {len(qc_jobs)} 个")
            print(f"  QC任务总核时: {qc_hours:.2f}h")
            
            if qc_hours == 0:
                print(f"  ⚠️  QC任务核时为0，跳过结算")
                continue
            
            # 调用结算服务
            success, message = BillingService.settle_cluster_analysis_job(db, job)
            
            if success:
                print(f"  ✅ 结算成功: {message}")
                print(f"  更新后核时: {job.cpu_hours_used:.2f}h")
                print(f"  任务计数: {job.task_count}")
                settled_count += 1
                total_hours += job.cpu_hours_used
            else:
                print(f"  ❌ 结算失败: {message}")
        
        print(f"\n=== 结算完成 ===")
        print(f"成功结算任务数: {settled_count}")
        print(f"总核时: {total_hours:.2f}h")
        
        # 验证结果
        print(f"\n=== 验证结果 ===")
        remaining_jobs = db.query(AdvancedClusterJob).filter(
            AdvancedClusterJob.status == 'COMPLETED',
            AdvancedClusterJob.cpu_hours_used == 0.0
        ).count()
        print(f"剩余未结算任务: {remaining_jobs}")
        
        # 统计所有已结算任务
        all_settled = db.query(AdvancedClusterJob).filter(
            AdvancedClusterJob.status == 'COMPLETED',
            AdvancedClusterJob.cpu_hours_used > 0.0
        ).all()
        
        if all_settled:
            print(f"已结算任务: {len(all_settled)} 个")
            total_settled_hours = sum(job.cpu_hours_used for job in all_settled)
            print(f"已结算总核时: {total_settled_hours:.2f}h")
            
            # 按用户统计
            user_stats = {}
            for job in all_settled:
                user = db.query(User).filter(User.id == job.user_id).first()
                username = user.username if user else f"User#{job.user_id}"
                if username not in user_stats:
                    user_stats[username] = {'count': 0, 'hours': 0.0}
                user_stats[username]['count'] += 1
                user_stats[username]['hours'] += job.cpu_hours_used
            
            print(f"\n按用户统计:")
            for username, stats in sorted(user_stats.items()):
                print(f"  {username}: {stats['count']} 个任务, {stats['hours']:.2f}h")
        
    except Exception as e:
        print(f"❌ 错误: {e}")
        import traceback
        traceback.print_exc()
    finally:
        db.close()

if __name__ == "__main__":
    main()
