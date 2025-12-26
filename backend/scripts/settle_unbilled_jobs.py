#!/usr/bin/env python3
"""
补结算脚本 - 处理未计费的已完成和已取消任务

业务规则:
- COMPLETED任务: 必须计费
- CANCELLED任务: 必须计费(用户主动取消)
- FAILED任务: 不计费(可能是系统问题)

使用方法:
    python scripts/settle_unbilled_jobs.py [--dry-run]
"""

import sys
import os
from pathlib import Path

# 添加项目根目录到路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from app.database import SessionLocal
from app.models.job import MDJob, JobStatus
from app.services.billing import BillingService
import argparse


def main():
    parser = argparse.ArgumentParser(description='补结算未计费任务')
    parser.add_argument('--dry-run', action='store_true', 
                       help='只显示将要结算的任务,不实际执行')
    args = parser.parse_args()
    
    db = SessionLocal()
    
    try:
        # 查找未结算的COMPLETED和CANCELLED任务
        unbilled_jobs = db.query(MDJob).filter(
            MDJob.billed == False,
            MDJob.status.in_([JobStatus.COMPLETED, JobStatus.CANCELLED])
        ).order_by(MDJob.user_id, MDJob.id).all()
        
        print(f"\n{'='*80}")
        print(f"未计费任务补结算脚本")
        print(f"{'='*80}\n")
        print(f"找到 {len(unbilled_jobs)} 个未结算任务")
        print(f"模式: {'DRY RUN (不实际执行)' if args.dry_run else '实际执行'}\n")
        
        if not unbilled_jobs:
            print("✅ 没有需要补结算的任务")
            return
        
        # 按用户分组统计
        user_stats = {}
        for job in unbilled_jobs:
            if job.user_id not in user_stats:
                user_stats[job.user_id] = {
                    'jobs': [],
                    'total_hours': 0.0
                }
            
            total_hours = (job.actual_cpu_hours or 0.0) + (job.resp_cpu_hours or 0.0)
            user_stats[job.user_id]['jobs'].append(job)
            user_stats[job.user_id]['total_hours'] += total_hours
        
        # 显示统计信息
        print(f"{'用户ID':<10} {'任务数':<10} {'总核时':<15} {'状态分布'}")
        print(f"{'-'*80}")
        
        for user_id, stats in sorted(user_stats.items()):
            status_counts = {}
            for job in stats['jobs']:
                status = job.status.value
                status_counts[status] = status_counts.get(status, 0) + 1
            
            status_str = ', '.join([f"{s}: {c}" for s, c in status_counts.items()])
            print(f"{user_id:<10} {len(stats['jobs']):<10} {stats['total_hours']:<15.2f} {status_str}")
        
        print(f"\n{'='*80}\n")
        
        # 执行结算
        settled_count = 0
        failed_count = 0
        total_hours_settled = 0.0
        
        for job in unbilled_jobs:
            total_hours = (job.actual_cpu_hours or 0.0) + (job.resp_cpu_hours or 0.0)
            
            print(f"\n任务 #{job.id} (用户 {job.user_id}, 状态: {job.status.value})")
            print(f"  MD核时: {job.actual_cpu_hours:.2f}h, RESP核时: {job.resp_cpu_hours:.2f}h")
            print(f"  总计: {total_hours:.2f}h")
            
            if args.dry_run:
                print(f"  [DRY RUN] 将会结算")
                settled_count += 1
                total_hours_settled += total_hours
            else:
                try:
                    success, message = BillingService.settle_job(db, job)
                    if success:
                        settled_count += 1
                        total_hours_settled += total_hours
                        print(f"  ✅ 结算成功: {message}")
                    else:
                        failed_count += 1
                        print(f"  ❌ 结算失败: {message}")
                except Exception as e:
                    failed_count += 1
                    print(f"  ❌ 异常: {str(e)}")
        
        # 最终统计
        print(f"\n{'='*80}")
        print(f"补结算完成")
        print(f"{'='*80}\n")
        print(f"总任务数: {len(unbilled_jobs)}")
        print(f"成功结算: {settled_count}")
        print(f"失败: {failed_count}")
        print(f"总核时: {total_hours_settled:.2f}h")
        
        if args.dry_run:
            print(f"\n⚠️  这是DRY RUN模式,没有实际执行结算")
            print(f"   要实际执行,请运行: python {sys.argv[0]}")
        else:
            print(f"\n✅ 补结算已完成")
        
    finally:
        db.close()


if __name__ == '__main__':
    main()
