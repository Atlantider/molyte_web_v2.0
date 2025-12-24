#!/usr/bin/env python3
"""
回填已完成 MD 任务的系统结构数据

这个脚本用于处理在修复前已完成的 MD 任务，
为它们提取系统结构并保存到 SystemStructure 表，
使得前端能够显示整体溶液结构 (System)。

使用方法:
    python backfill_system_structures.py [--job-id JOB_ID] [--limit LIMIT]
    
    --job-id: 指定单个任务 ID（可选）
    --limit: 限制处理的任务数量（默认：所有）
"""

import sys
import os
from pathlib import Path
from typing import Optional, Dict, Any
import argparse
import logging

# 添加后端路径
backend_path = Path(__file__).parent.parent
sys.path.insert(0, str(backend_path))

from app.database import SessionLocal
from app.models.job import MDJob, JobStatus
from app.models.result import SystemStructure
from app.services.solvation import get_system_structure

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def backfill_system_structure(job_id: int, db) -> bool:
    """
    为单个任务回填系统结构
    
    Args:
        job_id: MD 任务 ID
        db: 数据库会话
        
    Returns:
        True 如果成功，False 如果失败
    """
    try:
        # 获取任务
        job = db.query(MDJob).filter(MDJob.id == job_id).first()
        if not job:
            logger.warning(f"任务 {job_id} 不存在")
            return False
        
        # 检查任务状态
        if job.status != JobStatus.COMPLETED:
            logger.warning(f"任务 {job_id} 状态为 {job.status}，跳过")
            return False
        
        # 检查工作目录
        if not job.work_dir:
            logger.warning(f"任务 {job_id} 没有工作目录，跳过")
            return False
        
        work_dir = Path(job.work_dir)
        if not work_dir.exists():
            logger.warning(f"任务 {job_id} 的工作目录不存在: {work_dir}")
            return False
        
        # 检查是否已有系统结构
        existing = db.query(SystemStructure).filter(
            SystemStructure.md_job_id == job_id
        ).first()
        
        if existing:
            logger.info(f"任务 {job_id} 已有系统结构，跳过")
            return True
        
        # 提取系统结构
        logger.info(f"正在处理任务 {job_id}...")
        system_result = get_system_structure(str(work_dir), frame_idx=-1)
        
        if not system_result or 'xyz_content' not in system_result:
            logger.warning(f"任务 {job_id} 无法提取系统结构")
            return False
        
        # 创建 SystemStructure 记录
        system_structure = SystemStructure(
            md_job_id=job_id,
            frame_index=system_result.get('frame_index', 0),
            total_frames=system_result.get('total_frames', 1),
            atom_count=system_result.get('atom_count', 0),
            box=system_result.get('box', [0, 0, 0]),
            xyz_content=system_result['xyz_content'],
        )
        
        db.add(system_structure)
        db.commit()
        
        logger.info(
            f"✅ 任务 {job_id} 系统结构已保存 "
            f"({system_result.get('atom_count', 0)} 原子, "
            f"帧 {system_result.get('frame_index', 0)}/{system_result.get('total_frames', 1)})"
        )
        return True
        
    except Exception as e:
        logger.error(f"❌ 任务 {job_id} 处理失败: {e}", exc_info=True)
        db.rollback()
        return False


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='回填已完成 MD 任务的系统结构数据'
    )
    parser.add_argument(
        '--job-id',
        type=int,
        help='指定单个任务 ID'
    )
    parser.add_argument(
        '--limit',
        type=int,
        help='限制处理的任务数量'
    )
    
    args = parser.parse_args()
    
    db = SessionLocal()
    
    try:
        if args.job_id:
            # 处理单个任务
            logger.info(f"处理单个任务: {args.job_id}")
            success = backfill_system_structure(args.job_id, db)
            sys.exit(0 if success else 1)
        else:
            # 处理所有已完成的任务
            query = db.query(MDJob).filter(
                MDJob.status == JobStatus.COMPLETED
            ).order_by(MDJob.id.desc())
            
            if args.limit:
                query = query.limit(args.limit)
            
            jobs = query.all()
            logger.info(f"找到 {len(jobs)} 个已完成的任务")
            
            success_count = 0
            failed_count = 0
            skipped_count = 0
            
            for job in jobs:
                if backfill_system_structure(job.id, db):
                    success_count += 1
                else:
                    failed_count += 1
            
            logger.info(
                f"\n处理完成: "
                f"成功 {success_count}, "
                f"失败 {failed_count}"
            )
            
    finally:
        db.close()


if __name__ == '__main__':
    main()

