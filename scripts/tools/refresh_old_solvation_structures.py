#!/usr/bin/env python3
"""
重新生成旧的溶剂化结构，添加 mol_order 字段

这个脚本会：
1. 查找所有 mol_order 为 None 的溶剂化结构
2. 找到对应的 MD 任务
3. 调用 refresh_solvation_structures API 重新分析
"""

import sys
import os
from pathlib import Path

# 添加项目路径
backend_dir = Path(__file__).parent.parent / "backend"
sys.path.insert(0, str(backend_dir))

# 加载环境变量
from dotenv import load_dotenv
env_file = backend_dir / ".env"
if env_file.exists():
    load_dotenv(env_file)
else:
    print(f"警告: 未找到 .env 文件: {env_file}")

from app.database import SessionLocal
from app.models.result import SolvationStructure
from app.models.job import MDJob, JobStatus
from app.services.solvation import analyze_solvation_structures
from app.models.electrolyte import ElectrolyteSystem
from sqlalchemy import func
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def refresh_solvation_for_md_job(db: SessionLocal, md_job_id: int, dry_run: bool = False):
    """
    重新生成指定 MD 任务的溶剂化结构
    
    Args:
        db: 数据库会话
        md_job_id: MD 任务 ID
        dry_run: 如果为 True，只打印信息，不实际修改数据库
    """
    # 获取 MD 任务
    job = db.query(MDJob).filter(MDJob.id == md_job_id).first()
    if not job:
        logger.error(f"MD 任务 {md_job_id} 不存在")
        return False
    
    # 检查任务状态
    if job.status not in [JobStatus.COMPLETED, JobStatus.POSTPROCESSING]:
        logger.warning(f"MD 任务 {md_job_id} 状态为 {job.status}，跳过")
        return False
    
    # 检查工作目录
    if not job.work_dir:
        logger.warning(f"MD 任务 {md_job_id} 没有工作目录，跳过")
        return False
    
    work_dir = Path(job.work_dir)
    if not work_dir.exists():
        logger.warning(f"MD 任务 {md_job_id} 工作目录不存在: {work_dir}，跳过")
        return False
    
    # 获取电解液配置
    snapshot = job.config.get("system_snapshot") if job.config else None
    
    if snapshot:
        electrolyte_data = {
            "cations": snapshot.get("cations"),
            "anions": snapshot.get("anions"),
            "solvents": snapshot.get("solvents"),
        }
    else:
        # 兼容旧任务
        electrolyte = db.query(ElectrolyteSystem).filter(
            ElectrolyteSystem.id == job.system_id
        ).first()
        
        if not electrolyte:
            logger.warning(f"MD 任务 {md_job_id} 没有电解液配置，跳过")
            return False
        
        electrolyte_data = {
            "cations": electrolyte.cations,
            "anions": electrolyte.anions,
            "solvents": electrolyte.solvents,
        }
    
    logger.info(f"开始重新分析 MD 任务 {md_job_id} 的溶剂化结构...")
    
    if dry_run:
        logger.info(f"  [DRY RUN] 将重新分析工作目录: {work_dir}")
        return True
    
    try:
        # 执行溶剂化分析
        results = analyze_solvation_structures(
            work_dir=str(work_dir),
            electrolyte_data=electrolyte_data,
            cutoff=3.0,
        )
        
        if not results:
            logger.warning(f"MD 任务 {md_job_id} 未找到溶剂化结构数据")
            return False
        
        # 删除旧的溶剂化结构记录
        old_count = db.query(SolvationStructure).filter(
            SolvationStructure.md_job_id == md_job_id
        ).count()
        
        db.query(SolvationStructure).filter(
            SolvationStructure.md_job_id == md_job_id
        ).delete()
        
        # 插入新的溶剂化结构记录
        for result in results:
            structure = SolvationStructure(
                md_job_id=md_job_id,
                center_ion=result['center_ion'],
                structure_type=result['structure_type'],
                coordination_num=result['coordination_num'],
                composition=result['composition'],
                mol_order=result.get('mol_order'),
                file_path=result['file_path'],
                xyz_content=result.get('xyz_content'),
                snapshot_frame=result['snapshot_frame'],
                description=result['description'],
            )
            db.add(structure)
        
        db.commit()
        
        logger.info(f"✅ MD 任务 {md_job_id} 溶剂化结构更新成功: {old_count} → {len(results)} 个结构")
        
        # 检查 mol_order
        new_structures = db.query(SolvationStructure).filter(
            SolvationStructure.md_job_id == md_job_id
        ).all()
        
        with_mol_order = sum(1 for s in new_structures if s.mol_order)
        logger.info(f"   包含 mol_order 的结构: {with_mol_order}/{len(new_structures)}")
        
        return True
    
    except Exception as e:
        logger.error(f"❌ MD 任务 {md_job_id} 溶剂化结构分析失败: {e}")
        db.rollback()
        return False


def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(description='重新生成旧的溶剂化结构')
    parser.add_argument('--dry-run', action='store_true', help='只打印信息，不实际修改数据库')
    parser.add_argument('--md-job-id', type=int, help='只处理指定的 MD 任务 ID')
    args = parser.parse_args()
    
    db = SessionLocal()
    
    try:
        if args.md_job_id:
            # 只处理指定的 MD 任务
            logger.info(f"处理 MD 任务 {args.md_job_id}")
            success = refresh_solvation_for_md_job(db, args.md_job_id, dry_run=args.dry_run)
            if success:
                logger.info("✅ 完成")
            else:
                logger.error("❌ 失败")
        else:
            # 查找所有需要更新的 MD 任务
            # 找到所有有溶剂化结构但 mol_order 为 None 的 MD 任务
            md_jobs_with_old_structures = db.query(SolvationStructure.md_job_id).filter(
                SolvationStructure.mol_order.is_(None)
            ).distinct().all()
            
            md_job_ids = [job_id for (job_id,) in md_jobs_with_old_structures]
            
            logger.info(f"找到 {len(md_job_ids)} 个需要更新的 MD 任务")
            
            if args.dry_run:
                logger.info("DRY RUN 模式，不会实际修改数据库")
            
            success_count = 0
            fail_count = 0
            skip_count = 0
            
            for i, md_job_id in enumerate(md_job_ids, 1):
                logger.info(f"\n[{i}/{len(md_job_ids)}] 处理 MD 任务 {md_job_id}")
                
                result = refresh_solvation_for_md_job(db, md_job_id, dry_run=args.dry_run)
                
                if result:
                    success_count += 1
                elif result is False:
                    fail_count += 1
                else:
                    skip_count += 1
            
            logger.info(f"\n=== 总结 ===")
            logger.info(f"成功: {success_count}")
            logger.info(f"失败: {fail_count}")
            logger.info(f"跳过: {skip_count}")
            logger.info(f"总计: {len(md_job_ids)}")
    
    finally:
        db.close()


if __name__ == '__main__':
    main()

