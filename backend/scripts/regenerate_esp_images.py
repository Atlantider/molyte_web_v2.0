#!/usr/bin/env python3
"""
批量重新生成 QC 任务的 ESP 图像并上传到 COS

使用方法：
    python regenerate_esp_images.py [--job-ids 1,2,3] [--all] [--upload]
"""

import os
import sys
import argparse
import logging
from pathlib import Path
from datetime import datetime
import base64
import yaml

# 添加项目路径
backend_path = Path(__file__).parent.parent
sys.path.insert(0, str(backend_path))

# 设置 PYTHONPATH 环境变量
os.environ['PYTHONPATH'] = str(backend_path)

from app.database import SessionLocal
from app.models.qc import QCJob, QCResult
from sqlalchemy import and_

# 延迟导入 generate_esp_visualization，避免循环导入
def get_generate_esp_visualization():
    from app.tasks.qc_postprocess import generate_esp_visualization
    return generate_esp_visualization

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)s | %(name)s | %(message)s'
)
logger = logging.getLogger(__name__)


def init_cos_client():
    """初始化 COS 客户端"""
    try:
        from qcloud_cos import CosConfig, CosS3Client

        # 加载配置
        config_path = Path(__file__).parent.parent.parent / "deployment" / "polling_worker_config.yaml"
        if not config_path.exists():
            logger.warning(f"Config file not found: {config_path}")
            return None, None, None

        with open(config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)

        if 'cos' not in config:
            logger.warning("COS config not found in polling_worker_config.yaml")
            return None, None, None

        cos_config = config['cos']
        cos_cfg = CosConfig(
            Region=cos_config['region'],
            SecretId=cos_config['secret_id'],
            SecretKey=cos_config['secret_key'],
            Scheme='https'
        )
        cos_client = CosS3Client(cos_cfg)
        bucket = cos_config['bucket']

        logger.info(f"COS client initialized (Bucket: {bucket})")
        return cos_client, bucket, config
    except Exception as e:
        logger.error(f"Failed to initialize COS client: {e}")
        return None, None, None


def regenerate_esp_for_job(job_id: int, db, cos_client=None, cos_bucket=None, config=None):
    """为单个 QC 任务重新生成 ESP 图像"""
    try:
        # 获取 QC 任务
        job = db.query(QCJob).filter(QCJob.id == job_id).first()
        if not job:
            logger.error(f"QC Job {job_id} not found")
            return False

        # 检查任务状态（使用字符串比较避免枚举导入问题）
        job_status = str(job.status) if hasattr(job.status, 'value') else job.status
        if job_status != "COMPLETED":
            logger.warning(f"QC Job {job_id} status is {job_status}, skipping")
            return False

        # 检查工作目录
        work_dir = Path(job.work_dir)
        if not work_dir.exists():
            logger.error(f"Work directory not found: {work_dir}")
            return False

        # 查找 .fchk 文件
        fchk_files = list(work_dir.glob("*.fchk"))
        if not fchk_files:
            logger.error(f"No .fchk file found in {work_dir}")
            return False

        fchk_file = fchk_files[0]
        logger.info(f"Found fchk file: {fchk_file.name}")

        # 生成 ESP 图像
        logger.info(f"Regenerating ESP image for job {job_id}...")
        generate_esp_visualization = get_generate_esp_visualization()
        esp_image_path = generate_esp_visualization(work_dir, job.molecule_name, str(fchk_file))

        if not esp_image_path:
            logger.error(f"Failed to generate ESP image for job {job_id}")
            return False

        logger.info(f"✅ ESP image generated: {esp_image_path}")

        # 获取或创建 QC 结果记录
        qc_result = db.query(QCResult).filter(QCResult.qc_job_id == job_id).first()
        if not qc_result:
            logger.warning(f"No QC result found for job {job_id}, creating new one")
            qc_result = QCResult(qc_job_id=job_id, smiles=job.smiles)
            db.add(qc_result)

        # 更新 ESP 图像路径（本地路径）
        qc_result.esp_image_path = str(esp_image_path)

        # 如果需要上传到 COS
        if cos_client and cos_bucket and config:
            try:
                # 获取 QC 结果前缀
                qc_result_prefix = config.get('cos', {}).get('qc_result_prefix', 'QC_results/')
                cos_key = f"{qc_result_prefix}{job_id}/ESP.png"

                # 上传到 COS
                with open(esp_image_path, 'rb') as f:
                    cos_client.put_object(
                        Bucket=cos_bucket,
                        Body=f,
                        Key=cos_key
                    )
                logger.info(f"✅ Uploaded to COS: {cos_key}")

                # 更新数据库中的 COS 路径
                qc_result.esp_image_path = cos_key
            except Exception as e:
                logger.warning(f"Failed to upload to COS: {e}")

        # 保存 base64 编码的图像内容（用于混合云架构）
        try:
            with open(esp_image_path, 'rb') as f:
                image_content = base64.b64encode(f.read()).decode('utf-8')
                qc_result.esp_image_content = image_content
                logger.info(f"Saved base64 encoded ESP image content")
        except Exception as e:
            logger.warning(f"Failed to encode ESP image: {e}")

        db.commit()
        logger.info(f"✅ Updated QC result for job {job_id}")
        return True

    except Exception as e:
        logger.error(f"Error regenerating ESP for job {job_id}: {e}", exc_info=True)
        return False


def main():
    print("Parsing arguments...", flush=True)
    parser = argparse.ArgumentParser(description="Regenerate ESP images for QC jobs")
    parser.add_argument("--job-ids", type=str, help="Comma-separated job IDs (e.g., 1,2,3)")
    parser.add_argument("--all", action="store_true", help="Regenerate for all completed jobs")
    parser.add_argument("--upload", action="store_true", help="Upload to COS")
    parser.add_argument("--limit", type=int, default=None, help="Limit number of jobs to process")

    args = parser.parse_args()
    print(f"Arguments: {args}", flush=True)

    print("Connecting to database...", flush=True)
    db = SessionLocal()
    print("Database connected", flush=True)

    # 初始化 COS 客户端（如果需要上传）
    cos_client = None
    cos_bucket = None
    config = None
    if args.upload:
        cos_client, cos_bucket, config = init_cos_client()
        if not cos_client:
            logger.warning("COS upload requested but client initialization failed, will skip COS upload")

    try:
        # 确定要处理的任务 ID
        job_ids = []

        print("Determining job IDs...", flush=True)
        if args.job_ids:
            job_ids = [int(x.strip()) for x in args.job_ids.split(",")]
            print(f"Processing specified jobs: {job_ids}", flush=True)
        elif args.all:
            # 获取所有已完成的 QC 任务（使用字符串比较避免枚举导入问题）
            print("Querying completed QC jobs...", flush=True)
            try:
                completed_jobs = db.query(QCJob).filter(
                    QCJob.status == "COMPLETED"
                ).all()
                job_ids = [job.id for job in completed_jobs]
                print(f"Found {len(job_ids)} completed QC jobs", flush=True)
            except Exception as e:
                print(f"Error querying jobs: {e}", flush=True)
                raise
        else:
            logger.error("Please specify --job-ids or --all")
            return 1

        # 限制处理数量
        if args.limit:
            job_ids = job_ids[:args.limit]
            logger.info(f"Limited to {len(job_ids)} jobs")

        # 处理每个任务
        success_count = 0
        failed_count = 0

        for i, job_id in enumerate(job_ids, 1):
            logger.info(f"\n[{i}/{len(job_ids)}] Processing job {job_id}...")
            if regenerate_esp_for_job(job_id, db, cos_client, cos_bucket, config):
                success_count += 1
            else:
                failed_count += 1

        # 总结
        logger.info(f"\n{'='*60}")
        logger.info(f"✅ Successfully regenerated: {success_count}")
        logger.info(f"❌ Failed: {failed_count}")
        logger.info(f"{'='*60}")

        return 0 if failed_count == 0 else 1

    finally:
        db.close()


if __name__ == "__main__":
    print("Starting regenerate_esp_images.py...", flush=True)
    try:
        sys.exit(main())
    except Exception as e:
        print(f"Error: {e}", flush=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)

