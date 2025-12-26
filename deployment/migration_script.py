"""
数据存储迁移脚本

用途: 将旧的扁平目录结构迁移到新的层次化结构

使用方法:
    python migration_script.py --dry-run  # 预览迁移
    python migration_script.py --execute  # 执行迁移
    python migration_script.py --rollback # 回滚迁移

警告: 执行前请务必备份数据!
"""
import shutil
import logging
from pathlib import Path
from datetime import datetime
from typing import List, Dict
import argparse
import json

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class DataMigration:
    """数据迁移工具"""
    
    def __init__(self, base_dir: str = "/public/home/xiaoji/molyte_web"):
        self.base_dir = Path(base_dir)
        self.old_data_dir = self.base_dir / "data"
        self.new_data_dir = self.base_dir / "data_new"
        self.backup_dir = self.base_dir / "data_backup"
        self.migration_log = []
    
    def create_new_structure(self):
        """创建新的目录结构"""
        logger.info("Creating new directory structure...")
        
        directories = [
            # Jobs
            "jobs/qc/active",
            "jobs/qc/completed",
            "jobs/qc/archived",
            "jobs/md/active",
            "jobs/md/completed",
            "jobs/md/archived",
            "jobs/anion/active",
            "jobs/anion/completed",
            "jobs/anion/archived",
            "jobs/cluster/active",
            "jobs/cluster/completed",
            "jobs/cluster/archived",
            # Library
            "library/force_fields/cations",
            "library/force_fields/anions",
            "library/charges/anions",
            "library/charges/solvents",
            "library/molecules/common",
            # Users
            "users",
            # Temp
            "temp/uploads_pending",
            "temp/coord_gen",
            "temp/calculations",
            # Backup
            "backup/database",
            "backup/critical_jobs",
            # Logs
            "logs/api",
            "logs/workers",
            "logs/slurm",
        ]
        
        for dir_path in directories:
            full_path = self.new_data_dir / dir_path
            full_path.mkdir(parents=True, exist_ok=True)
            logger.info(f"✓ Created: {full_path}")
    
    def backup_current_data(self):
        """备份当前数据"""
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        backup_path = self.backup_dir / f"pre_migration_{timestamp}"
        
        logger.info(f"Backing up data to {backup_path}...")
        
        try:
            shutil.copytree(self.old_data_dir, backup_path, symlinks=True)
            logger.info(f"✓ Backup completed: {backup_path}")
            return backup_path
        except Exception as e:
            logger.error(f"✗ Backup failed: {e}")
            return None
    
    def migrate_qc_jobs(self, dry_run: bool = True):
        """迁移QC任务"""
        logger.info("Migrating QC jobs...")
        
        old_qc_work = self.old_data_dir / "qc_work"
        if not old_qc_work.exists():
            logger.warning(f"Old QC work directory not found: {old_qc_work}")
            return
        
        # 查找所有QC任务目录
        qc_jobs = list(old_qc_work.glob("QC-*"))
        logger.info(f"Found {len(qc_jobs)} QC jobs to migrate")
        
        for old_job_dir in qc_jobs:
            job_id = old_job_dir.name
            
            # 创建新的任务目录结构
            new_job_dir = self.new_data_dir / "jobs" / "qc" / "active" / job_id
            
            if dry_run:
                logger.info(f"[DRY RUN] Would migrate: {old_job_dir} → {new_job_dir}")
                continue
            
            # 创建子目录
            input_dir = new_job_dir / "input"
            work_dir = new_job_dir / "work"
            output_dir = new_job_dir / "output"
            logs_dir = new_job_dir / "logs"
            
            for d in [input_dir, work_dir, output_dir, logs_dir]:
                d.mkdir(parents=True, exist_ok=True)
            
            # 移动文件到对应目录
            try:
                # Input files (.gjf)
                for file in old_job_dir.glob("*.gjf"):
                    shutil.move(str(file), str(input_dir / file.name))
                
                # Output files (.log, .fchk)
                for file in old_job_dir.glob("*.log"):
                    if not file.name.startswith("qc_"):  # 排除slurm日志
                        shutil.move(str(file), str(output_dir / file.name))
                
                for file in old_job_dir.glob("*.fchk"):
                    shutil.move(str(file), str(output_dir / file.name))
                
                # Work files (.chk)
                for file in old_job_dir.glob("*.chk"):
                    shutil.move(str(file), str(work_dir / file.name))
                
                # Logs (qc_*.log)
                for file in old_job_dir.glob("qc_*.log"):
                    shutil.move(str(file), str(logs_dir / file.name))
                
                # 移动剩余文件到work
                for file in old_job_dir.iterdir():
                    if file.is_file():
                        shutil.move(str(file), str(work_dir / file.name))
                
                self.migration_log.append({
                    'job_id': job_id,
                    'type': 'qc',
                    'status': 'success',
                    'timestamp': datetime.now().isoformat()
                })
                
                logger.info(f"✓ Migrated: {job_id}")
                
            except Exception as e:
                logger.error(f"✗ Error migrating {job_id}: {e}")
                self.migration_log.append({
                    'job_id': job_id,
                    'type': 'qc',
                    'status': 'error',
                    'error': str(e),
                    'timestamp': datetime.now().isoformat()
                })
    
    def migrate_md_jobs(self, dry_run: bool = True):
        """迁移MD任务"""
        logger.info("Migrating MD jobs...")
        
        old_md_work = self.old_data_dir / "md_work"
        if not old_md_work.exists():
            logger.warning(f"Old MD work directory not found: {old_md_work}")
            return
        
        md_jobs = list(old_md_work.glob("MD-*"))
        logger.info(f"Found {len(md_jobs)} MD jobs to migrate")
        
        for old_job_dir in md_jobs:
            job_id = old_job_dir.name
            new_job_dir = self.new_data_dir / "jobs" / "md" / "active" / job_id
            
            if dry_run:
                logger.info(f"[DRY RUN] Would migrate: {old_job_dir} → {new_job_dir}")
                continue
            
            try:
                # MD任务结构类似，创建子目录
                for subdir in ["input", "work", "output", "logs"]:
                    (new_job_dir / subdir).mkdir(parents=True, exist_ok=True)
                
                # 简单移动整个目录
                shutil.copytree(old_job_dir, new_job_dir / "work", dirs_exist_ok=True)
                
                logger.info(f"✓ Migrated: {job_id}")
                
            except Exception as e:
                logger.error(f"✗ Error migrating {job_id}: {e}")
    
    def migrate_library_files(self, dry_run: bool = True):
        """迁移库文件（力场、电荷等）"""
        logger.info("Migrating library files...")
        
        # 迁移initial_salts到library/force_fields
        old_salts = self.old_data_dir / "initial_salts"
        if old_salts.exists():
            new_ff_dir = self.new_data_dir / "library" / "force_fields"
            
            if dry_run:
                logger.info(f"[DRY RUN] Would migrate: {old_salts} → {new_ff_dir}")
            else:
                try:
                    shutil.copytree(old_salts, new_ff_dir, dirs_exist_ok=True)
                    logger.info(f"✓ Migrated force fields")
                except Exception as e:
                    logger.error(f"✗ Error migrating force fields: {e}")
        
        # 迁移charges
        old_charges = self.old_data_dir / "charges"
        if old_charges.exists():
            new_charges_dir = self.new_data_dir / "library" / "charges"
            
            if dry_run:
                logger.info(f"[DRY RUN] Would migrate: {old_charges} → {new_charges_dir}")
            else:
                try:
                    shutil.copytree(old_charges, new_charges_dir, dirs_exist_ok=True)
                    logger.info(f"✓ Migrated charges")
                except Exception as e:
                    logger.error(f"✗ Error migrating charges: {e}")
    
    def create_symlinks(self):
        """创建软链接以保持向后兼容"""
        logger.info("Creating compatibility symlinks...")
        
        symlinks = [
            (self.old_data_dir / "qc_work", self.new_data_dir / "jobs" / "qc" / "active"),
            (self.old_data_dir / "md_work", self.new_data_dir / "jobs" / "md" / "active"),
            (self.old_data_dir / "initial_salts", self.new_data_dir / "library" / "force_fields"),
            (self.old_data_dir / "charges", self.new_data_dir / "library" / "charges"),
        ]
        
        for link_path, target_path in symlinks:
            try:
                if link_path.exists() and not link_path.is_symlink():
                    # 重命名原目录
                    link_path.rename(link_path.with_suffix('.old'))
                
                if not link_path.exists():
                    link_path.symlink_to(target_path)
                    logger.info(f"✓ Created symlink: {link_path} → {target_path}")
            except Exception as e:
                logger.error(f"✗ Error creating symlink: {e}")
    
    def save_migration_log(self):
        """保存迁移日志"""
        log_file = self.backup_dir / f"migration_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        log_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(log_file, 'w') as f:
            json.dump(self.migration_log, f, indent=2)
        
        logger.info(f"Migration log saved to: {log_file}")
    
    def run_migration(self, dry_run: bool = True):
        """执行完整迁移"""
        logger.info("="* 60)
        logger.info("Starting data migration...")
        logger.info(f"Mode: {'DRY RUN' if dry_run else 'EXECUTE'}")
        logger.info("="* 60)
        
        # 1. 创建新结构
        self.create_new_structure()
        
        if not dry_run:
            # 2. 备份
            backup_path = self.backup_current_data()
            if not backup_path:
                logger.error("Backup failed! Aborting migration.")
                return False
        
        # 3. 迁移数据
        self.migrate_qc_jobs(dry_run)
        self.migrate_md_jobs(dry_run)
        self.migrate_library_files(dry_run)
        
        if not dry_run:
            # 4. 创建软链接
            self.create_symlinks()
            
            # 5. 保存日志
            self.save_migration_log()
        
        logger.info("="* 60)
        logger.info("Migration completed!")
        logger.info("="* 60)
        
        return True


def main():
    parser = argparse.ArgumentParser(description="Migrate data to new directory structure")
    parser.add_argument(
        '--mode',
        choices=['dry-run', 'execute', 'rollback'],
        default='dry-run',
        help='Migration mode'
    )
    parser.add_argument(
        '--base-dir',
        default='/public/home/xiaoji/molyte_web',
        help='Base directory path'
    )
    
    args = parser.parse_args()
    
    migration = DataMigration(base_dir=args.base_dir)
    
    if args.mode == 'dry-run':
        migration.run_migration(dry_run=True)
    elif args.mode == 'execute':
        print("⚠️  WARNING: This will migrate your data!")
        print("⚠️  Make sure you have a backup!")
        response = input("Continue? (yes/no): ")
        if response.lower() == 'yes':
            migration.run_migration(dry_run=False)
        else:
            print("Migration cancelled.")
    elif args.mode == 'rollback':
        print("Rollback functionality not yet implemented")
        print("Please manually restore from backup")


if __name__ == "__main__":
    main()
