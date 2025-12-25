"""
数据清理和归档自动化工具

功能:
- 自动清理completed目录中的过期任务
- 归档任务到压缩文件
- 磁盘使用监控
- 任务统计

使用方法:
    python cleanup_automation.py --check    # 检查可清理的任务
    python cleanup_automation.py --clean    # 执行清理
    python cleanup_automation.py --archive  # 归档任务
"""
import shutil
import tarfile
import logging
from pathlib import Path
from datetime import datetime, timedelta
from typing import Dict, List
import argparse

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class DataCleanup:
    """数据清理和归档工具"""
    
    def __init__(self, base_dir: str = "/public/home/xiaoji/molyte_web"):
        self.base_dir = Path(base_dir)
        self.data_dir = self.base_dir / "data"
        self.jobs_dir = self.data_dir / "jobs"
    
    def find_completedtasks(self, days_old: int = 7) -> List[Path]:
        """
        查找超过N天的completed任务
        
        Args:
            days_old: 保留天数
            
        Returns:
            需要清理的任务目录列表
        """
        cutoff_date = datetime.now() - timedelta(days=days_old)
        tasks_to_clean = []
        
        for job_type in ['qc', 'md', 'anion', 'cluster']:
            completed_dir = self.jobs_dir / job_type / "completed"
            
            if not completed_dir.exists():
                continue
            
            # 遍历日期目录
            for date_dir in completed_dir.glob("20*"):  # 2025-12-24 格式
                try:
                    # 解析日期
                    dir_date = datetime.strptime(date_dir.name, '%Y-%m-%d')
                    
                    if dir_date < cutoff_date:
                        # 遍历该日期下的所有任务
                        for task_dir in date_dir.iterdir():
                            if task_dir.is_dir():
                                tasks_to_clean.append(task_dir)
                except (ValueError, OSError) as e:
                    logger.warning(f"Skipping invalid directory {date_dir}: {e}")
                    continue
        
        return tasks_to_clean
    
    def archive_task(self, task_dir: Path, job_type: str) -> Path:
        """
        归档单个任务
        
        Args:
            task_dir: 任务目录
            job_type: 任务类型
            
        Returns:
            归档文件路径
        """
        # 确定归档月份
        date_dir = task_dir.parent
        month = date_dir.name[:7]  # 2025-12
        
        # 归档目录
        archive_dir = self.jobs_dir / job_type / "archived" / month
        archive_dir.mkdir(parents=True, exist_ok=True)
        
        # 归档文件名
        task_id = task_dir.name
        archive_file = archive_dir / f"{date_dir.name}_{task_id}.tar.gz"
        
        # 创建tar.gz归档
        with tarfile.open(archive_file, "w:gz") as tar:
            tar.add(task_dir, arcname=task_id)
        
        logger.info(f"✓ Archived: {task_dir} → {archive_file}")
        return archive_file
    
    def cleanup_completed(self, days_to_keep: int = 7, dry_run: bool = True):
        """
        清理completed目录
        
        Args:
            days_to_keep: 保留天数
            dry_run: 是否只是预览
        """
        logger.info(f"{'[DRY RUN] ' if dry_run else ''}Cleaning completed tasks older than {days_to_keep} days...")
        
        tasks_to_clean = self.find_completed_tasks(days_old=days_to_keep)
        
        if not tasks_to_clean:
            logger.info("No tasks to clean")
            return
        
        logger.info(f"Found {len(tasks_to_clean)} tasks to clean")
        
        for task_dir in tasks_to_clean:
            # 确定任务类型
            job_type = task_dir.parts[-4]  # jobs/qc/completed/date/task
            task_id = task_dir.name
            
            if dry_run:
                size = self._get_dir_size(task_dir)
                logger.info(f"[DRY RUN] Would clean: {task_id} ({size / 1e6:.1f} MB)")
                continue
            
            try:
                # 归档
                archive_file = self.archive_task(task_dir, job_type)
                
                # 删除原目录
                shutil.rmtree(task_dir)
                
                logger.info(f"✓ Cleaned: {task_id}")
                
            except Exception as e:
                logger.error(f"✗ Error cleaning {task_id}: {e}")
        
        # 清理空的日期目录
        if not dry_run:
            self._cleanup_empty_date_dirs()
    
    def _cleanup_empty_date_dirs(self):
        """清理空的日期目录"""
        for job_type in ['qc', 'md', 'anion', 'cluster']:
            completed_dir = self.jobs_dir / job_type / "completed"
            
            if not completed_dir.exists():
                continue
            
            for date_dir in completed_dir.iterdir():
                if date_dir.is_dir() and not list(date_dir.iterdir()):
                    date_dir.rmdir()
                    logger.info(f"✓ Removed empty directory: {date_dir}")
    
    def _get_dir_size(self, path: Path) -> int:
        """获取目录大小（字节）"""
        total_size = 0
        for item in path.rglob('*'):
            if item.is_file():
                total_size += item.stat().st_size
        return total_size
    
    def monitor_disk_usage(self) -> Dict[str, dict]:
        """
        监控磁盘使用情况
        
        Returns:
            各目录的使用统计
        """
        stats = {}
        
        directories_to_monitor = {
            'QC Active': self.jobs_dir / "qc" / "active",
            'QC Completed': self.jobs_dir / "qc" / "completed",
            'QC Archived': self.jobs_dir / "qc" / "archived",
            'MD Active': self.jobs_dir / "md" / "active",
            'MD Completed': self.jobs_dir / "md" / "completed",
            'Library': self.data_dir / "library",
            'Temp': self.data_dir / "temp",
        }
        
        for name, path in directories_to_monitor.items():
            if not path.exists():
                stats[name] = {'exists': False}
                continue
            
            # 磁盘使用
            total, used, free = shutil.disk_usage(path)
            
            # 目录大小
            dir_size = self._get_dir_size(path)
            
            # 文件/目录数量
            file_count = sum(1 for _ in path.rglob('*') if _.is_file())
            dir_count = sum(1 for _ in path.rglob('*') if _.is_dir())
            
            stats[name] = {
                'exists': True,
                'path': str(path),
                'size_gb': dir_size / 1e9,
                'file_count': file_count,
                'dir_count': dir_count,
                'disk_total_gb': total / 1e9,
                'disk_used_gb': used / 1e9,
                'disk_free_gb': free / 1e9,
                'disk_usage_pct': (used / total) * 100,
            }
        
        return stats
    
    def print_disk_usage_report(self):
        """打印磁盘使用报告"""
        stats = self.monitor_disk_usage()
        
        print("="* 80)
        print("Disk Usage Report")
        print("="* 80)
        print(f"{'Directory':<20} {'Size':<12} {'Files':<10} {'Dirs':<10} {'Disk %':<10}")
        print("-"* 80)
        
        for name, info in stats.items():
            if not info.get('exists'):
                print(f"{name:<20} {'N/A':<12} {'-':<10} {'-':<10} {'-':<10}")
                continue
            
            size_str = f"{info['size_gb']:.2f} GB"
            files_str = str(info['file_count'])
            dirs_str = str(info['dir_count'])
            disk_pct = f"{info['disk_usage_pct']:.1f}%"
            
            print(f"{name:<20} {size_str:<12} {files_str:<10} {dirs_str:<10} {disk_pct:<10}")
        
        print("="* 80)
    
    def count_jobs_by_status(self) -> Dict[str, int]:
        """统计各状态的任务数量"""
        counts = {}
        
        for job_type in ['qc', 'md', 'anion', 'cluster']:
            for status in ['active', 'completed', 'archived']:
                dir_path = self.jobs_dir / job_type / status
                
                if not dir_path.exists():
                    counts[f"{job_type}_{status}"] = 0
                    continue
                
                if status == 'archived':
                    # 归档文件数量
                    count = len(list(dir_path.rglob("*.tar.gz")))
                else:
                    # 目录数量
                    if status == 'completed':
                        # 遍历日期子目录
                        count = 0
                        for date_dir in dir_path.iterdir():
                            if date_dir.is_dir():
                                count += len([d for d in date_dir.iterdir() if d.is_dir()])
                    else:
                        count = len([d for d in dir_path.iterdir() if d.is_dir()])
                
                counts[f"{job_type}_{status}"] = count
        
        return counts
    
    def print_job_statistics(self):
        """打印任务统计"""
        counts = self.count_jobs_by_status()
        
        print("="* 60)
        print("Job Statistics")
        print("="* 60)
        print(f"{'Type':<15} {'Active':<10} {'Completed':<12} {'Archived':<10} {'Total':<10}")
        print("-"* 60)
        
        for job_type in ['qc', 'md', 'anion', 'cluster']:
            active = counts.get(f"{job_type}_active", 0)
            completed = counts.get(f"{job_type}_completed", 0)
            archived = counts.get(f"{job_type}_archived", 0)
            total = active + completed + archived
            
            print(f"{job_type.upper():<15} {active:<10} {completed:<12} {archived:<10} {total:<10}")
        
        print("="* 60)


def main():
    parser = argparse.ArgumentParser(description="Data cleanup and monitoring tool")
    parser.add_argument(
        '--action',
        choices=['check', 'clean', 'archive', 'stats', 'monitor'],
        required=True,
        help='Action to perform'
    )
    parser.add_argument(
        '--days',
        type=int,
        default=7,
        help='Days to keep completed tasks (default: 7)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Preview without making changes'
    )
    parser.add_argument(
        '--base-dir',
        default='/public/home/xiaoji/molyte_web',
        help='Base directory path'
    )
    
    args = parser.parse_args()
    
    cleanup = DataCleanup(base_dir=args.base_dir)
    
    if args.action == 'check':
        tasks = cleanup.find_completed_tasks(days_old=args.days)
        print(f"Found {len(tasks)} tasks to clean (older than {args.days} days)")
        for task in tasks[:10]:  # 只显示前10个
            print(f"  - {task}")
        if len(tasks) > 10:
            print(f"  ... and {len(tasks) - 10} more")
    
    elif args.action == 'clean':
        cleanup.cleanup_completed(days_to_keep=args.days, dry_run=args.dry_run)
    
    elif args.action == 'archive':
        print("Archive functionality is included in clean action")
    
    elif args.action == 'stats':
        cleanup.print_job_statistics()
    
    elif args.action == 'monitor':
        cleanup.print_disk_usage_report()
        print()
        cleanup.print_job_statistics()


if __name__ == "__main__":
    main()
