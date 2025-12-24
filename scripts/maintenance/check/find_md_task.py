#!/usr/bin/env python3
"""查找特定的MD任务"""
import sys
import os
import sqlite3
from pathlib import Path

def find_md_task():
    # 查找数据库文件
    db_path = Path('/public/home/xiaoji/molyte_web/backend/molyte.db')
    if not db_path.exists():
        print(f"数据库文件不存在: {db_path}")
        return []

    try:
        conn = sqlite3.connect(str(db_path))
        cursor = conn.cursor()

        # 查找包含 EL-20251205-0002 的任务
        query = """
        SELECT id, status, work_dir, config, created_at, updated_at, error_message, slurm_job_id
        FROM md_jobs
        WHERE work_dir LIKE '%EL-20251205-0002%'
           OR config LIKE '%EL-20251205-0002%'
        """

        cursor.execute(query)
        tasks = cursor.fetchall()

        print(f"找到 {len(tasks)} 个匹配的任务:")
        for task in tasks:
            task_id, status, work_dir, config, created_at, updated_at, error_message, slurm_job_id = task
            print(f"  ID: {task_id}")
            print(f"  状态: {status}")
            print(f"  工作目录: {work_dir}")
            print(f"  Slurm Job ID: {slurm_job_id}")
            print(f"  创建时间: {created_at}")
            print(f"  更新时间: {updated_at}")
            print(f"  错误信息: {error_message}")
            print("-" * 80)

        return tasks

    except Exception as e:
        print(f"查询数据库失败: {e}")
        return []
    finally:
        if 'conn' in locals():
            conn.close()

if __name__ == "__main__":
    tasks = find_md_task()
