#!/usr/bin/env python3
"""
清理缺失文件的阴离子记录

这个脚本会检查数据库中的阴离子记录，如果对应的文件不存在，
则将记录标记为已删除（软删除）。
"""

import os
import sys
from datetime import datetime, timezone
from pathlib import Path

# 添加项目根目录到 Python 路径
sys.path.append('/opt/molyte_web_v1.0/backend')

from app.models.forcefield import AnionLibrary
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine

def cleanup_missing_anion_files():
    """清理缺失文件的阴离子记录"""
    
    # 数据库连接
    engine = create_engine('postgresql://molyte:molyte2025@127.0.0.1:5432/molyte_db')
    db = sessionmaker(bind=engine)()
    
    try:
        # 查找所有未删除的阴离子记录
        anions = db.query(AnionLibrary).filter(AnionLibrary.is_deleted == False).all()
        
        print(f"检查 {len(anions)} 个阴离子记录...")
        print("=" * 80)
        
        missing_files = []
        
        for anion in anions:
            # 检查文件是否存在
            lt_exists = os.path.exists(anion.lt_path) if anion.lt_path else False
            pdb_exists = os.path.exists(anion.pdb_path) if anion.pdb_path else False
            
            status = "✓" if (lt_exists and pdb_exists) else "✗"
            print(f"{status} {anion.anion_name}: lt={lt_exists}, pdb={pdb_exists}")
            
            # 如果文件不存在，标记为需要清理
            if not (lt_exists and pdb_exists):
                missing_files.append(anion)
        
        print("=" * 80)
        print(f"发现 {len(missing_files)} 个记录的文件缺失")
        
        if missing_files:
            print("\n准备清理以下记录:")
            for anion in missing_files:
                print(f"  - {anion.anion_name} (ID: {anion.id})")
            
            # 确认清理
            confirm = input(f"\n是否要软删除这 {len(missing_files)} 个记录? (y/N): ")
            
            if confirm.lower() == 'y':
                # 执行软删除
                for anion in missing_files:
                    anion.is_deleted = True
                    anion.deleted_at = datetime.now(timezone.utc)
                
                db.commit()
                print(f"✓ 已软删除 {len(missing_files)} 个记录")
            else:
                print("取消清理操作")
        else:
            print("所有记录的文件都存在，无需清理")
            
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()
        db.rollback()
    finally:
        db.close()

if __name__ == "__main__":
    cleanup_missing_anion_files()
