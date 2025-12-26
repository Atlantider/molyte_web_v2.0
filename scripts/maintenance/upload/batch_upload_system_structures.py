#!/usr/bin/env python3
"""
批量上传系统结构脚本
自动发现所有已完成的MD任务并上传系统结构到数据库
"""

import os
import sys
import requests
import json
import glob
from pathlib import Path

# 配置
API_BASE_URL = "https://www.molyte.xyz/api/v1"
WORKER_TOKEN = "worker_secret_token_2024"

def discover_completed_md_tasks():
    """
    自动发现所有已完成的MD任务
    通过查找包含 after_nvt.lammpstrj 文件的目录
    """
    tasks = []
    md_work_dir = Path("data/md_work")
    
    # 查找所有 after_nvt.lammpstrj 文件
    lammpstrj_files = list(md_work_dir.glob("*/*after_nvt.lammpstrj"))
    
    print(f"🔍 发现 {len(lammpstrj_files)} 个包含 after_nvt.lammpstrj 的任务目录")
    
    for file_path in lammpstrj_files:
        job_dir = file_path.parent
        job_name = job_dir.name
        
        # 提取任务信息
        task_info = {
            'job_name': job_name,
            'work_dir': str(job_dir),
            'lammpstrj_file': str(file_path),
            'task_id': None  # 稍后从API获取
        }
        
        tasks.append(task_info)
        print(f"  📁 {job_name}")
    
    return tasks

def get_job_id_by_name(job_name):
    """
    通过任务名称从API获取任务ID
    """
    try:
        # 尝试从API获取任务信息
        response = requests.get(
            f"{API_BASE_URL}/jobs",
            headers={"Authorization": f"Bearer {WORKER_TOKEN}"},
            params={"limit": 1000}  # 获取更多任务
        )
        
        if response.status_code == 200:
            jobs = response.json()
            for job in jobs:
                if job.get('name') == job_name:
                    return job.get('id')
        
        print(f"⚠️  无法找到任务 {job_name} 的ID")
        return None
        
    except Exception as e:
        print(f"❌ 获取任务ID失败: {e}")
        return None

def extract_system_structure_from_lammpstrj(lammpstrj_file, job_name):
    """
    从 LAMMPS 轨迹文件中提取系统结构
    """
    try:
        with open(lammpstrj_file, 'r') as f:
            lines = f.readlines()
        
        # 解析 LAMMPS 轨迹文件
        atoms = []
        in_atoms_section = False
        current_timestep = None
        atom_count = 0
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            if line.startswith("ITEM: TIMESTEP"):
                current_timestep = int(lines[i+1].strip())
                i += 2
                continue
            elif line.startswith("ITEM: NUMBER OF ATOMS"):
                atom_count = int(lines[i+1].strip())
                i += 2
                continue
            elif line.startswith("ITEM: ATOMS"):
                # 开始原子数据部分
                in_atoms_section = True
                atoms = []  # 重置原子列表（取最后一帧）
                i += 1
                continue
            elif in_atoms_section and line:
                # 解析原子行: id element mol type x y z q
                parts = line.split()
                if len(parts) >= 8:
                    atom_id = int(parts[0])
                    element = parts[1]
                    x, y, z = float(parts[4]), float(parts[5]), float(parts[6])
                    atoms.append((atom_id, element, x, y, z))
            
            i += 1
        
        if not atoms:
            print(f"❌ 无法从 {lammpstrj_file} 提取原子数据")
            return None
        
        # 转换为 XYZ 格式
        xyz_lines = [
            str(len(atoms)),
            f"Frame from {job_name} at timestep {current_timestep}"
        ]
        
        # 按原子ID排序
        for atom_id, element, x, y, z in sorted(atoms):
            xyz_lines.append(f"{element} {x:.6f} {y:.6f} {z:.6f}")
        
        xyz_content = "\n".join(xyz_lines)
        
        print(f"✅ 成功提取系统结构: {len(atoms)} 个原子, 时间步 {current_timestep}")
        return {
            'xyz_content': xyz_content,
            'atom_count': len(atoms),
            'timestep': current_timestep
        }
        
    except Exception as e:
        print(f"❌ 提取系统结构失败: {e}")
        return None

def upload_system_structure(job_id, structure_data):
    """
    上传系统结构到API
    """
    try:
        url = f"{API_BASE_URL}/workers/jobs/{job_id}/system_structure"
        
        payload = {
            'xyz_content': structure_data['xyz_content'],
            'atom_count': structure_data['atom_count'],
            'frame_info': f"Timestep {structure_data['timestep']}"
        }
        
        response = requests.post(
            url,
            json=payload,
            headers={
                "Authorization": f"Bearer {WORKER_TOKEN}",
                "Content-Type": "application/json"
            }
        )
        
        if response.status_code == 200:
            print(f"✅ 系统结构上传成功")
            return True
        else:
            print(f"❌ 系统结构上传失败: HTTP {response.status_code}")
            print(f"   响应: {response.text}")
            return False
            
    except Exception as e:
        print(f"❌ 系统结构上传失败: {e}")
        return False

def main():
    """
    主函数
    """
    print("🚀 开始批量上传系统结构...")
    
    # 发现所有已完成的MD任务
    tasks = discover_completed_md_tasks()
    
    if not tasks:
        print("❌ 没有发现任何已完成的MD任务")
        return
    
    print(f"\n📊 总共发现 {len(tasks)} 个已完成的MD任务")
    
    success_count = 0
    failed_count = 0
    
    for i, task in enumerate(tasks, 1):
        print(f"\n{'='*60}")
        print(f"🔄 处理任务 {i}/{len(tasks)}: {task['job_name']}")
        
        # 获取任务ID
        job_id = get_job_id_by_name(task['job_name'])
        if not job_id:
            print(f"❌ 跳过任务 {task['job_name']}: 无法获取任务ID")
            failed_count += 1
            continue
        
        print(f"📊 任务ID: {job_id}")
        
        # 提取系统结构
        structure_data = extract_system_structure_from_lammpstrj(
            task['lammpstrj_file'], 
            task['job_name']
        )
        
        if not structure_data:
            print(f"❌ 跳过任务 {task['job_name']}: 无法提取系统结构")
            failed_count += 1
            continue
        
        # 上传系统结构
        if upload_system_structure(job_id, structure_data):
            print(f"✅ 任务 {task['job_name']} 处理成功")
            success_count += 1
        else:
            print(f"❌ 任务 {task['job_name']} 处理失败")
            failed_count += 1
    
    print(f"\n{'='*60}")
    print(f"🎉 批量处理完成!")
    print(f"   ✅ 成功处理: {success_count} 个任务")
    print(f"   ❌ 处理失败: {failed_count} 个任务")
    print(f"   📊 总计: {len(tasks)} 个任务")
    
    if success_count > 0:
        print(f"\n💡 建议刷新前端页面查看更新后的系统结构可视化")

if __name__ == "__main__":
    main()
