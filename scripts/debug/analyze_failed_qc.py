#!/usr/bin/env python3
"""
分析失败的 QC 任务，查看 Gaussian 日志中的具体错误
"""

import sys
import subprocess
from pathlib import Path

# 失败任务的 Slurm ID
failed_jobs = {
    '7651': 'intermediate_Li_PF6_1_struct1466',
    '7880': 'intermediate_Li_EC_1_struct1466',
    '7505': 'intermediate_Li_EC_1_PF6_1_struct1466',
    '7316': 'intermediate_Li_EC_1_DEC_1_PF6_1_struct1466',
    '7325': 'intermediate_Li_DEC_1_PF6_1_struct1466',
    '7719': 'intermediate_Li_PF6_3_struct1466',
}

def find_work_dir(slurm_id, job_name):
    """查找工作目录"""
    qc_work_base = Path('data/qc_work')
    
    # 尝试多种可能的目录名
    possible_dirs = [
        qc_work_base / f"QC-{slurm_id}",
        qc_work_base / f"QC-{slurm_id}-{job_name}",
    ]
    
    # 也搜索包含 job_name 的目录
    if qc_work_base.exists():
        for d in qc_work_base.iterdir():
            if d.is_dir() and (slurm_id in d.name or job_name in d.name):
                possible_dirs.append(d)
    
    for d in possible_dirs:
        if d.exists():
            return d
    
    return None

def analyze_gaussian_log(log_file):
    """分析 Gaussian 日志文件"""
    if not log_file.exists():
        return "日志文件不存在"
    
    with open(log_file, 'r', errors='ignore') as f:
        content = f.read()
    
    errors = []
    
    # 检查是否正常结束
    if 'Normal termination' in content:
        errors.append("✅ 正常结束")
    elif 'Error termination' in content:
        errors.append("❌ 异常终止")
        
        # 提取错误行
        for line in content.split('\n'):
            if 'Error' in line and len(line) < 200:
                errors.append(f"  错误: {line.strip()}")
    
    # 检查常见错误
    if 'combination of multiplicity' in content and 'electrons is impossible' in content:
        errors.append("❌ 电荷/自旋多重度错误：电子数与自旋多重度不匹配")
        # 提取具体信息
        for line in content.split('\n'):
            if 'combination of multiplicity' in line:
                errors.append(f"  详情: {line.strip()}")
    
    if 'Convergence failure' in content or 'SCF has not converged' in content:
        errors.append("❌ SCF 不收敛")
    
    if 'segmentation violation' in content or 'segmentation fault' in content:
        errors.append("❌ 段错误（Gaussian 崩溃）")
    
    if 'Tors failed' in content or 'FormBX had a problem' in content:
        errors.append("❌ 内部坐标系统崩溃")
    
    if 'galloc' in content or 'Out of memory' in content:
        errors.append("❌ 内存不足")
    
    # 提取电荷和自旋多重度
    for line in content.split('\n'):
        if line.strip().startswith('Charge =') and 'Multiplicity =' in line:
            errors.append(f"  设置: {line.strip()}")
            break
    
    # 提取分子式
    for line in content.split('\n'):
        if 'Stoichiometry' in line:
            errors.append(f"  分子式: {line.strip()}")
            break
    
    return '\n'.join(errors) if errors else "未发现明显错误"

def main():
    print("=" * 80)
    print("分析失败的 QC 任务")
    print("=" * 80)
    
    for slurm_id, job_name in failed_jobs.items():
        print(f"\n{'='*80}")
        print(f"任务: {job_name}")
        print(f"Slurm ID: {slurm_id}")
        print(f"{'='*80}")
        
        # 查找工作目录
        work_dir = find_work_dir(slurm_id, job_name)
        
        if not work_dir:
            print(f"⚠️  未找到工作目录")
            continue
        
        print(f"工作目录: {work_dir}")
        
        # 查找 Gaussian 日志文件
        log_files = list(work_dir.glob("*_out.log")) + list(work_dir.glob("*.log"))
        log_files = [f for f in log_files if 'qc_out' not in f.name and 'qc_err' not in f.name]
        
        if not log_files:
            print(f"⚠️  未找到 Gaussian 日志文件")
            continue
        
        log_file = log_files[0]
        print(f"日志文件: {log_file.name}")
        print(f"\n分析结果:")
        print(analyze_gaussian_log(log_file))

if __name__ == '__main__':
    main()

