#!/usr/bin/env python3
"""
修复 PDB 文件格式，确保残基号在正确的列位置（23-26）
"""

from pathlib import Path
import re

def fix_pdb_file(pdb_path: Path) -> bool:
    """
    修复 PDB 文件格式，参考 FSI.pdb 的格式

    PDB 格式规范（参考 FSI.pdb）：
    - 列 1-6: 记录名称（ATOM）
    - 列 7-11: 原子序号
    - 列 13-16: 原子名称
    - 列 18-20: 残基名称
    - 列 23-26: 残基序号（整数，右对齐）
    - 列 31-38: X 坐标
    - 列 39-46: Y 坐标
    - 列 47-54: Z 坐标
    """
    try:
        with open(pdb_path, 'r') as f:
            lines = f.readlines()

        fixed_lines = []
        for line in lines:
            if line.startswith('ATOM'):
                # 解析现有的 ATOM 行
                parts = line.split()
                if len(parts) < 6:
                    fixed_lines.append(line)
                    continue

                # 提取信息
                atom_num = int(parts[1])
                atom_name = parts[2]
                res_name = parts[3]
                res_num = 1  # 默认残基号为 1

                # 提取坐标
                try:
                    x = float(parts[5])
                    y = float(parts[6])
                    z = float(parts[7])
                except (ValueError, IndexError):
                    # 如果无法解析，尝试从原始行提取
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                    except:
                        fixed_lines.append(line)
                        continue

                # 提取元素符号（取原子名称的第一个字母）
                element = atom_name.strip()[0] if atom_name.strip() else 'C'

                # 重新构建 ATOM 行，参考 FSI.pdb 的格式
                # ATOM      1  S1  MOL     1      -1.492  -0.011   0.142  1.00  0.00           S
                # 使用精确的列位置
                fixed_line = (
                    f"ATOM  {atom_num:5d}  {atom_name:3s} {res_name:3s}     {res_num:1d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}  \n"
                )
                fixed_lines.append(fixed_line)
            else:
                fixed_lines.append(line)

        # 写入修复后的文件
        with open(pdb_path, 'w') as f:
            f.writelines(fixed_lines)

        print(f"✅ 修复成功: {pdb_path}")
        return True

    except Exception as e:
        print(f"❌ 修复失败 {pdb_path}: {e}")
        return False


if __name__ == '__main__':
    # 修复 NFBS 和 DFBOP 的 PDB 文件
    pdb_files = [
        Path('data/initial_salts/NFBS.pdb'),
        Path('data/initial_salts/DFBOP.pdb'),
    ]
    
    print("="*80)
    print("修复 PDB 文件格式")
    print("="*80)
    
    for pdb_file in pdb_files:
        if pdb_file.exists():
            fix_pdb_file(pdb_file)
        else:
            print(f"⚠️  文件不存在: {pdb_file}")
    
    print("\n" + "="*80)
    print("修复完成！")
    print("="*80)

