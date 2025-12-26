#!/usr/bin/env python3
"""
从头重新生成 atom_mapping.json

从 job 数据和 data 文件重新生成正确的 atom_mapping.json
"""

import sys
import json
from pathlib import Path

# 添加项目根目录到 Python 路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.app.workers.atom_mapping_generator import AtomMappingGenerator


def regenerate_from_scratch(work_dir: Path):
    """从头重新生成 atom_mapping.json"""

    # 从 data 文件中读取实际的分子数量
    data_file = work_dir / f"{work_dir.name}.data"
    mol_counts = {}

    with open(data_file, 'r') as f:
        content = f.read()

    # 解析 Atoms 部分
    if 'Atoms' in content:
        atoms_section = content.split('Atoms')[1]
        if 'Bonds' in atoms_section:
            atoms_section = atoms_section.split('Bonds')[0]

        for line in atoms_section.strip().split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    mol_id = int(parts[1])
                    if mol_id not in mol_counts:
                        mol_counts[mol_id] = 0
                    mol_counts[mol_id] += 1
                except ValueError:
                    continue

    max_mol_id = max(mol_counts.keys()) if mol_counts else 0
    print(f"从 data 文件读取: 最大 mol_id={max_mol_id}")

    # 构建 job_data（使用实际的分子数量）
    job_data = {
        "name": work_dir.name,
        "cations": [{"name": "Li", "number": max_mol_id}],
        "anions": [
            {"name": "DFBOP", "number": max_mol_id},
            {"name": "FSI", "number": max_mol_id}
        ],
        "solvents": [
            {"name": "TTE", "number": max_mol_id},
            {"name": "DME", "number": max_mol_id}
        ]
    }
    
    print(f"Job data:")
    print(f"  Cations: {job_data['cations']}")
    print(f"  Anions: {job_data['anions']}")
    print(f"  Solvents: {job_data['solvents']}")
    
    # 创建 AtomMappingGenerator 实例
    generator = AtomMappingGenerator(Path("data/initial_salts"))
    
    # 生成 atom_mapping
    print(f"\n生成 atom_mapping...")
    mapping = generator.generate_atom_mapping(job_data, work_dir)
    
    print(f"\n✓ 生成完成: {len(mapping['molecules'])} 个分子")
    
    # 验证修正结果
    print("\n验证修正结果（FSI 分子）:")
    for mol in mapping['molecules']:
        if mol['molecule_name'] == 'FSI':
            print(f"\nFSI 分子 (mol_id={mol['molecule_id']}):")
            print(f"  原子数: {len(mol['atoms'])}")
            for atom in mol['atoms']:
                print(f"    atom_id={atom['atom_id']:3d}, atom_index={atom['atom_index']:2d}, "
                      f"element={atom['element']:2s}, atom_name={atom['atom_name']:10s}")
            
            elements = [atom['element'] for atom in mol['atoms']]
            print(f"  元素顺序: {' '.join(elements)}")
            print(f"  期望顺序: S S F F O O O O N")
            
            if elements == ['S', 'S', 'F', 'F', 'O', 'O', 'O', 'O', 'N']:
                print("  ✅ FSI 顺序正确")
            else:
                print("  ❌ FSI 顺序错误")
            
            break
    
    return True


if __name__ == "__main__":
    work_dir = Path("data/md_work/EL-20251209-0001-Li-FSI-DFBOP-TTE-DME-EL-20251209-0001-Li-FSI-DFBOP-TTE-DME-MD1")
    
    print(f"从头重新生成 atom_mapping.json for: {work_dir.name}\n")
    
    if regenerate_from_scratch(work_dir):
        print("\n✅ 成功重新生成 atom_mapping.json")
    else:
        print("\n❌ 重新生成失败")
        sys.exit(1)

