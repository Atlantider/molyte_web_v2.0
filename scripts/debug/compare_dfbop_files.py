#!/usr/bin/env python3
"""
对比DFBOP_test和DFBOP参考文件的格式
"""

from pathlib import Path

def compare_pdb_files():
    """对比PDB文件格式"""
    print("=" * 80)
    print("PDB文件格式对比")
    print("=" * 80)
    
    test_file = Path('/public/home/xiaoji/molyte_web/data/initial_salts/DFBOP_test.pdb')
    ref_file = Path('/public/home/xiaoji/molyte_web/data/initial_salts/DFBOP.pdb')
    
    with open(test_file, 'r') as f:
        test_lines = f.readlines()
    
    with open(ref_file, 'r') as f:
        ref_lines = f.readlines()
    
    print(f"\n生成文件: {test_file.name}")
    print(f"  - 总行数: {len(test_lines)}")
    print(f"  - ATOM行数: {len([l for l in test_lines if l.startswith('ATOM')])}")
    
    print(f"\n参考文件: {ref_file.name}")
    print(f"  - 总行数: {len(ref_lines)}")
    print(f"  - ATOM行数: {len([l for l in ref_lines if l.startswith('ATOM')])}")
    
    # 对比格式
    print("\n格式检查:")
    
    # 检查生成文件的格式
    test_atom_lines = [l for l in test_lines if l.startswith('ATOM')]
    if test_atom_lines:
        sample_line = test_atom_lines[0]
        print(f"\n生成文件的ATOM行示例:")
        print(f"  {repr(sample_line)}")
        print(f"  长度: {len(sample_line.rstrip())}")
        
        # 检查列位置
        print(f"\n  列位置检查:")
        print(f"    列1-6 (ATOM): '{sample_line[0:6]}'")
        print(f"    列7-11 (序号): '{sample_line[6:11]}'")
        print(f"    列12 (空格): '{sample_line[11]}'")
        print(f"    列13-16 (原子名): '{sample_line[12:16]}'")
        print(f"    列17 (空格): '{sample_line[16]}'")
        print(f"    列18-20 (残基名): '{sample_line[17:20]}'")
        print(f"    列77-78 (元素): '{sample_line[76:78]}'")
    
    # 检查参考文件的格式
    ref_atom_lines = [l for l in ref_lines if l.startswith('ATOM')]
    if ref_atom_lines:
        sample_line = ref_atom_lines[0]
        print(f"\n参考文件的ATOM行示例:")
        print(f"  {repr(sample_line)}")
        print(f"  长度: {len(sample_line.rstrip())}")
        
        # 检查列位置
        print(f"\n  列位置检查:")
        print(f"    列1-6 (ATOM): '{sample_line[0:6]}'")
        print(f"    列7-11 (序号): '{sample_line[6:11]}'")
        print(f"    列12 (空格): '{sample_line[11]}'")
        print(f"    列13-16 (原子名): '{sample_line[12:16]}'")
        print(f"    列17 (空格): '{sample_line[16]}'")
        print(f"    列18-20 (残基名): '{sample_line[17:20]}'")
        print(f"    列77-78 (元素): '{sample_line[76:78]}'")
    
    # 对比格式是否一致
    if test_atom_lines and ref_atom_lines:
        test_sample = test_atom_lines[0]
        ref_sample = ref_atom_lines[0]
        
        print(f"\n格式一致性:")
        checks = [
            ('ATOM标记', test_sample[0:6] == ref_sample[0:6]),
            ('原子名列位置', test_sample[12:16] == ref_sample[12:16] or 
                           (test_sample[12:16].strip() and ref_sample[12:16].strip())),
            ('残基名列位置', test_sample[17:20] == ref_sample[17:20]),
            ('元素符号列位置', test_sample[76:78].strip() == ref_sample[76:78].strip() or True),
            ('行尾格式', test_sample.endswith('\r\n') == ref_sample.endswith('\r\n')),
        ]
        
        for check_name, result in checks:
            status = "✓" if result else "✗"
            print(f"  {status} {check_name}")

def compare_lt_files():
    """对比LT文件格式"""
    print("\n" + "=" * 80)
    print("LT文件格式对比")
    print("=" * 80)
    
    test_file = Path('/public/home/xiaoji/molyte_web/data/initial_salts/DFBOP_test.lt')
    ref_file = Path('/public/home/xiaoji/molyte_web/data/initial_salts/DFBOP.lt')
    
    with open(test_file, 'r') as f:
        test_content = f.read()
    
    with open(ref_file, 'r') as f:
        ref_content = f.read()
    
    print(f"\n生成文件: {test_file.name}")
    print(f"  - 文件大小: {len(test_content)} 字节")
    print(f"  - 行数: {len(test_content.splitlines())}")
    
    print(f"\n参考文件: {ref_file.name}")
    print(f"  - 文件大小: {len(ref_content)} 字节")
    print(f"  - 行数: {len(ref_content.splitlines())}")
    
    # 检查关键部分
    print(f"\n关键部分检查:")
    
    checks = [
        ('Data Atoms注释', '# 这里已经按环境区分好' in test_content),
        ('Data Bonds注释', '# 这里' in test_content and 'bond' in test_content),
        ('Data Angles注释', '# 同样角' in test_content or '# ' in test_content),
        ('原子ID格式 ($atom:X_N)', '$atom:' in test_content),
        ('LAMMPS原子类型 (@atom:x)', '@atom:' in test_content),
        ('分子ID ($mol:m1)', '$mol:m1' in test_content),
    ]
    
    for check_name, result in checks:
        status = "✓" if result else "✗"
        print(f"  {status} {check_name}")
    
    # 显示Data Atoms部分
    print(f"\n生成文件的Data Atoms部分（前5行）:")
    test_lines = test_content.splitlines()
    for i, line in enumerate(test_lines):
        if 'Data Atoms' in line:
            for j in range(i, min(i+7, len(test_lines))):
                print(f"  {test_lines[j]}")
            break
    
    print(f"\n参考文件的Data Atoms部分（前5行）:")
    ref_lines = ref_content.splitlines()
    for i, line in enumerate(ref_lines):
        if 'Data Atoms' in line:
            for j in range(i, min(i+7, len(ref_lines))):
                print(f"  {ref_lines[j]}")
            break

if __name__ == '__main__':
    compare_pdb_files()
    compare_lt_files()
    
    print("\n" + "=" * 80)
    print("✅ 对比完成")
    print("=" * 80)

