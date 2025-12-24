#!/usr/bin/env python3
"""
为完整数据集预计算分子指纹
"""
import pandas as pd
import numpy as np
import pickle
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
import logging
import time

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

DATA_PATH = Path("/opt/molyte_web_v1.0/backend/data")

def calculate_molecular_fingerprint(smiles: str) -> np.ndarray:
    """计算分子指纹"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # 使用Morgan指纹
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        return np.array(fp)
    except Exception as e:
        logger.error(f"Failed to calculate fingerprint for {smiles}: {e}")
        return None

def precompute_fingerprints_for_dataset(csv_file: Path, output_file: Path, batch_size: int = 5000):
    """为数据集预计算指纹（分批处理以节省内存）"""
    logger.info(f"开始预计算指纹: {csv_file.name}")

    df = pd.read_csv(csv_file)
    logger.info(f"加载数据: {len(df)} 行")

    all_fingerprints = []
    all_valid_indices = []

    start_time = time.time()

    # 分批处理以节省内存
    num_batches = (len(df) + batch_size - 1) // batch_size
    logger.info(f"分 {num_batches} 批处理，每批 {batch_size} 行")

    for batch_idx in range(num_batches):
        batch_start = batch_idx * batch_size
        batch_end = min((batch_idx + 1) * batch_size, len(df))
        batch_df = df.iloc[batch_start:batch_end]

        batch_fingerprints = []
        batch_valid_indices = []

        for idx, row in batch_df.iterrows():
            smiles = row.get('SMILES')
            if pd.notna(smiles):
                fp = calculate_molecular_fingerprint(smiles)
                if fp is not None:
                    batch_fingerprints.append(fp)
                    batch_valid_indices.append(idx)

        all_fingerprints.extend(batch_fingerprints)
        all_valid_indices.extend(batch_valid_indices)

        elapsed = time.time() - start_time
        rate = (batch_end) / elapsed if elapsed > 0 else 0
        logger.info(f"进度: {batch_end}/{len(df)} ({batch_end/len(df)*100:.1f}%) - 速度: {rate:.0f} 分子/秒 - 有效指纹: {len(all_fingerprints)}")

    elapsed = time.time() - start_time
    logger.info(f"指纹计算完成: {len(all_fingerprints)} 个有效指纹，耗时 {elapsed:.1f} 秒")

    # 保存指纹
    logger.info(f"开始保存指纹到: {output_file}")
    fingerprint_data = {
        'fingerprints': np.array(all_fingerprints),
        'valid_indices': all_valid_indices,
        'total_samples': len(df)
    }

    try:
        with open(output_file, 'wb') as f:
            pickle.dump(fingerprint_data, f)

        file_size = output_file.stat().st_size / 1024 / 1024
        logger.info(f"✓ 指纹已保存到: {output_file}")
        logger.info(f"✓ 文件大小: {file_size:.1f} MB")
        logger.info(f"✓ 指纹数量: {len(all_fingerprints)}")
        logger.info(f"✓ 有效索引数量: {len(all_valid_indices)}")
    except Exception as e:
        logger.error(f"✗ 保存指纹失败: {e}")
        raise

def main():
    """主函数"""
    import sys

    # 如果提供了参数，只处理指定的数据集
    if len(sys.argv) > 1:
        csv_name = sys.argv[1]
        pkl_name = csv_name.replace('.csv', '_fingerprints.pkl')
    else:
        # 默认处理 BP 数据集
        csv_name = 'bp_data.csv'
        pkl_name = 'bp_fingerprints.pkl'

    csv_file = DATA_PATH / csv_name
    pkl_file = DATA_PATH / pkl_name

    if csv_file.exists():
        logger.info(f"\n{'='*60}")
        logger.info(f"处理: {csv_name}")
        logger.info(f"{'='*60}")

        try:
            precompute_fingerprints_for_dataset(csv_file, pkl_file)
            logger.info(f"\n✓ {csv_name} 指纹预计算完成！")
        except Exception as e:
            logger.error(f"✗ 处理 {csv_name} 时出错: {e}")
            import traceback
            traceback.print_exc()
    else:
        logger.warning(f"文件不存在: {csv_file}")

if __name__ == "__main__":
    main()
