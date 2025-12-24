#!/usr/bin/env python3
"""
预计算分子指纹脚本

这个脚本会预先计算所有分子的指纹并保存到文件中，
这样后续的预测就不需要重新计算指纹了。
"""

import pandas as pd
import numpy as np
import pickle
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# 数据文件路径
BACKEND_DATA_PATH = Path("/opt/molyte_web_v1.0/backend/data")

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

def precompute_fingerprints_for_dataset(csv_file: Path, output_file: Path):
    """为数据集预计算指纹"""
    if not csv_file.exists():
        logger.warning(f"Dataset file not found: {csv_file}")
        return
    
    logger.info(f"Loading dataset: {csv_file}")
    df = pd.read_csv(csv_file)
    
    fingerprints = []
    valid_indices = []
    
    logger.info(f"Computing fingerprints for {len(df)} molecules...")
    for idx, row in df.iterrows():
        if idx % 100 == 0:
            logger.info(f"Progress: {idx}/{len(df)} ({idx/len(df)*100:.1f}%)")
        
        smiles = row.get('SMILES')
        if pd.notna(smiles):
            fp = calculate_molecular_fingerprint(smiles)
            if fp is not None:
                fingerprints.append(fp)
                valid_indices.append(idx)
    
    # 保存指纹和有效索引
    fingerprint_data = {
        'fingerprints': np.array(fingerprints),
        'valid_indices': valid_indices,
        'total_samples': len(df)
    }
    
    with open(output_file, 'wb') as f:
        pickle.dump(fingerprint_data, f)
    
    logger.info(f"Saved {len(fingerprints)} fingerprints to {output_file}")

def main():
    """主函数"""
    logger.info("Starting fingerprint precomputation...")
    
    # 为小型数据集预计算指纹
    datasets = [
        ('bp_data_small.csv', 'bp_fingerprints_small.pkl'),
        ('mp_data_small.csv', 'mp_fingerprints_small.pkl'),
        ('fp_data_small.csv', 'fp_fingerprints_small.pkl'),
        ('cluster_data_small.csv', 'cluster_fingerprints_small.pkl')
    ]
    
    for csv_name, pkl_name in datasets:
        csv_file = BACKEND_DATA_PATH / csv_name
        output_file = BACKEND_DATA_PATH / pkl_name
        
        if output_file.exists():
            logger.info(f"Fingerprints already exist: {output_file}")
            continue
        
        precompute_fingerprints_for_dataset(csv_file, output_file)
    
    logger.info("Fingerprint precomputation completed!")

if __name__ == "__main__":
    main()
