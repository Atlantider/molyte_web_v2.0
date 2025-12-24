#!/usr/bin/env python3
"""
聚类功能性能调试脚本
"""
import time
import sys
import os
sys.path.append('/opt/molyte_web_v1.0/backend')

from app.services.predict_block_service import (
    find_similar_molecules,
    calculate_molecular_fingerprint,
    load_precomputed_fingerprints,
    load_bp_data
)

def test_similar_molecules():
    """测试相似分子查找速度"""
    print("=== 测试相似分子查找速度 ===")
    smiles = "CCO"
    
    print(f"查找 {smiles} 的相似分子...")
    
    start_time = time.time()
    similar_molecules = find_similar_molecules(smiles, n_similar=10)
    end_time = time.time()
    
    print(f"查找相似分子耗时: {end_time - start_time:.4f} 秒")
    print(f"找到相似分子数量: {len(similar_molecules)}")
    
    for i, mol in enumerate(similar_molecules[:3]):  # 只显示前3个
        print(f"  {i+1}. {mol['smiles']} (相似度: {mol['similarity']:.4f})")
        print(f"     属性: {mol['properties']}")
    
    return similar_molecules

def test_fingerprint_loading():
    """测试指纹加载"""
    print("\n=== 测试指纹加载 ===")
    
    start_time = time.time()
    fingerprint_data = load_precomputed_fingerprints('bp')
    end_time = time.time()
    
    print(f"加载指纹耗时: {end_time - start_time:.4f} 秒")
    if fingerprint_data:
        print(f"指纹数量: {len(fingerprint_data['fingerprints'])}")
        print(f"有效索引数量: {len(fingerprint_data['valid_indices'])}")
        print(f"指纹形状: {fingerprint_data['fingerprints'].shape}")
    else:
        print("指纹加载失败")
    
    return fingerprint_data

def test_data_loading():
    """测试数据加载"""
    print("\n=== 测试数据加载 ===")
    
    start_time = time.time()
    bp_data = load_bp_data()
    end_time = time.time()
    
    print(f"加载BP数据耗时: {end_time - start_time:.4f} 秒")
    if bp_data is not None:
        print(f"数据行数: {len(bp_data)}")
        print(f"数据列: {list(bp_data.columns)}")
        print(f"前3行SMILES: {bp_data['SMILES'].head(3).tolist()}")
    else:
        print("数据加载失败")
    
    return bp_data

def test_fingerprint_calculation():
    """测试指纹计算"""
    print("\n=== 测试指纹计算 ===")
    
    test_smiles = ["CCO", "CC(=O)O", "C1COC(=O)O1"]
    
    for smiles in test_smiles:
        start_time = time.time()
        fp = calculate_molecular_fingerprint(smiles)
        end_time = time.time()
        
        print(f"{smiles}: {end_time - start_time:.4f} 秒 (长度: {len(fp) if fp is not None else 'None'})")

def test_step_by_step():
    """逐步测试相似分子查找的每个步骤"""
    print("\n=== 逐步测试相似分子查找 ===")
    smiles = "CCO"
    
    # 步骤1：计算查询分子指纹
    print("步骤1: 计算查询分子指纹")
    start_time = time.time()
    query_fp = calculate_molecular_fingerprint(smiles)
    end_time = time.time()
    print(f"  耗时: {end_time - start_time:.4f} 秒")
    
    if query_fp is None:
        print("  查询分子指纹计算失败")
        return
    
    # 步骤2：加载预计算指纹
    print("步骤2: 加载预计算指纹")
    start_time = time.time()
    fingerprint_data = load_precomputed_fingerprints('bp')
    end_time = time.time()
    print(f"  耗时: {end_time - start_time:.4f} 秒")
    
    if not fingerprint_data:
        print("  预计算指纹加载失败")
        return
    
    # 步骤3：加载数据
    print("步骤3: 加载BP数据")
    start_time = time.time()
    bp_data = load_bp_data()
    end_time = time.time()
    print(f"  耗时: {end_time - start_time:.4f} 秒")
    
    if bp_data is None:
        print("  BP数据加载失败")
        return
    
    # 步骤4：KNN搜索
    print("步骤4: KNN相似性搜索")
    from sklearn.neighbors import NearestNeighbors
    import numpy as np
    
    fingerprints = fingerprint_data['fingerprints']
    valid_indices = fingerprint_data['valid_indices']
    
    start_time = time.time()
    knn = NearestNeighbors(n_neighbors=min(10, len(fingerprints)), metric='euclidean')
    knn.fit(fingerprints)
    distances, indices = knn.kneighbors([query_fp])
    end_time = time.time()
    print(f"  耗时: {end_time - start_time:.4f} 秒")
    print(f"  找到 {len(indices[0])} 个邻居")
    
    # 步骤5：构建结果
    print("步骤5: 构建结果")
    start_time = time.time()
    results = []
    for i, idx in enumerate([valid_indices[j] for j in indices[0]]):
        if idx >= len(bp_data):
            continue
        row = bp_data.iloc[idx]
        results.append({
            'smiles': row['SMILES'],
            'similarity': float(1.0 / (distances[0][i] + 1e-6)),
            'properties': {'bp': float(row['BoilingPoint'])}
        })
    end_time = time.time()
    print(f"  耗时: {end_time - start_time:.4f} 秒")
    print(f"  构建了 {len(results)} 个结果")

def main():
    print("开始聚类功能性能诊断...")
    
    # 测试各个组件
    test_fingerprint_calculation()
    test_fingerprint_loading()
    test_data_loading()
    test_step_by_step()
    test_similar_molecules()
    
    print("\n诊断完成!")

if __name__ == "__main__":
    main()
