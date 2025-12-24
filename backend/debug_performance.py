#!/usr/bin/env python3
"""
性能调试脚本
"""
import time
import sys
import os
sys.path.append('/opt/molyte_web_v1.0/backend')

from app.services.predict_block_service import (
    predict_property, 
    calculate_molecular_fingerprint,
    load_precomputed_fingerprints,
    load_bp_data
)

def test_fingerprint_calculation():
    """测试指纹计算速度"""
    print("=== 测试指纹计算速度 ===")
    smiles = "CCO"
    
    start_time = time.time()
    fp = calculate_molecular_fingerprint(smiles)
    end_time = time.time()
    
    print(f"计算 {smiles} 的指纹耗时: {end_time - start_time:.4f} 秒")
    print(f"指纹长度: {len(fp) if fp is not None else 'None'}")
    return fp

def test_precomputed_loading():
    """测试预计算指纹加载速度"""
    print("\n=== 测试预计算指纹加载速度 ===")
    
    start_time = time.time()
    fingerprint_data = load_precomputed_fingerprints('bp')
    end_time = time.time()
    
    print(f"加载预计算指纹耗时: {end_time - start_time:.4f} 秒")
    if fingerprint_data:
        print(f"指纹数量: {len(fingerprint_data['fingerprints'])}")
        print(f"有效索引数量: {len(fingerprint_data['valid_indices'])}")
    else:
        print("加载失败")
    
    return fingerprint_data

def test_data_loading():
    """测试数据加载速度"""
    print("\n=== 测试数据加载速度 ===")
    
    start_time = time.time()
    data = load_bp_data()
    end_time = time.time()
    
    print(f"加载 BP 数据耗时: {end_time - start_time:.4f} 秒")
    if data is not None:
        print(f"数据行数: {len(data)}")
    else:
        print("加载失败")
    
    return data

def test_prediction():
    """测试完整预测速度"""
    print("\n=== 测试完整预测速度 ===")
    smiles = "CCO"
    
    start_time = time.time()
    result = predict_property(smiles, 'bp')
    end_time = time.time()
    
    print(f"预测 {smiles} 的 BP 耗时: {end_time - start_time:.4f} 秒")
    print(f"预测结果: {result}")

def main():
    print("开始性能诊断...")
    
    # 测试各个组件
    test_fingerprint_calculation()
    test_precomputed_loading()
    test_data_loading()
    test_prediction()
    
    print("\n诊断完成!")

if __name__ == "__main__":
    main()
