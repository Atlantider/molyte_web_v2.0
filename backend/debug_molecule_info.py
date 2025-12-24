#!/usr/bin/env python3
"""
分子信息服务性能调试脚本
"""
import time
import sys
import os
sys.path.append('/opt/molyte_web_v1.0/backend')

from app.services.molecule_info_service import (
    get_molecule_info_fast,
    get_pubchem_info,
    COMMON_MOLECULES
)

def test_common_molecules():
    """测试常见分子查询速度"""
    print("=== 测试常见分子查询速度 ===")
    
    test_smiles = ["CCO", "CC(=O)O", "C1COC(=O)O1"]
    
    for smiles in test_smiles:
        print(f"\n查询分子: {smiles}")
        
        # 检查是否在常见分子列表中
        if smiles in COMMON_MOLECULES:
            print(f"  在常见分子列表中: {COMMON_MOLECULES[smiles]}")
        else:
            print("  不在常见分子列表中")
        
        start_time = time.time()
        mol_info = get_molecule_info_fast(smiles)
        end_time = time.time()
        
        print(f"  查询耗时: {end_time - start_time:.4f} 秒")
        print(f"  结果: {mol_info}")

def test_pubchem_api():
    """测试PubChem API速度"""
    print("\n=== 测试PubChem API速度 ===")
    
    # 测试一个不在常见分子列表中的分子
    test_smiles = "CC(C)O"  # 异丙醇
    
    print(f"查询分子: {test_smiles}")
    
    start_time = time.time()
    pubchem_info = get_pubchem_info(test_smiles, timeout=5)
    end_time = time.time()
    
    print(f"PubChem API 耗时: {end_time - start_time:.4f} 秒")
    print(f"结果: {pubchem_info}")

def test_batch_molecule_info():
    """测试批量分子信息查询"""
    print("\n=== 测试批量分子信息查询 ===")
    
    test_smiles = [
        "CCO",           # 乙醇 (常见)
        "CC(=O)O",       # 乙酸 (常见)
        "C1COC(=O)O1",   # 碳酸乙烯酯 (常见)
        "CC(C)O",        # 异丙醇 (不常见)
        "CCCCCCCC"       # 辛烷 (不常见)
    ]
    
    total_start = time.time()
    
    for i, smiles in enumerate(test_smiles):
        print(f"\n分子 {i+1}: {smiles}")
        
        start_time = time.time()
        mol_info = get_molecule_info_fast(smiles)
        end_time = time.time()
        
        print(f"  耗时: {end_time - start_time:.4f} 秒")
        print(f"  名称: {mol_info.get('name', 'Unknown')}")
        print(f"  CAS: {mol_info.get('cas_number', 'N/A')}")
        print(f"  分子式: {mol_info.get('molecular_formula', 'N/A')}")
    
    total_end = time.time()
    print(f"\n总耗时: {total_end - total_start:.4f} 秒")

def test_timeout_behavior():
    """测试超时行为"""
    print("\n=== 测试超时行为 ===")
    
    # 使用一个可能导致超时的SMILES
    test_smiles = "CCCCCCCCCCCCCCCCCCCC"  # 长链烷烃
    
    print(f"查询分子: {test_smiles}")
    
    # 测试短超时
    start_time = time.time()
    pubchem_info = get_pubchem_info(test_smiles, timeout=1)  # 1秒超时
    end_time = time.time()
    
    print(f"短超时(1秒) 耗时: {end_time - start_time:.4f} 秒")
    print(f"结果: {pubchem_info}")
    
    # 测试快速查询
    start_time = time.time()
    mol_info = get_molecule_info_fast(test_smiles)
    end_time = time.time()
    
    print(f"快速查询耗时: {end_time - start_time:.4f} 秒")
    print(f"结果: {mol_info}")

def main():
    print("开始分子信息服务性能诊断...")
    
    test_common_molecules()
    test_pubchem_api()
    test_batch_molecule_info()
    test_timeout_behavior()
    
    print("\n诊断完成!")

if __name__ == "__main__":
    main()
