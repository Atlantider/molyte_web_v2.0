#!/usr/bin/env python3
"""
测试相似分子API性能
"""
import time
import requests
import json

def test_similar_molecules_api():
    """测试相似分子API"""
    url = "http://localhost:8000/api/v1/ai-discovery/similar-molecules"
    
    params = {
        "smiles": "CCO",
        "n_similar": 10
    }
    
    headers = {
        "Authorization": "Bearer test-token"
    }
    
    print("=== 测试相似分子API ===")
    print(f"查询分子: {params['smiles']}")
    print(f"请求相似分子数量: {params['n_similar']}")
    
    start_time = time.time()
    
    try:
        response = requests.get(url, params=params, headers=headers, timeout=30)
        end_time = time.time()
        
        print(f"API请求耗时: {end_time - start_time:.4f} 秒")
        print(f"HTTP状态码: {response.status_code}")
        
        if response.status_code == 200:
            data = response.json()
            print(f"返回相似分子数量: {data['count']}")
            
            for i, mol in enumerate(data['similar_molecules'][:3]):  # 只显示前3个
                print(f"\n相似分子 {i+1}:")
                print(f"  SMILES: {mol['smiles']}")
                print(f"  相似度: {mol['similarity']:.4f}")
                
                mol_info = mol.get('molecule_info', {})
                print(f"  名称: {mol_info.get('name', 'Unknown')}")
                print(f"  CAS: {mol_info.get('cas_number', 'N/A')}")
                print(f"  分子式: {mol_info.get('molecular_formula', 'N/A')}")
                
                properties = mol.get('properties', {})
                is_real = mol.get('is_real_data', {})
                
                for prop, value in properties.items():
                    real_tag = " (真实)" if is_real.get(prop, False) else " (预测)"
                    print(f"  {prop.upper()}: {value:.2f}{real_tag}")
        else:
            print(f"请求失败: {response.text}")
            
    except requests.exceptions.Timeout:
        print("请求超时")
    except Exception as e:
        print(f"请求错误: {e}")

def test_multiple_requests():
    """测试多次请求的性能"""
    print("\n=== 测试多次请求性能 ===")
    
    test_smiles = ["CCO", "CC(=O)O", "C1COC(=O)O1"]
    url = "http://localhost:8000/api/v1/ai-discovery/similar-molecules"
    headers = {"Authorization": "Bearer test-token"}
    
    total_start = time.time()
    
    for i, smiles in enumerate(test_smiles):
        print(f"\n请求 {i+1}: {smiles}")
        
        start_time = time.time()
        try:
            response = requests.get(
                url, 
                params={"smiles": smiles, "n_similar": 5}, 
                headers=headers, 
                timeout=30
            )
            end_time = time.time()
            
            print(f"  耗时: {end_time - start_time:.4f} 秒")
            print(f"  状态: {response.status_code}")
            
            if response.status_code == 200:
                data = response.json()
                print(f"  找到: {data['count']} 个相似分子")
            
        except Exception as e:
            print(f"  错误: {e}")
    
    total_end = time.time()
    print(f"\n总耗时: {total_end - total_start:.4f} 秒")

if __name__ == "__main__":
    print("开始API性能测试...")
    test_similar_molecules_api()
    test_multiple_requests()
    print("\n测试完成!")
