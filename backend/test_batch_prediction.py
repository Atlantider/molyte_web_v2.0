#!/usr/bin/env python3
"""
测试批量预测性能
"""
import time
import requests
import json

def test_batch_prediction():
    """测试批量预测API"""
    url = "http://localhost:8000/api/v1/ai-discovery/predict-batch"
    
    # 测试数据
    test_smiles = [
        "CCO",           # 乙醇
        "CC(=O)O",       # 乙酸
        "C1COC(=O)O1"    # 碳酸乙烯酯
    ]
    
    payload = {
        "smiles_list": test_smiles,
        "properties": ["bp", "mp", "fp"],
        "include_qm9": True
    }
    
    headers = {
        "Content-Type": "application/json",
        "Authorization": "Bearer test-token"  # 需要有效的token
    }
    
    print("=== 测试批量预测API ===")
    print(f"测试分子数量: {len(test_smiles)}")
    print(f"SMILES: {test_smiles}")
    
    start_time = time.time()
    
    try:
        response = requests.post(url, json=payload, headers=headers, timeout=30)
        end_time = time.time()
        
        print(f"请求耗时: {end_time - start_time:.4f} 秒")
        print(f"HTTP状态码: {response.status_code}")
        
        if response.status_code == 200:
            data = response.json()
            print(f"成功预测: {data['success_count']}/{data['total_count']}")
            print(f"失败分子: {data['failed_smiles']}")
            
            for i, result in enumerate(data['results']):
                print(f"\n分子 {i+1}: {result['smiles']}")
                mol_info = result.get('molecule_info', {})
                print(f"  名称: {mol_info.get('name', 'Unknown')}")
                print(f"  CAS: {mol_info.get('cas_number', 'N/A')}")
                print(f"  分子式: {mol_info.get('molecular_formula', 'N/A')}")
                
                properties = result.get('predicted_properties', {})
                is_real = result.get('is_real_data', {})
                
                for prop, value in properties.items():
                    real_tag = " (真实)" if is_real.get(prop, False) else " (预测)"
                    print(f"  {prop.upper()}: {value:.2f}{real_tag}")
        else:
            print(f"请求失败: {response.text}")
            
    except requests.exceptions.Timeout:
        print("请求超时")
    except Exception as e:
        print(f"请求错误: {e}")

def test_single_prediction():
    """测试单个预测API（对比）"""
    url = "http://localhost:8000/api/v1/ai-discovery/predict"
    
    test_smiles = [
        "CCO",           # 乙醇
        "CC(=O)O",       # 乙酸
        "C1COC(=O)O1"    # 碳酸乙烯酯
    ]
    
    headers = {
        "Content-Type": "application/json",
        "Authorization": "Bearer test-token"
    }
    
    print("\n=== 测试单个预测API（串行）===")
    print(f"测试分子数量: {len(test_smiles)}")
    
    start_time = time.time()
    success_count = 0
    
    for smiles in test_smiles:
        payload = {
            "smiles": smiles,
            "properties": ["bp", "mp", "fp"],
            "include_qm9": True
        }
        
        try:
            response = requests.post(url, json=payload, headers=headers, timeout=30)
            if response.status_code == 200:
                success_count += 1
                print(f"  {smiles}: 成功")
            else:
                print(f"  {smiles}: 失败 - {response.status_code}")
        except Exception as e:
            print(f"  {smiles}: 错误 - {e}")
    
    end_time = time.time()
    print(f"串行预测耗时: {end_time - start_time:.4f} 秒")
    print(f"成功预测: {success_count}/{len(test_smiles)}")

if __name__ == "__main__":
    print("开始性能测试...")
    test_batch_prediction()
    test_single_prediction()
    print("\n测试完成!")
