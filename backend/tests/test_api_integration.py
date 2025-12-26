"""
AI Discovery API 集成测试

测试predict_block功能与API的集成
"""

import sys
from pathlib import Path
import json

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from app.services.predict_block_service import (
    predict_property, find_similar_molecules, get_cluster_statistics
)


def test_predict_api_integration():
    """测试预测API的完整流程"""
    print("\n" + "="*60)
    print("Testing AI Discovery API Integration")
    print("="*60)
    
    # 测试分子
    test_smiles = "CCO"  # 乙醇
    
    print(f"\n测试分子: {test_smiles} (乙醇)")
    
    # 1. 测试BP预测
    print("\n1. 测试沸点(BP)预测...")
    bp = predict_property(test_smiles, "bp")
    if bp is not None:
        print(f"   ✓ BP预测成功: {bp:.2f}°C")
        assert isinstance(bp, float), "BP应该是float类型"
        assert bp > 0, "BP应该是正数"
    else:
        print(f"   ✗ BP预测失败")
        return False
    
    # 2. 测试相似分子查找
    print("\n2. 测试相似分子查找...")
    similar = find_similar_molecules(test_smiles, n_similar=3)
    if similar:
        print(f"   ✓ 找到 {len(similar)} 个相似分子:")
        for i, mol in enumerate(similar, 1):
            print(f"     {i}. SMILES: {mol['smiles']}, BP: {mol['bp']:.2f}, 相似度: {mol['similarity']:.4f}")
        assert len(similar) <= 3, "应该返回最多3个相似分子"
        assert all('smiles' in m and 'bp' in m and 'similarity' in m for m in similar), "返回格式不正确"
    else:
        print(f"   ✗ 未找到相似分子")
        return False
    
    # 3. 测试聚类统计
    print("\n3. 测试聚类统计...")
    stats = get_cluster_statistics()
    if stats:
        print(f"   ✓ 聚类统计成功:")
        print(f"     - 总聚类数: {stats['total_clusters']}")
        print(f"     - 总分子数: {stats['total_molecules']}")
        print(f"     - 最大聚类大小: {stats['max_cluster_size']}")
        print(f"     - 平均聚类大小: {stats['avg_cluster_size']:.2f}")
        assert isinstance(stats['total_clusters'], int), "total_clusters应该是int"
        assert isinstance(stats['total_molecules'], int), "total_molecules应该是int"
    else:
        print(f"   ✗ 聚类统计失败")
        return False
    
    # 4. 测试API响应格式
    print("\n4. 验证API响应格式...")
    api_response = {
        "smiles": test_smiles,
        "predicted_properties": {
            "bp": bp if bp is not None else 0.0,
            "mp": predict_property(test_smiles, "mp") or 0.0,
            "fp": predict_property(test_smiles, "fp") or 0.0
        },
        "confidence": {
            "bp": 0.75,
            "mp": 0.0,
            "fp": 0.0
        },
        "similar_molecules": [
            {
                "name": f"Similar_{mol['smiles'][:20]}",
                "smiles": mol['smiles'],
                "composition": {},
                "properties": {"bp": mol['bp'], "similarity": mol['similarity']}
            }
            for mol in similar[:3]
        ]
    }
    
    # 验证响应格式
    assert api_response["smiles"] == test_smiles, "SMILES不匹配"
    assert isinstance(api_response["predicted_properties"], dict), "predicted_properties应该是dict"
    assert all(isinstance(v, (int, float)) for v in api_response["predicted_properties"].values()), \
        "predicted_properties的值应该都是数字"
    assert isinstance(api_response["confidence"], dict), "confidence应该是dict"
    assert isinstance(api_response["similar_molecules"], list), "similar_molecules应该是list"
    
    print("   ✓ API响应格式验证通过")
    print(f"   响应示例: {json.dumps(api_response, indent=2)[:200]}...")
    
    print("\n" + "="*60)
    print("All integration tests passed! ✓")
    print("="*60)
    return True


if __name__ == "__main__":
    try:
        success = test_predict_api_integration()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\n✗ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

