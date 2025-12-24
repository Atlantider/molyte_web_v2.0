"""
Test predict endpoint
"""

import sys
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from app.services.predict_block_service import predict_property, find_similar_molecules


def test_predict_endpoint():
    """测试预测端点"""
    print("\n" + "="*60)
    print("Testing Predict Endpoint")
    print("="*60)
    
    test_smiles = "CCO"  # 乙醇
    
    print(f"\nTesting prediction for SMILES: {test_smiles}")
    
    # 测试 BP 预测
    print("\n1. Testing BP prediction...")
    bp = predict_property(test_smiles, "bp")
    if bp is not None:
        print(f"✓ BP prediction: {bp:.2f}")
    else:
        print("✗ BP prediction failed")
    
    # 测试 MP 预测
    print("\n2. Testing MP prediction...")
    mp = predict_property(test_smiles, "mp")
    if mp is not None:
        print(f"✓ MP prediction: {mp:.2f}")
    else:
        print("✗ MP prediction failed")
    
    # 测试 FP 预测
    print("\n3. Testing FP prediction...")
    fp = predict_property(test_smiles, "fp")
    if fp is not None:
        print(f"✓ FP prediction: {fp:.2f}")
    else:
        print("✗ FP prediction failed")
    
    # 测试相似分子查找
    print("\n4. Testing similar molecules search...")
    similar = find_similar_molecules(test_smiles, n_similar=3)
    if similar:
        print(f"✓ Found {len(similar)} similar molecules")
        for i, mol in enumerate(similar, 1):
            print(f"  {i}. SMILES: {mol['smiles']}, BP: {mol['bp']:.2f}, Similarity: {mol['similarity']:.4f}")
    else:
        print("✗ Similar molecules search failed")
    
    print("\n" + "="*60)
    print("✓ All tests completed!")
    print("="*60)


if __name__ == "__main__":
    try:
        test_predict_endpoint()
    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

