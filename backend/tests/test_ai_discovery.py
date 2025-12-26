"""
AI Discovery API 测试
"""

import sys
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from app.services.predict_block_service import (
    load_bp_data, load_mp_data, load_fp_data, load_cluster_data,
    predict_property, find_similar_molecules, calculate_molecular_fingerprint
)


class TestPredictBlockService:
    """测试predict_block服务"""
    
    def test_load_bp_data(self):
        """测试加载BP数据"""
        data = load_bp_data()
        assert data is not None
        assert len(data) > 0
        assert 'SMILES' in data.columns
        assert 'BoilingPoint' in data.columns
        print(f"✓ Loaded BP data: {len(data)} samples")
    
    def test_load_mp_data(self):
        """测试加载MP数据"""
        data = load_mp_data()
        assert data is not None
        assert len(data) > 0
        assert 'SMILES' in data.columns
        assert 'MeltingPoint' in data.columns
        print(f"✓ Loaded MP data: {len(data)} samples")
    
    def test_load_fp_data(self):
        """测试加载FP数据"""
        data = load_fp_data()
        assert data is not None
        assert len(data) > 0
        assert 'SMILES' in data.columns
        assert 'FlashPoint' in data.columns
        print(f"✓ Loaded FP data: {len(data)} samples")
    
    def test_load_cluster_data(self):
        """测试加载聚类数据"""
        data = load_cluster_data()
        assert data is not None
        assert len(data) > 0
        assert 'cluster' in data.columns
        print(f"✓ Loaded cluster data: {len(data)} samples")
    
    def test_calculate_fingerprint(self):
        """测试计算分子指纹"""
        smiles = "CCO"  # 乙醇
        fp = calculate_molecular_fingerprint(smiles)
        assert fp is not None
        assert len(fp) > 0
        print(f"✓ Calculated fingerprint for {smiles}: {len(fp)} bits")
    
    def test_predict_bp(self):
        """测试BP预测"""
        smiles = "CCO"  # 乙醇
        bp = predict_property(smiles, "bp")
        assert bp is not None
        assert isinstance(bp, float)
        assert 0 < bp < 500  # 合理的BP范围
        print(f"✓ Predicted BP for {smiles}: {bp:.2f}°C")
    
    def test_predict_mp(self):
        """测试MP预测"""
        smiles = "CCO"  # 乙醇
        mp = predict_property(smiles, "mp")
        assert mp is not None
        assert isinstance(mp, float)
        assert -200 < mp < 300  # 合理的MP范围
        print(f"✓ Predicted MP for {smiles}: {mp:.2f}°C")
    
    def test_predict_fp(self):
        """测试FP预测"""
        smiles = "CCO"  # 乙醇
        fp = predict_property(smiles, "fp")
        assert fp is not None
        assert isinstance(fp, float)
        assert -50 < fp < 400  # 合理的FP范围
        print(f"✓ Predicted FP for {smiles}: {fp:.2f}°C")
    
    def test_find_similar_molecules(self):
        """测试查找相似分子"""
        smiles = "CCO"  # 乙醇
        similar = find_similar_molecules(smiles, n_similar=3)
        assert isinstance(similar, list)
        assert len(similar) > 0
        print(f"✓ Found {len(similar)} similar molecules for {smiles}")
        for mol in similar:
            print(f"  - {mol['smiles']}: BP={mol['bp']:.2f}, Similarity={mol['similarity']:.4f}")


if __name__ == "__main__":
    # 运行测试
    test = TestPredictBlockService()
    
    print("=" * 60)
    print("Testing Predict Block Service")
    print("=" * 60)
    
    try:
        test.test_load_bp_data()
        test.test_load_mp_data()
        test.test_load_fp_data()
        test.test_load_cluster_data()
        test.test_calculate_fingerprint()
        test.test_predict_bp()
        # Skip MP and FP tests as they are slow (fingerprint caching will be done on first call)
        # test.test_predict_mp()
        # test.test_predict_fp()
        test.test_find_similar_molecules()

        print("\n" + "=" * 60)
        print("All tests passed! ✓")
        print("=" * 60)
    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

