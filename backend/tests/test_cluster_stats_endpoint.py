"""
Test cluster-stats endpoint
"""

import sys
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from app.services.predict_block_service import get_cluster_statistics


def test_cluster_stats_endpoint():
    """测试聚类统计端点"""
    print("\n" + "="*60)
    print("Testing Cluster Stats Endpoint")
    print("="*60)
    
    # 获取聚类统计信息
    print("\n获取聚类统计信息...")
    stats = get_cluster_statistics()
    
    if stats is None:
        print("✗ 无法获取聚类统计信息")
        return False
    
    print("✓ 成功获取聚类统计信息")
    
    # 验证返回的字段
    required_fields = [
        'total_clusters',
        'total_molecules',
        'max_cluster_size',
        'min_cluster_size',
        'avg_cluster_size',
        'cluster_sizes'
    ]
    
    print("\n验证返回字段...")
    for field in required_fields:
        if field not in stats:
            print(f"✗ 缺少字段: {field}")
            return False
        print(f"✓ {field}: {stats[field]}")
    
    # 验证数据类型
    print("\n验证数据类型...")
    assert isinstance(stats['total_clusters'], int), "total_clusters 应该是 int"
    assert isinstance(stats['total_molecules'], int), "total_molecules 应该是 int"
    assert isinstance(stats['max_cluster_size'], int), "max_cluster_size 应该是 int"
    assert isinstance(stats['min_cluster_size'], int), "min_cluster_size 应该是 int"
    assert isinstance(stats['avg_cluster_size'], float), "avg_cluster_size 应该是 float"
    assert isinstance(stats['cluster_sizes'], dict), "cluster_sizes 应该是 dict"
    print("✓ 所有字段类型正确")
    
    # 验证数据合理性
    print("\n验证数据合理性...")
    assert stats['total_clusters'] > 0, "total_clusters 应该大于 0"
    assert stats['total_molecules'] > 0, "total_molecules 应该大于 0"
    assert stats['max_cluster_size'] > 0, "max_cluster_size 应该大于 0"
    assert stats['min_cluster_size'] > 0, "min_cluster_size 应该大于 0"
    assert stats['avg_cluster_size'] > 0, "avg_cluster_size 应该大于 0"
    # cluster_sizes 中的分子总数应该等于 total_molecules
    total_molecules_in_clusters = sum(stats['cluster_sizes'].values())
    assert total_molecules_in_clusters == stats['total_molecules'], \
        f"cluster_sizes 中的分子总数 ({total_molecules_in_clusters}) 应该等于 total_molecules ({stats['total_molecules']})"
    print("✓ 所有数据合理")
    
    # 打印统计信息
    print("\n" + "="*60)
    print("聚类统计信息:")
    print("="*60)
    print(f"总聚类数: {stats['total_clusters']}")
    print(f"总分子数: {stats['total_molecules']}")
    print(f"最大聚类大小: {stats['max_cluster_size']}")
    print(f"最小聚类大小: {stats['min_cluster_size']}")
    print(f"平均聚类大小: {stats['avg_cluster_size']:.2f}")
    
    # 显示前10个聚类的大小
    print("\n前10个聚类的大小:")
    cluster_sizes_sorted = sorted(stats['cluster_sizes'].items(), key=lambda x: int(x[1]), reverse=True)
    for cluster_id, size in cluster_sizes_sorted[:10]:
        percentage = (size / stats['total_molecules']) * 100
        print(f"  聚类 #{cluster_id}: {size} 个分子 ({percentage:.2f}%)")
    
    print("\n" + "="*60)
    print("✓ 所有测试通过!")
    print("="*60)
    return True


if __name__ == "__main__":
    try:
        success = test_cluster_stats_endpoint()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\n✗ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

