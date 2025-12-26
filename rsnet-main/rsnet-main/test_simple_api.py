"""
RSNet Simple API - 综合测试套件

测试新的简化API的各项功能
"""

import pytest
import tempfile
import shutil
from pathlib import Path

from rsnet_simple_api import RSNetSimpleAPI, generate_reaction_network


class TestRSNetSimpleAPI:
    """测试RSNet Simple API"""
    
    @pytest.fixture
    def temp_output_dir(self):
        """创建临时输出目录"""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)
    
    @pytest.fixture
    def api(self, temp_output_dir):
        """创建API实例"""
        return RSNetSimpleAPI(output_dir=temp_output_dir)
    
    def test_basic_network_generation(self, api):
        """测试基本网络生成"""
        result = api.generate_reaction_network(
            smiles_list=['CCO'],  # 乙醇
            temperature=300.0,
            max_generations=1,
            max_species=10,
            visualize=False,
            save_results=False
        )
        
        assert 'network' in result
        assert 'molecules' in result
        assert 'reactions' in result
        assert 'statistics' in result
        assert len(result['molecules']) >= 1
    
    def test_multiple_molecules(self, api):
        """测试多个分子输入"""
        result = api.generate_reaction_network(
            smiles_list=['C', 'C=C', 'CCO'],
            temperature=400.0,
            max_generations=2,
            max_species=20,
            visualize=False,
            save_results=False
        )
        
        assert len(result['molecules']) >= 3
        assert result['statistics']['num_molecules'] >= 3
    
    def test_visualization_generation(self, api, temp_output_dir):
        """测试可视化生成"""
        result = api.generate_reaction_network(
            smiles_list=['CCO'],
            temperature=300.0,
            max_generations=1,
            visualize=True,
            save_results=False
        )
        
        assert 'visualization_path' in result
        viz_path = Path(result['visualization_path'])
        assert viz_path.exists()
        assert viz_path.suffix == '.png'
    
    def test_results_saving(self, api, temp_output_dir):
        """测试结果保存"""
        result = api.generate_reaction_network(
            smiles_list=['C'],
            temperature=300.0,
            max_generations=1,
            visualize=False,
            save_results=True
        )
        
        assert 'json_path' in result
        json_path = Path(result['json_path'])
        assert json_path.exists()
        assert json_path.suffix == '.json'
        
        # 验证JSON内容
        import json
        with open(json_path, 'r') as f:
            data = json.load(f)
        
        assert 'statistics' in data
        assert 'molecules' in data
        assert 'reactions' in data
    
    def test_different_temperatures(self, api):
        """测试不同温度"""
        temps = [300.0, 500.0, 800.0]
        
        for temp in temps:
            result = api.generate_reaction_network(
                smiles_list=['C'],
                temperature=temp,
                max_generations=1,
                max_species=10,
                visualize=False,
                save_results=False
            )
            
            assert result['environment']['temperature'] == temp
            assert 'network' in result
    
    def test_electrode_types(self, api):
        """测试不同电极类型"""
        for electrode in ['anode', 'cathode']:
            result = api.generate_reaction_network(
                smiles_list=['CCO'],
                temperature=300.0,
                electrode_type=electrode,
                max_generations=1,
                visualize=False,
                save_results=False
            )
            
            assert result['environment']['electrode_type'] == electrode
    
    def test_generation_limits(self, api):
        """测试代数限制"""
        result = api.generate_reaction_network(
            smiles_list=['CCO'],
            temperature=300.0,
            max_generations=2,
            max_species=15,
            visualize=False,
            save_results=False
        )
        
        stats = result['statistics']
        assert stats['max_generation'] <= 2
        assert stats['num_molecules'] <= 15
    
    def test_invalid_smiles(self, api):
        """测试无效SMILES处理"""
        # 应该能处理部分无效的SMILES
        result = api.generate_reaction_network(
            smiles_list=['CCO', 'INVALID_SMILES', 'C'],
            temperature=300.0,
            max_generations=1,
            visualize=False,
            save_results=False
        )
        
        # 应该至少有有效的分子
        assert len(result['molecules']) >= 2
    
    def test_convenience_function(self, temp_output_dir):
        """测试便捷函数"""
        result = generate_reaction_network(
            smiles_list=['C'],
            temperature=300.0,
            max_generations=1,
            output_dir=temp_output_dir,
            visualize=False,
            save_results=False
        )
        
        assert 'network' in result
        assert 'molecules' in result
        assert 'reactions' in result


class TestNetworkProperties:
    """测试网络属性"""
    
    def test_reaction_energy(self):
        """测试反应能量计算"""
        result = generate_reaction_network(
            smiles_list=['CCO'],
            temperature=300.0,
            max_generations=1,
            visualize=False,
            save_results=False
        )
        
        # 检查反应是否有能量信息
        for rxn in result['reactions']:
            # 能量可能存在或不存在，取决于计算
            if hasattr(rxn, 'reaction_energy'):
                assert isinstance(rxn.reaction_energy, (int, float))
    
    def test_molecule_properties(self):
        """测试分子属性"""
        result = generate_reaction_network(
            smiles_list=['CCO', 'C=C'],
            temperature=300.0,
            max_generations=1,
            visualize=False,
            save_results=False
        )
        
        for mol in result['molecules']:
            assert hasattr(mol, 'smiles')
            assert hasattr(mol, 'name')
            assert mol.smiles is not None
            assert mol.name is not None
    
    def test_statistics_completeness(self):
        """测试统计信息完整性"""
        result = generate_reaction_network(
            smiles_list=['C'],
            temperature=300.0,
            max_generations=2,
            visualize=False,
            save_results=False
        )
        
        stats = result['statistics']
        required_keys = ['num_molecules', 'num_reactions', 'max_generation']
        
        for key in required_keys:
            assert key in stats
            assert isinstance(stats[key], int)
            assert stats[key] >= 0


class TestPerformance:
    """测试性能"""
    
    def test_generation_time(self):
        """测试生成时间记录"""
        result = generate_reaction_network(
            smiles_list=['C'],
            temperature=300.0,
            max_generations=1,
            visualize=False,
            save_results=False
        )
        
        assert 'generation_time' in result
        assert result['generation_time'] > 0
        assert result['generation_time'] < 300  # 应该在5分钟内完成
    
    def test_small_network_fast(self):
        """测试小网络快速生成"""
        import time
        
        start = time.time()
        result = generate_reaction_network(
            smiles_list=['C'],
            temperature=300.0,
            max_generations=1,
            max_species=5,
            visualize=False,
            save_results=False
        )
        elapsed = time.time() - start
        
        # 小网络应该很快
        assert elapsed < 60  # 1分钟内


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
