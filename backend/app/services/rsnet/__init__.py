"""
RSNet Service Module
反应网络生成服务 - 统一接口

this module provides a clean interface to the RSNet reaction network generation system
"""

from .core.molecule import Molecule
from .core.environment import Environment  
from .core.reaction import Reaction

__all__ = ['Molecule', 'Environment', 'Reaction', 'RSNetService']

class RSNetService:
    """
    RSNet服务统一接口
    
    提供简化的API调用反应网络生成功能
    """
    
    def __init__(self):
        """初始化RSNet服务"""
        # Lazy import to avoid circular dependencies
        pass
    
    def generate_network(self, **kwargs):
        """
        生成反应网络
        
        Args:
            smiles_list: 初始分子SMILES列表
            temperature: 温度(K)
            electrode_type: 电极类型
            voltage: 电压(V)
            max_generations: 最大代数
            max_species: 最大分子数
            energy_cutoff: 能量截断值
            
        Returns:
            Dict containing network, molecules, reactions, statistics
        """
        # Import here to avoid loading heavy modules at startup
        from app.services.rsnet.api import UniversalRSNetAPI
        
        api = UniversalRSNetAPI()
        return api.generate_reaction_network(**kwargs)
