"""
引擎工厂 - 管理QC引擎实例
"""
from typing import Dict, Type
from app.services.qc_engines import QCEngine


class QCEngineFactory:
    """QC引擎工厂"""
    
    _engines: Dict[str, Type[QCEngine]] = {}
    
    @classmethod
    def register_engine(cls, name: str, engine_class: Type[QCEngine]):
        """注册引擎"""
        cls._engines[name.lower()] = engine_class
    
    @classmethod
    def get_engine(cls, engine_name: str) -> QCEngine:
        """
        获取QC引擎实例
        
        Args:
            engine_name: 引擎名称 (gaussian, pyscf)
            
        Returns:
            引擎实例
            
        Raises:
            ValueError: 未知引擎
        """
        engine_class = cls._engines.get(engine_name.lower())
        if not engine_class:
            raise ValueError(
                f"Unknown QC engine: {engine_name}. "
                f"Available engines: {list(cls._engines.keys())}"
            )
        return engine_class()
    
    @classmethod
    def list_engines(cls) -> list:
        """列出所有可用引擎"""
        return list(cls._engines.keys())


# 自动注册引擎
def _auto_register_engines():
    """自动注册所有可用引擎"""
    try:
        from app.services.qc_engines.gaussian import GaussianEngine
        QCEngineFactory.register_engine('gaussian', GaussianEngine)
    except ImportError:
        pass
    
    try:
        from app.services.qc_engines.pyscf import PySCFEngine  
        QCEngineFactory.register_engine('pyscf', PySCFEngine)
    except ImportError:
        pass


# 在模块加载时注册引擎
_auto_register_engines()
