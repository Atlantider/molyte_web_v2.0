"""
分子相似性搜索 API
支持基于中心分子的加权相似性搜索
"""
import logging
from typing import List, Dict, Optional
from fastapi import APIRouter, Depends, HTTPException, Query
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session

from app.database import get_db
from app.dependencies import get_current_active_user
from app.models.user import User
from app.services.predict_block_service import (
    find_similar_molecules, 
    predict_property,
    check_if_real_value
)
from app.services.molecule_info_service import get_molecule_info_fast

logger = logging.getLogger(__name__)
router = APIRouter()


class PropertyWeights(BaseModel):
    """属性权重配置"""
    melting_point: float = Field(default=1.0, ge=0, description="熔点权重")
    boiling_point: float = Field(default=1.0, ge=0, description="沸点权重") 
    flash_point: float = Field(default=1.0, ge=0, description="闪点权重")
    dipole_moment: float = Field(default=1.0, ge=0, description="偶极矩权重")
    polarizability: float = Field(default=1.0, ge=0, description="极化率权重")
    homo: float = Field(default=1.0, ge=0, description="HOMO权重")
    lumo: float = Field(default=1.0, ge=0, description="LUMO权重")
    gap: float = Field(default=1.0, ge=0, description="HOMO-LUMO能隙权重")


class SimilaritySearchRequest(BaseModel):
    """相似性搜索请求"""
    center_smiles: str = Field(..., description="中心分子的SMILES")
    n_similar: int = Field(default=20, ge=1, le=100, description="返回相似分子数量")
    weights: Optional[PropertyWeights] = Field(default=None, description="属性权重")
    include_predictions: bool = Field(default=True, description="是否包含预测值")


class MoleculeProperties(BaseModel):
    """分子属性"""
    melting_point: Optional[float] = None
    boiling_point: Optional[float] = None
    flash_point: Optional[float] = None
    dipole_moment: Optional[float] = None
    polarizability: Optional[float] = None
    homo: Optional[float] = None
    lumo: Optional[float] = None
    gap: Optional[float] = None
    
    # 标记哪些是真实值，哪些是预测值
    real_values: Dict[str, bool] = Field(default_factory=dict)


class SimilarMolecule(BaseModel):
    """相似分子信息"""
    smiles: str
    name: Optional[str] = None
    similarity_score: float = Field(..., description="相似度分数 (0-1)")
    weighted_distance: float = Field(..., description="加权距离")
    properties: MoleculeProperties
    is_center: bool = Field(default=False, description="是否为中心分子")


class SimilaritySearchResponse(BaseModel):
    """相似性搜索响应"""
    center_molecule: SimilarMolecule
    similar_molecules: List[SimilarMolecule]
    weights_used: PropertyWeights
    total_found: int


def _get_molecule_properties(smiles: str, include_predictions: bool = True) -> MoleculeProperties:
    """获取分子的所有属性（真实值+预测值）"""
    properties = MoleculeProperties()
    
    # 物理性质列表
    physical_props = ['bp', 'mp', 'fp']  # boiling_point, melting_point, flash_point
    
    for prop in physical_props:
        # 检查是否有真实值
        is_real, real_value = check_if_real_value(smiles, prop)
        
        if is_real and real_value is not None:
            # 使用真实值
            if prop == 'bp':
                properties.boiling_point = real_value
                properties.real_values['boiling_point'] = True
            elif prop == 'mp':
                properties.melting_point = real_value
                properties.real_values['melting_point'] = True
            elif prop == 'fp':
                properties.flash_point = real_value
                properties.real_values['flash_point'] = True
        elif include_predictions:
            # 使用预测值
            predicted_value = predict_property(smiles, prop)
            if predicted_value is not None:
                if prop == 'bp':
                    properties.boiling_point = predicted_value
                    properties.real_values['boiling_point'] = False
                elif prop == 'mp':
                    properties.melting_point = predicted_value
                    properties.real_values['melting_point'] = False
                elif prop == 'fp':
                    properties.flash_point = predicted_value
                    properties.real_values['flash_point'] = False
    
    # TODO: 添加从QC数据库获取HOMO/LUMO/GAP等量化性质的逻辑
    # 这里暂时使用模拟数据
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # 使用RDKit计算一些描述符作为占位符
            properties.dipole_moment = round(Descriptors.BertzCT(mol) / 100, 2)  # 占位符
            properties.polarizability = round(Descriptors.TPSA(mol) / 10, 2)  # 占位符
            properties.real_values['dipole_moment'] = False
            properties.real_values['polarizability'] = False
    except:
        pass
    
    return properties


def _calculate_weighted_distance(props1: MoleculeProperties, props2: MoleculeProperties, 
                                weights: PropertyWeights) -> float:
    """计算两个分子之间的加权欧几里得距离"""
    import numpy as np
    
    # 提取属性值，缺失值用0填充
    def get_prop_value(props: MoleculeProperties, prop_name: str) -> float:
        value = getattr(props, prop_name, None)
        return float(value) if value is not None else 0.0
    
    # 属性名称映射
    prop_names = [
        'melting_point', 'boiling_point', 'flash_point',
        'dipole_moment', 'polarizability', 'homo', 'lumo', 'gap'
    ]
    
    # 计算加权距离
    weighted_diff_squared = 0.0
    total_weight = 0.0
    
    for prop_name in prop_names:
        val1 = get_prop_value(props1, prop_name)
        val2 = get_prop_value(props2, prop_name)
        weight = getattr(weights, prop_name, 1.0)
        
        if val1 != 0.0 or val2 != 0.0:  # 只有当至少一个值非零时才计算
            diff = val1 - val2
            weighted_diff_squared += weight * weight * diff * diff
            total_weight += weight
    
    if total_weight == 0:
        return float('inf')
    
    return np.sqrt(weighted_diff_squared)


@router.post("/search", response_model=SimilaritySearchResponse)
def search_similar_molecules(
    request: SimilaritySearchRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    基于中心分子搜索相似分子
    
    支持多属性加权相似性计算，可以调整不同属性的重要性
    """
    try:
        # 设置默认权重
        weights = request.weights or PropertyWeights()
        
        # 获取中心分子的属性
        center_props = _get_molecule_properties(request.center_smiles, request.include_predictions)
        center_info = get_molecule_info_fast(request.center_smiles)
        
        # 使用现有的相似性搜索找到候选分子
        similar_candidates = find_similar_molecules(request.center_smiles, request.n_similar * 3)
        
        if not similar_candidates:
            raise HTTPException(status_code=404, detail="未找到相似分子")
        
        # 计算每个候选分子的加权距离
        similar_molecules = []
        
        for candidate in similar_candidates:
            candidate_smiles = candidate['smiles']
            candidate_props = _get_molecule_properties(candidate_smiles, request.include_predictions)
            candidate_info = get_molecule_info_fast(candidate_smiles)
            
            # 计算加权距离
            weighted_distance = _calculate_weighted_distance(center_props, candidate_props, weights)
            
            # 转换为相似度分数 (0-1)
            similarity_score = 1.0 / (1.0 + weighted_distance)
            
            similar_molecules.append(SimilarMolecule(
                smiles=candidate_smiles,
                name=candidate_info.get('name'),
                similarity_score=similarity_score,
                weighted_distance=weighted_distance,
                properties=candidate_props,
                is_center=candidate_smiles == request.center_smiles
            ))
        
        # 按相似度排序并取前N个
        similar_molecules.sort(key=lambda x: x.weighted_distance)
        similar_molecules = similar_molecules[:request.n_similar]
        
        # 创建中心分子对象
        center_molecule = SimilarMolecule(
            smiles=request.center_smiles,
            name=center_info.get('name'),
            similarity_score=1.0,
            weighted_distance=0.0,
            properties=center_props,
            is_center=True
        )
        
        return SimilaritySearchResponse(
            center_molecule=center_molecule,
            similar_molecules=similar_molecules,
            weights_used=weights,
            total_found=len(similar_molecules)
        )
        
    except Exception as e:
        logger.error(f"Similarity search failed: {e}")
        raise HTTPException(status_code=500, detail=f"相似性搜索失败: {str(e)}")


@router.get("/weights/presets")
def get_weight_presets():
    """获取预设的权重配置"""
    return {
        "thermal_properties": PropertyWeights(
            melting_point=3.0,
            boiling_point=3.0,
            flash_point=2.0,
            dipole_moment=1.0,
            polarizability=1.0,
            homo=0.5,
            lumo=0.5,
            gap=0.5
        ),
        "electronic_properties": PropertyWeights(
            melting_point=0.5,
            boiling_point=0.5,
            flash_point=0.5,
            dipole_moment=2.0,
            polarizability=2.0,
            homo=3.0,
            lumo=3.0,
            gap=3.0
        ),
        "balanced": PropertyWeights(
            melting_point=1.0,
            boiling_point=1.0,
            flash_point=1.0,
            dipole_moment=1.0,
            polarizability=1.0,
            homo=1.0,
            lumo=1.0,
            gap=1.0
        )
    }
