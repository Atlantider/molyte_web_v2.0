"""
溶元AI挖掘 API 端点

提供分子预测和聚类浏览功能
"""

from fastapi import APIRouter, Depends, HTTPException, Query
from typing import List, Optional, Dict, Any
from pydantic import BaseModel
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import json
import os
from pathlib import Path
import logging

from app.models.user import User
from app.dependencies import get_current_active_user
from app.utils.permissions import require_module_access, MODULE_AI_DISCOVERY
from app.services.predict_block_service import (
    load_bp_data, load_mp_data, load_fp_data,
    predict_property, predict_property_detailed, predict_qm9_properties, find_similar_molecules, get_cluster_statistics,
    get_cluster_molecules, get_cluster_center_molecule, get_qm9_properties,
    get_similar_molecules_in_cluster, get_similar_molecules_with_weights, generate_molecule_image
)
# 使用predict_block_service中的load_cluster_data
from app.services.predict_block_service import load_cluster_data as load_cluster_data_from_service

router = APIRouter()
logger = logging.getLogger(__name__)

# 临时测试端点，不需要认证
@router.get("/test-similar-molecules")
async def test_similar_molecules(
    smiles: str = Query(..., description="SMILES字符串"),
    n_similar: int = Query(10, description="返回相似分子数量"),
    scope: str = Query("global", description="搜索范围")
):
    """测试相似分子搜索功能（无需认证）"""
    try:
        from app.services.predict_block_service import find_similar_molecules_by_qm9_distance

        similar_molecules = find_similar_molecules_by_qm9_distance(
            smiles=smiles,
            n_similar=n_similar,
            scope=scope
        )

        return {
            "query_smiles": smiles,
            "similar_molecules": similar_molecules,
            "count": len(similar_molecules)
        }
    except Exception as e:
        logger.error(f"Test similar molecules failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))

# 数据文件路径
PREDICT_BLOCK_PATH = Path("/opt/molyte_web_v1.0/predict_block-master")
DATASETS_PATH = PREDICT_BLOCK_PATH / "datasets"

# 常见溶剂分子的原子组成
SOLVENT_COMPOSITIONS = {
    "EC": {"C": 3, "H": 4, "O": 3},
    "DMC": {"C": 3, "H": 6, "O": 3},
    "DEC": {"C": 5, "H": 10, "O": 3},
    "EMC": {"C": 4, "H": 8, "O": 3},
    "PC": {"C": 4, "H": 6, "O": 3},
    "VC": {"C": 3, "H": 2, "O": 3},
    "FEC": {"C": 3, "H": 3, "O": 3, "F": 1},
    "PP": {"C": 6, "H": 12, "O": 2},
    "MP": {"C": 4, "H": 8, "O": 2},
    "ACN": {"C": 2, "H": 3, "N": 1},
    "DMF": {"C": 3, "H": 7, "N": 1, "O": 1},
    "NMP": {"C": 5, "H": 9, "N": 1, "O": 1},
    "DME": {"C": 4, "H": 10, "O": 2},
    "DOL": {"C": 3, "H": 6, "O": 2},
    "THF": {"C": 4, "H": 8, "O": 1},
    "DEE": {"C": 4, "H": 10, "O": 1},
    "1,4-Dioxane": {"C": 4, "H": 8, "O": 2},
    "SN": {"C": 4, "H": 4, "N": 2},
    "ADN": {"C": 6, "H": 8, "N": 2},
    "Acrylonitrile": {"C": 3, "H": 3, "N": 1},
    "DMAc": {"C": 4, "H": 9, "N": 1, "O": 1},
    "GBL": {"C": 4, "H": 6, "O": 2}
}

# ============================================================================
# Pydantic Models
# ============================================================================

class MoleculeInfo(BaseModel):
    """分子信息"""
    name: str
    smiles: Optional[str] = None
    composition: Dict[str, int]
    cluster_id: Optional[int] = None
    properties: Optional[Dict[str, Any]] = None
    cas_number: Optional[str] = None  # CAS 号
    molecular_formula: Optional[str] = None  # 分子式
    molecular_weight: Optional[float] = None  # 分子量
    is_real_data: Optional[Dict[str, bool]] = None  # 每个属性是否为真实值
    similarity: Optional[float] = None  # 相似度分数
    image: Optional[str] = None  # base64编码的分子结构图像

class ClusterInfo(BaseModel):
    """聚类信息"""
    cluster_id: int
    size: int
    molecules: List[str]
    center_properties: Optional[Dict[str, float]] = None
    description: Optional[str] = None

class PredictionRequest(BaseModel):
    """分子预测请求"""
    smiles: str
    properties: List[str] = ["bp", "mp", "fp", "alpha", "mu", "gap", "homo", "lumo"]  # 8种属性
    include_qm9: bool = True  # 是否包含 QM9 量子化学属性
    dark_mode: bool = False  # 是否为暗色模式（影响分子图像生成）

class BatchPredictionRequest(BaseModel):
    """批量分子预测请求"""
    smiles_list: List[str]
    properties: List[str] = ["bp", "mp", "fp", "alpha", "mu", "gap", "homo", "lumo"]  # 8种属性
    include_qm9: bool = True  # 是否包含 QM9 量子化学属性
    dark_mode: bool = False  # 是否为暗色模式（影响分子图像生成）

class PredictionResponse(BaseModel):
    """分子预测响应"""
    smiles: str
    predicted_properties: Dict[str, float]
    confidence: Optional[Dict[str, float]] = None
    similar_molecules: Optional[List[MoleculeInfo]] = None
    molecule_info: Optional[Dict[str, Any]] = None  # 主分子的信息
    is_real_data: Optional[Dict[str, bool]] = None  # 每个属性是否为真实值
    image: Optional[str] = None  # base64编码的分子结构图像

class BatchPredictionResponse(BaseModel):
    """批量分子预测响应"""
    results: List[PredictionResponse]
    total_count: int
    success_count: int
    failed_smiles: List[str] = []

class CenterMolecule(BaseModel):
    """中心分子（用于聚类浏览）"""
    smiles: str
    name: Optional[str] = None
    cluster_id: Optional[int] = None
    properties: Optional[Dict[str, float]] = None

class TemplatesResponse(BaseModel):
    """中心分子列表响应"""
    templates: List[CenterMolecule]
    total_count: int = 0

class ClusterAnalysisResponse(BaseModel):
    """聚类分析响应"""
    total_clusters: int
    clusters: List[ClusterInfo]
    feature_importance: Optional[Dict[str, float]] = None

class WeightedSimilarityRequest(BaseModel):
    """基于权重的相似分子搜索请求"""
    cluster_id: int
    center_mol_index: int
    n_similar: int = 10
    weights: Optional[Dict[str, float]] = None  # 特征权重，如 {"mu": 1.0, "alpha": 0.5, ...}

# ============================================================================
# Helper Functions
# ============================================================================

def count_atoms_from_sdf(sdf_text: str) -> Dict[str, int]:
    """从SDF文本中统计原子数量"""
    lines = sdf_text.split("\n")
    atom_counts = {}
    for line in lines[3:]:  # SDF 块第 4 行开始是原子坐标
        if line.strip() == "" or line.startswith("M  END"):
            break
        parts = line.split()
        if len(parts) < 4:
            continue
        atom = parts[3]
        atom_counts[atom] = atom_counts.get(atom, 0) + 1
    return atom_counts



def find_molecule_cluster(composition: Dict[str, int]) -> Optional[int]:
    """根据原子组成查找分子所属的聚类"""
    try:
        df = load_cluster_data()
        if df is None:
            return None

        for idx, row in df.iterrows():
            counts = count_atoms_from_sdf(row["mol_expr"])
            if all(counts.get(atom, 0) == num for atom, num in composition.items()):
                return int(row['cluster'])
        return None
    except:
        return None

def get_mock_cluster_analysis() -> ClusterAnalysisResponse:
    """返回模拟的聚类分析数据"""
    clusters = []

    # 创建一些模拟聚类
    cluster_data = [
        {"id": 1, "size": 15, "molecules": ["EC", "DMC", "DEC"], "desc": "碳酸酯类溶剂"},
        {"id": 2, "size": 12, "molecules": ["EMC", "PC", "VC"], "desc": "环状碳酸酯类"},
        {"id": 3, "size": 8, "molecules": ["FEC", "DFEC", "TTE"], "desc": "氟化溶剂"},
        {"id": 4, "size": 6, "molecules": ["DME", "DOL", "TEGDME"], "desc": "醚类溶剂"},
        {"id": 5, "size": 4, "molecules": ["AN", "SN", "GBL"], "desc": "其他有机溶剂"}
    ]

    for cluster in cluster_data:
        clusters.append(ClusterInfo(
            cluster_id=cluster["id"],
            size=cluster["size"],
            molecules=cluster["molecules"],
            center_properties={
                "mu": np.random.uniform(1.0, 5.0),
                "alpha": np.random.uniform(5.0, 15.0),
                "homo": np.random.uniform(-0.3, -0.1),
                "lumo": np.random.uniform(0.1, 0.3),
                "gap": np.random.uniform(0.2, 0.5)
            },
            description=cluster["desc"]
        ))

    return ClusterAnalysisResponse(
        total_clusters=len(clusters),
        clusters=clusters,
        feature_importance={
            "mu": 0.2,
            "alpha": 0.15,
            "homo": 0.25,
            "lumo": 0.25,
            "gap": 0.15
        }
    )

# ============================================================================
# API Endpoints
# ============================================================================

@router.get("/templates", response_model=TemplatesResponse)
def get_molecule_templates(
    limit: int = Query(50, description="返回的中心分子数，默认50个"),
    search: Optional[str] = Query(None, description="搜索关键词（可选）"),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取中心分子列表（用于聚类浏览）

    返回聚类数据中的代表性分子，用户可以选择其中一个作为中心分子进行聚类探索
    """
    require_module_access(current_user, MODULE_AI_DISCOVERY)

    try:
        cluster_data = load_cluster_data_from_service()
        if cluster_data is None or cluster_data.empty:
            logger.warning("Cluster data is empty")
            return TemplatesResponse(templates=[], total_count=0)

        logger.info(f"Loaded cluster data with {len(cluster_data)} rows")

        # 从每个聚类中选择一个代表分子（第一个）
        templates = []
        seen_clusters = set()

        for idx, (_, row) in enumerate(cluster_data.iterrows()):
            cluster_id = row.get('cluster')

            # 每个聚类只选择一个代表分子
            if pd.notna(cluster_id) and cluster_id in seen_clusters:
                continue

            if pd.notna(cluster_id):
                seen_clusters.add(cluster_id)

            smiles = row.get('SMILES')
            if not smiles:
                continue

            try:
                # 提取QM9属性
                qm9_features = ['mu', 'alpha', 'homo', 'lumo', 'gap']
                properties = {}
                for feature in qm9_features:
                    if feature in row and pd.notna(row[feature]):
                        properties[feature] = float(row[feature])

                # 使用SMILES的前30个字符作为名称（避免调用get_molecule_info_fast导致超时）
                name = smiles[:30] + ('...' if len(smiles) > 30 else '')

                template = CenterMolecule(
                    smiles=smiles,
                    name=name,
                    cluster_id=int(cluster_id) if pd.notna(cluster_id) else None,
                    properties=properties if properties else None
                )

                templates.append(template)
                logger.debug(f"Added template: {smiles}")

                # 达到限制数量后停止
                if len(templates) >= limit:
                    break
            except Exception as e:
                logger.warning(f"Failed to process SMILES {smiles}: {e}")
                continue

        logger.info(f"Returning {len(templates)} templates")
        return TemplatesResponse(
            templates=templates,
            total_count=len(templates)
        )

    except Exception as e:
        logger.error(f"Failed to get center molecules: {e}", exc_info=True)
        raise HTTPException(
            status_code=500,
            detail=f"获取中心分子列表失败: {str(e)}"
        )

@router.get("/clusters", response_model=ClusterAnalysisResponse)
def get_cluster_analysis(
    n_clusters: int = Query(50, description="聚类数量"),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取分子聚类分析结果

    基于分子的物理化学性质进行聚类分析
    """
    require_module_access(current_user, MODULE_AI_DISCOVERY)

    try:
        df = load_cluster_data()
        if df is None:
            # 如果数据文件不存在，返回模拟数据
            logger.warning("Cluster data not found, returning mock data")
            return get_mock_cluster_analysis()
    except Exception as e:
        logger.error(f"Failed to load cluster data: {e}")
        # 如果加载失败，返回模拟数据
        return get_mock_cluster_analysis()

    try:
        # 统计每个聚类的信息
        clusters = []
        cluster_stats = df.groupby('cluster').size().reset_index(name='size')
        cluster_stats = cluster_stats.sort_values('size', ascending=False)

        # 显示前n_clusters个聚类
        for _, row in cluster_stats.head(n_clusters).iterrows():
            cluster_id = int(row['cluster'])
            size = int(row['size'])

            # 获取该聚类中的分子
            cluster_molecules = df[df['cluster'] == cluster_id]
            molecules = [f"Mol_{i+1}" for i in range(min(size, 5))]  # 最多显示5个

            # 计算聚类中心的性质
            features = ["mu", "alpha", "homo", "lumo", "gap"]
            center_props = {}

            for feature in features:
                if feature in cluster_molecules.columns:
                    center_props[feature] = float(cluster_molecules[feature].mean())

            clusters.append(ClusterInfo(
                cluster_id=cluster_id,
                size=size,
                molecules=molecules,
                center_properties=center_props,
                description=f"聚类 {cluster_id} - {size} 个分子"
            ))

        return ClusterAnalysisResponse(
            total_clusters=int(df['cluster'].max()),
            clusters=clusters,
            feature_importance={
                "mu": 0.2,
                "alpha": 0.15,
                "homo": 0.25,
                "lumo": 0.25,
                "gap": 0.15
            }
        )

    except Exception as e:
        logger.error(f"Error in cluster analysis: {e}")
        # 如果处理过程中出错，返回模拟数据
        return get_mock_cluster_analysis()

@router.get("/cluster-stats")
def get_cluster_stats(
    current_user: User = Depends(get_current_active_user)
):
    """
    获取聚类统计信息

    返回聚类的总数、分子总数、最大聚类大小、平均聚类大小等统计信息
    """
    require_module_access(current_user, MODULE_AI_DISCOVERY)

    try:
        stats = get_cluster_statistics()
        if stats is None:
            raise HTTPException(
                status_code=500,
                detail="无法获取聚类统计信息"
            )
        return stats
    except Exception as e:
        logger.error(f"Failed to get cluster statistics: {e}")
        raise HTTPException(
            status_code=500,
            detail=f"获取聚类统计信息失败: {str(e)}"
        )


@router.get("/cluster/{cluster_id}/molecules")
def get_cluster_molecules_endpoint(
    cluster_id: int,
    limit: int = Query(100, description="返回的最大分子数"),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取聚类中的分子列表

    返回指定聚类中的分子及其属性信息
    """
    require_module_access(current_user, MODULE_AI_DISCOVERY)

    try:
        molecules = get_cluster_molecules(cluster_id, limit=limit)
        if molecules is None:
            raise HTTPException(
                status_code=404,
                detail=f"聚类 {cluster_id} 不存在或无分子数据"
            )

        return {
            'cluster_id': cluster_id,
            'molecule_count': len(molecules),
            'molecules': molecules
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to get cluster molecules for cluster {cluster_id}: {e}")
        raise HTTPException(
            status_code=500,
            detail=f"获取聚类分子失败: {str(e)}"
        )


@router.get("/cluster/{cluster_id}/center")
def get_cluster_center_endpoint(
    cluster_id: int,
    current_user: User = Depends(get_current_active_user)
):
    """
    获取聚类的中心分子

    返回聚类中所有分子属性的平均值作为中心分子
    """
    require_module_access(current_user, MODULE_AI_DISCOVERY)

    try:
        center = get_cluster_center_molecule(cluster_id)
        if center is None:
            raise HTTPException(
                status_code=404,
                detail=f"聚类 {cluster_id} 不存在"
            )

        return center
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to get cluster center for cluster {cluster_id}: {e}")
        raise HTTPException(
            status_code=500,
            detail=f"获取聚类中心失败: {str(e)}"
        )


def validate_smiles(smiles: str) -> bool:
    """简单的SMILES格式验证"""
    if not smiles or len(smiles.strip()) == 0:
        return False

    # 基本的SMILES字符检查
    valid_chars = set("CNOPSFClBrI()[]=-#@+123456789")
    return all(c in valid_chars for c in smiles.upper())

@router.post("/predict", response_model=PredictionResponse)
def predict_molecule_properties(
    request: PredictionRequest,
    current_user: User = Depends(get_current_active_user)
):
    """
    预测分子性质

    基于SMILES输入预测分子的物理化学性质
    """
    require_module_access(current_user, MODULE_AI_DISCOVERY)

    # 验证SMILES格式
    if not validate_smiles(request.smiles):
        raise HTTPException(
            status_code=400,
            detail="无效的SMILES格式，请检查输入"
        )

    try:
        # 使用predict_block的数据进行预测
        predicted_props = {}
        confidence = {}
        is_real_data = {}

        for prop in request.properties:
            # 使用详细预测函数
            detailed_result = predict_property_detailed(request.smiles, prop)

            if detailed_result is not None:
                predicted_props[prop] = float(detailed_result['value'])
                is_real_data[prop] = detailed_result['is_real']
                # 真实值的置信度更高
                confidence[prop] = 0.95 if detailed_result['is_real'] else 0.75
            else:
                # 如果预测失败，返回0.0作为默认值
                logger.warning(f"Failed to predict {prop} for {request.smiles}")
                predicted_props[prop] = 0.0
                confidence[prop] = 0.0
                is_real_data[prop] = False

        # 如果请求包含 QM9 属性，添加到预测结果中
        if request.include_qm9:
            qm9_props = predict_qm9_properties(request.smiles)
            if qm9_props:
                predicted_props.update(qm9_props)
                # 为 QM9 属性设置置信度
                for qm9_prop in qm9_props.keys():
                    confidence[qm9_prop] = 0.5  # QM9 属性的置信度较低

        # 查找相似分子（预测时返回10个相似分子）
        similar_molecules = []
        similar_data = find_similar_molecules(request.smiles, n_similar=10)

        for sim_mol in similar_data:
            # 获取分子信息
            from app.services.molecule_info_service import get_molecule_info_fast
            mol_info = get_molecule_info_fast(sim_mol['smiles'])

            # 检查相似分子的属性是否为真实值
            sim_is_real_data = {}
            for prop in ['bp', 'mp', 'fp']:
                if prop in sim_mol['properties']:
                    from app.services.predict_block_service import check_if_real_value
                    is_real, _ = check_if_real_value(sim_mol['smiles'], prop)
                    sim_is_real_data[prop] = is_real

            similar_molecules.append(MoleculeInfo(
                name=mol_info.get('name', f"Molecule_{sim_mol['smiles'][:10]}"),
                smiles=sim_mol['smiles'],
                composition={},  # 简化处理
                properties=sim_mol['properties'],  # 包含所有可用的属性
                similarity=sim_mol['similarity'],
                cas_number=mol_info.get('cas_number'),
                molecular_formula=mol_info.get('molecular_formula'),
                molecular_weight=mol_info.get('molecular_weight'),
                is_real_data=sim_is_real_data
            ))

        # 获取主分子的信息
        from app.services.molecule_info_service import get_molecule_info_fast
        main_molecule_info = get_molecule_info_fast(request.smiles)

        # 生成分子结构图像（支持dark模式）
        dark_mode = request.dark_mode if hasattr(request, 'dark_mode') else False
        molecule_image = generate_molecule_image(request.smiles, dark_mode=dark_mode)

        return PredictionResponse(
            smiles=request.smiles,
            predicted_properties=predicted_props,
            confidence=confidence,
            similar_molecules=similar_molecules,
            molecule_info=main_molecule_info,
            is_real_data=is_real_data,
            image=molecule_image  # 添加分子结构图像
        )

    except Exception as e:
        logger.error(f"Prediction error: {e}")
        raise HTTPException(
            status_code=500,
            detail=f"预测过程中发生错误: {str(e)}"
        )

@router.post("/predict-batch", response_model=BatchPredictionResponse)
def predict_molecules_batch(
    request: BatchPredictionRequest,
    current_user: User = Depends(get_current_active_user)
):
    """
    批量预测分子性质

    一次性预测多个分子的物理化学性质，提高效率
    """
    require_module_access(current_user, MODULE_AI_DISCOVERY)

    results = []
    failed_smiles = []
    success_count = 0

    for smiles in request.smiles_list:
        try:
            # 验证SMILES格式
            if not validate_smiles(smiles):
                failed_smiles.append(smiles)
                continue

            # 创建单个预测请求
            single_request = PredictionRequest(
                smiles=smiles,
                properties=request.properties,
                include_qm9=request.include_qm9
            )

            # 使用现有的预测逻辑
            predicted_props = {}
            confidence = {}
            is_real_data = {}

            for prop in request.properties:
                detailed_result = predict_property_detailed(smiles, prop)
                if detailed_result is not None:
                    predicted_props[prop] = float(detailed_result['value'])
                    is_real_data[prop] = detailed_result['is_real']
                    confidence[prop] = 0.95 if detailed_result['is_real'] else 0.75
                else:
                    predicted_props[prop] = 0.0
                    confidence[prop] = 0.0
                    is_real_data[prop] = False

            # QM9 属性
            if request.include_qm9:
                qm9_props = predict_qm9_properties(smiles)
                if qm9_props:
                    predicted_props.update(qm9_props)
                    for qm9_prop in qm9_props.keys():
                        confidence[qm9_prop] = 0.5

            # 获取分子信息（不包含相似分子以提高速度）
            from app.services.molecule_info_service import get_molecule_info_fast, COMMON_MOLECULES

            # 批量预测时只为常见分子查询详细信息，提高速度
            if smiles in COMMON_MOLECULES:
                main_molecule_info = get_molecule_info_fast(smiles)
            else:
                # 为不常见分子生成基本信息
                from rdkit import Chem
                try:
                    rdkit_mol = Chem.MolFromSmiles(smiles)
                    if rdkit_mol:
                        mol_formula = Chem.rdMolDescriptors.CalcMolFormula(rdkit_mol)
                        mol_weight = Chem.rdMolDescriptors.CalcExactMolWt(rdkit_mol)
                        main_molecule_info = {
                            'smiles': smiles,
                            'name': mol_formula,  # 使用分子式作为名称
                            'molecular_formula': mol_formula,
                            'molecular_weight': mol_weight,
                            'cas_number': None
                        }
                    else:
                        main_molecule_info = {
                            'smiles': smiles,
                            'name': f"Molecule_{smiles[:10]}",
                            'molecular_formula': 'Unknown',
                            'molecular_weight': None,
                            'cas_number': None
                        }
                except:
                    main_molecule_info = {
                        'smiles': smiles,
                        'name': f"Molecule_{smiles[:10]}",
                        'molecular_formula': 'Unknown',
                        'molecular_weight': None,
                        'cas_number': None
                    }

            # 生成分子结构图像（支持dark模式）
            dark_mode = request.dark_mode if hasattr(request, 'dark_mode') else False
            molecule_image = generate_molecule_image(smiles, dark_mode=dark_mode)

            results.append(PredictionResponse(
                smiles=smiles,
                predicted_properties=predicted_props,
                confidence=confidence,
                similar_molecules=[],  # 批量预测时不包含相似分子
                molecule_info=main_molecule_info,
                is_real_data=is_real_data,
                image=molecule_image  # 添加分子结构图像
            ))
            success_count += 1

        except Exception as e:
            logger.error(f"Failed to predict {smiles}: {e}")
            failed_smiles.append(smiles)

    return BatchPredictionResponse(
        results=results,
        total_count=len(request.smiles_list),
        success_count=success_count,
        failed_smiles=failed_smiles
    )

class PropertyOptimizationRequest(BaseModel):
    """属性优化搜索请求"""
    smiles: str
    n_results: int = 50
    scope: str = "global"  # 'cluster' 或 'global'
    # 优化目标：属性名称 -> (目标值, 权重)
    # 例如: {"gap": (5.0, 1.0), "bp": (100, 0.5)} 表示希望gap >= 5.0，bp <= 100
    optimization_targets: Optional[Dict[str, Dict[str, float]]] = None
    # 相似度权重 (0-1)
    similarity_weight: float = 0.5

class MultiObjectiveOptimizationRequest(BaseModel):
    """多目标优化请求"""
    smiles: str
    n_results: int = 100  # 增加候选分子数量
    scope: str = "global"  # 'cluster' 或 'global'
    optimization_targets: Dict[str, Dict[str, Any]]  # 优化目标配置
    similarity_weight: float = 0.3
    optimization_method: str = "pareto"  # "pareto" 或 "weighted"

@router.post("/optimize-molecules")
def optimize_molecules_search(
    request: PropertyOptimizationRequest,
    current_user: User = Depends(get_current_active_user)
):
    """
    属性优化搜索

    找到既与中心分子相似，又具有更好属性的分子
    支持多个属性的同时优化

    所有认证用户都可以使用此功能
    """

    try:
        from app.services.predict_block_service import find_similar_molecules_by_qm9_distance

        # 首先找到相似分子
        similar_molecules = find_similar_molecules_by_qm9_distance(
            request.smiles,
            n_similar=request.n_results * 2,  # 获取更多候选
            scope=request.scope
        )

        if not similar_molecules:
            return {
                'center_molecule': {'smiles': request.smiles},
                'optimized_molecules': [],
                'optimization_targets': request.optimization_targets or {}
            }

        # 获取中心分子的属性
        center_mol = similar_molecules[0] if similar_molecules else None
        center_props = center_mol['properties'] if center_mol else {}

        # 根据优化目标筛选和排序
        optimized_molecules = []

        for mol in similar_molecules:
            mol_props = mol['properties']
            optimization_score = 0.0
            improvements = {}

            if request.optimization_targets:
                for prop_name, target_info in request.optimization_targets.items():
                    if isinstance(target_info, dict):
                        target_value = target_info.get('target', 0)
                        weight = target_info.get('weight', 1.0)
                        direction = target_info.get('direction', 'higher')  # 'higher' 或 'lower'
                    else:
                        target_value = target_info
                        weight = 1.0
                        direction = 'higher'

                    if prop_name in mol_props and prop_name in center_props:
                        mol_val = mol_props[prop_name]
                        center_val = center_props[prop_name]

                        if direction == 'higher':
                            # 目标是更高的值
                            if mol_val > center_val:
                                improvement = (mol_val - center_val) / abs(center_val) if center_val != 0 else 0
                                improvements[prop_name] = improvement
                                optimization_score += improvement * weight
                        else:
                            # 目标是更低的值
                            if mol_val < center_val:
                                improvement = (center_val - mol_val) / abs(center_val) if center_val != 0 else 0
                                improvements[prop_name] = improvement
                                optimization_score += improvement * weight

            # 结合相似度和优化分数
            combined_score = (mol['similarity'] * request.similarity_weight +
                            optimization_score * (1 - request.similarity_weight))

            optimized_molecules.append({
                'smiles': mol['smiles'],
                'properties': mol_props,
                'similarity': mol['similarity'],
                'optimization_score': optimization_score,
                'combined_score': combined_score,
                'improvements': improvements
            })

        # 按综合分数排序
        optimized_molecules.sort(key=lambda x: x['combined_score'], reverse=True)
        optimized_molecules = optimized_molecules[:request.n_results]

        return {
            'center_molecule': {
                'smiles': request.smiles,
                'properties': center_props
            },
            'optimized_molecules': optimized_molecules,
            'optimization_targets': request.optimization_targets or {}
        }

    except Exception as e:
        logger.error(f"Optimization search error: {e}")
        raise HTTPException(
            status_code=500,
            detail=f"属性优化搜索失败: {str(e)}"
        )


@router.post("/multi-objective-optimize")
def multi_objective_optimization(
    request: MultiObjectiveOptimizationRequest,
    current_user: User = Depends(get_current_active_user)
):
    """
    多目标优化搜索

    支持Pareto前沿分析和加权优化两种方法
    返回更多候选分子以提供更好的优化选择
    """
    try:
        from app.services.predict_block_service import find_similar_molecules_by_qm9_distance
        import numpy as np

        # 获取候选分子 - 限制数量以提高性能
        max_candidates = min(request.n_results * 2, 100)  # 最多100个候选分子
        similar_molecules = find_similar_molecules_by_qm9_distance(
            request.smiles,
            n_similar=max_candidates,
            scope=request.scope
        )

        if not similar_molecules:
            return {
                'center_molecule': {'smiles': request.smiles},
                'candidates': [],
                'pareto_frontier': [],
                'optimization_targets': request.optimization_targets,
                'method': request.optimization_method
            }

        # 获取中心分子属性
        center_mol = similar_molecules[0] if similar_molecules else None
        center_props = center_mol['properties'] if center_mol else {}

        # 计算每个候选分子的目标函数值
        candidates = []
        for mol in similar_molecules:
            mol_props = mol['properties']
            objective_values = {}

            # 计算每个目标的值
            for prop_name, target_config in request.optimization_targets.items():
                if prop_name in mol_props:
                    direction = target_config.get('direction', 'up')
                    weight = target_config.get('weight', 1.0)

                    # 标准化目标值（相对于中心分子的改进）
                    if prop_name in center_props and center_props[prop_name] != 0:
                        center_val = center_props[prop_name]
                        mol_val = mol_props[prop_name]

                        if direction == 'up':
                            # 希望值更高
                            improvement = (mol_val - center_val) / abs(center_val)
                        else:
                            # 希望值更低
                            improvement = (center_val - mol_val) / abs(center_val)

                        objective_values[prop_name] = improvement
                    else:
                        objective_values[prop_name] = 0.0

            # 获取分子详细信息
            from app.services.molecule_info_service import get_molecule_info_fast
            molecule_info = get_molecule_info_fast(mol['smiles'])

            candidates.append({
                'smiles': mol['smiles'],
                'properties': mol_props,
                'similarity': mol['similarity'],
                'objective_values': objective_values,
                'molecule_info': molecule_info  # 添加分子信息
            })

        # 根据优化方法处理
        if request.optimization_method == "pareto":
            # Pareto前沿分析
            pareto_frontier = find_pareto_frontier(candidates, request.optimization_targets)
            result_molecules = pareto_frontier[:request.n_results]
        else:
            # 加权优化
            for candidate in candidates:
                weighted_score = 0.0
                for prop_name, target_config in request.optimization_targets.items():
                    weight = target_config.get('weight', 1.0)
                    if prop_name in candidate['objective_values']:
                        weighted_score += candidate['objective_values'][prop_name] * weight
                candidate['weighted_score'] = weighted_score

            # 按加权分数排序
            candidates.sort(key=lambda x: x.get('weighted_score', 0), reverse=True)
            result_molecules = candidates[:request.n_results]
            pareto_frontier = []

        return {
            'center_molecule': {
                'smiles': request.smiles,
                'properties': center_props
            },
            'candidates': result_molecules,
            'pareto_frontier': pareto_frontier if request.optimization_method == "pareto" else [],
            'optimization_targets': request.optimization_targets,
            'method': request.optimization_method,
            'total_candidates': len(similar_molecules)
        }

    except Exception as e:
        logger.error(f"Multi-objective optimization error: {e}")
        raise HTTPException(
            status_code=500,
            detail=f"多目标优化失败: {str(e)}"
        )


def find_pareto_frontier(candidates, optimization_targets):
    """找到Pareto前沿"""
    pareto_solutions = []

    for i, candidate in enumerate(candidates):
        is_dominated = False

        # 检查是否被其他解支配
        for j, other in enumerate(candidates):
            if i == j:
                continue

            dominates = True
            strictly_better = False

            # 检查所有目标
            for prop_name in optimization_targets.keys():
                if prop_name in candidate['objective_values'] and prop_name in other['objective_values']:
                    candidate_val = candidate['objective_values'][prop_name]
                    other_val = other['objective_values'][prop_name]

                    if other_val > candidate_val:
                        strictly_better = True
                    elif other_val < candidate_val:
                        dominates = False
                        break

            if dominates and strictly_better:
                is_dominated = True
                break

        if not is_dominated:
            pareto_solutions.append(candidate)

    return pareto_solutions


@router.get("/similar-molecules")
def get_similar_molecules_by_smiles(
    smiles: str = Query(..., description="查询分子的 SMILES"),
    n_similar: int = Query(50, description="返回的相似分子数，默认50个"),
    scope: str = Query("global", description="搜索范围: 'cluster' 或 'global'"),
    current_user: User = Depends(get_current_active_user)
):
    """
    根据 SMILES 获取相似分子

    基于QM9属性（alpha, mu, homo, lumo, gap）的欧氏距离查找相似分子
    支持两种scope：
    - 'cluster': 在同一聚类内查找相似分子
    - 'global': 在全局范围内查找相似分子
    """
    require_module_access(current_user, MODULE_AI_DISCOVERY)

    try:
        # 使用QM9属性距离查找相似分子
        from app.services.predict_block_service import find_similar_molecules_by_qm9_distance

        similar_molecules = find_similar_molecules_by_qm9_distance(
            smiles,
            n_similar=n_similar,
            scope=scope
        )

        # 为每个相似分子添加详细信息（优化版本）
        enhanced_molecules = []
        for mol in similar_molecules:
            from app.services.molecule_info_service import get_molecule_info_fast, COMMON_MOLECULES
            from app.services.predict_block_service import check_if_real_value

            # 只为常见分子查询详细信息，其他分子使用基本信息
            if mol['smiles'] in COMMON_MOLECULES:
                mol_info = get_molecule_info_fast(mol['smiles'])
            else:
                # 为不常见分子生成基本信息，避免慢速PubChem查询
                from rdkit import Chem
                try:
                    rdkit_mol = Chem.MolFromSmiles(mol['smiles'])
                    if rdkit_mol:
                        mol_formula = Chem.rdMolDescriptors.CalcMolFormula(rdkit_mol)
                        mol_weight = Chem.rdMolDescriptors.CalcExactMolWt(rdkit_mol)
                        mol_info = {
                            'smiles': mol['smiles'],
                            'name': f"Molecule_{mol['smiles'][:10]}",
                            'molecular_formula': mol_formula,
                            'molecular_weight': mol_weight,
                            'cas_number': None
                        }
                    else:
                        mol_info = {
                            'smiles': mol['smiles'],
                            'name': f"Molecule_{mol['smiles'][:10]}",
                            'molecular_formula': 'Unknown',
                            'molecular_weight': None,
                            'cas_number': None
                        }
                except:
                    mol_info = {
                        'smiles': mol['smiles'],
                        'name': f"Molecule_{mol['smiles'][:10]}",
                        'molecular_formula': 'Unknown',
                        'molecular_weight': None,
                        'cas_number': None
                    }

            # 检查属性是否为真实值
            is_real_data = {}
            for prop in ['bp', 'mp', 'fp', 'alpha', 'mu', 'gap', 'homo', 'lumo']:
                if prop in mol['properties']:
                    is_real, _ = check_if_real_value(mol['smiles'], prop)
                    is_real_data[prop] = is_real

            enhanced_molecules.append({
                'smiles': mol['smiles'],
                'properties': mol['properties'],
                'similarity': mol['similarity'],
                'molecule_info': mol_info,
                'is_real_data': is_real_data
            })

        return {
            "smiles": smiles,
            "similar_molecules": enhanced_molecules,
            "count": len(enhanced_molecules)
        }

    except Exception as e:
        logger.error(f"Failed to find similar molecules for {smiles}: {e}")
        raise HTTPException(
            status_code=500,
            detail=f"查找相似分子失败: {str(e)}"
        )

@router.get("/clusters/{cluster_id}")
def get_cluster_detail(
    cluster_id: int,
    current_user: User = Depends(get_current_active_user)
):
    """获取特定聚类的详细信息"""
    require_module_access(current_user, MODULE_AI_DISCOVERY)

    df = load_cluster_data()
    if df is None:
        raise HTTPException(status_code=404, detail="聚类数据文件不存在")

    cluster_data = df[df['cluster'] == cluster_id]
    if cluster_data.empty:
        raise HTTPException(status_code=404, detail=f"聚类 {cluster_id} 不存在")

    # 返回聚类详细信息
    features = ["mu", "alpha", "homo", "lumo", "gap"]
    molecules = []

    # 只返回前100个分子以避免数据过大
    for idx, (_, row) in enumerate(cluster_data.head(100).iterrows()):
        mol_props = {feature: float(row[feature]) for feature in features if feature in row and pd.notna(row[feature])}
        molecules.append({
            "index": int(idx),
            "properties": mol_props
        })

    return {
        "cluster_id": cluster_id,
        "size": len(cluster_data),
        "molecules": molecules,
        "statistics": {
            feature: {
                "mean": float(cluster_data[feature].mean()),
                "std": float(cluster_data[feature].std()),
                "min": float(cluster_data[feature].min()),
                "max": float(cluster_data[feature].max())
            } for feature in features if feature in cluster_data.columns
        }
    }


@router.get("/cluster/{cluster_id}/similar-molecules")
def get_cluster_similar_molecules(
    cluster_id: int,
    center_index: int = Query(0, description="中心分子在聚类中的索引"),
    n_similar: int = Query(10, description="返回的相似分子数"),
    current_user: User = Depends(get_current_active_user)
):
    """
    获取聚类中与指定分子相似的分子

    基于 QM9 属性计算相似性，返回最相似的分子列表
    """
    require_module_access(current_user, MODULE_AI_DISCOVERY)

    try:
        result = get_similar_molecules_in_cluster(cluster_id, center_index, n_similar=n_similar)
        if result is None:
            raise HTTPException(
                status_code=404,
                detail=f"聚类 {cluster_id} 或分子索引 {center_index} 不存在"
            )

        return result

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to get similar molecules in cluster {cluster_id}: {e}")
        raise HTTPException(
            status_code=500,
            detail=f"获取相似分子失败: {str(e)}"
        )


@router.post("/cluster/{cluster_id}/similar-molecules-weighted")
def get_cluster_similar_molecules_weighted(
    cluster_id: int,
    request: WeightedSimilarityRequest,
    current_user: User = Depends(get_current_active_user)
):
    """
    获取聚类中与指定分子相似的分子（支持自定义权重）

    基于 QM9 属性和自定义权重计算相似性，返回最相似的分子列表
    """
    require_module_access(current_user, MODULE_AI_DISCOVERY)

    try:
        result = get_similar_molecules_with_weights(
            request.cluster_id,
            request.center_mol_index,
            n_similar=request.n_similar,
            weights=request.weights
        )
        if result is None:
            raise HTTPException(
                status_code=404,
                detail=f"聚类 {request.cluster_id} 或分子索引 {request.center_mol_index} 不存在"
            )

        return result

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to get similar molecules with weights in cluster {cluster_id}: {e}")
        raise HTTPException(
            status_code=500,
            detail=f"获取相似分子失败: {str(e)}"
        )
