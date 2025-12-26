"""
Predict Block 集成服务

提供分子性质预测、聚类分析等功能
"""

import logging
import pandas as pd
import numpy as np
import pickle
import base64
import io
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw
from .molecule_info_service import get_molecule_info_fast

logger = logging.getLogger(__name__)

# 数据文件路径 - 使用后端数据目录
BACKEND_DATA_PATH = Path("/opt/molyte_web_v1.0/backend/data")
PREDICT_BLOCK_PATH = Path("/opt/molyte_web_v1.0/predict_block-master/predict_block-master")
DATASETS_PATH = PREDICT_BLOCK_PATH / "datasets"

# 数据集缓存
_bp_data_cache = None
_mp_data_cache = None
_fp_data_cache = None
_cluster_data_cache = None

# 指纹缓存
_fingerprint_cache = {}
_bp_fingerprints_cache = None
_mp_fingerprints_cache = None
_fp_fingerprints_cache = None
_cluster_fingerprints_cache = None


def load_bp_data() -> Optional[pd.DataFrame]:
    """加载沸点数据集"""
    global _bp_data_cache
    if _bp_data_cache is not None:
        return _bp_data_cache

    try:
        # 优先使用小型数据集以提高速度
        bp_file = BACKEND_DATA_PATH / "bp_data_small.csv"
        if bp_file.exists():
            _bp_data_cache = pd.read_csv(bp_file)
            logger.info(f"Loaded BP data (small) from backend: {len(_bp_data_cache)} samples")
            return _bp_data_cache

        # 备选：使用完整数据集
        bp_file = BACKEND_DATA_PATH / "bp_data.csv"
        if bp_file.exists():
            _bp_data_cache = pd.read_csv(bp_file)
            logger.info(f"Loaded BP data from backend: {len(_bp_data_cache)} samples")
            return _bp_data_cache

        # 备选：从 predict_block 的数据中加载
        predictions_file = PREDICT_BLOCK_PATH / "dm" / "predictions_all_properties.csv"
        if predictions_file.exists():
            df = pd.read_csv(predictions_file)
            # 过滤包含 BP 数据的行
            if 'BP' in df.columns:
                _bp_data_cache = df[['SMILES', 'BP']].dropna()
                _bp_data_cache.columns = ['SMILES', 'BoilingPoint']
                logger.info(f"Loaded BP data from predict_block: {len(_bp_data_cache)} samples")
                return _bp_data_cache
    except Exception as e:
        logger.error(f"Failed to load BP data: {e}")

    return None


def load_mp_data() -> Optional[pd.DataFrame]:
    """加载熔点数据集"""
    global _mp_data_cache
    if _mp_data_cache is not None:
        return _mp_data_cache

    try:
        # 优先使用小型数据集以提高速度
        mp_file = BACKEND_DATA_PATH / "mp_data_small.csv"
        if mp_file.exists():
            _mp_data_cache = pd.read_csv(mp_file)
            logger.info(f"Loaded MP data (small) from backend: {len(_mp_data_cache)} samples")
            return _mp_data_cache

        # 备选：使用完整数据集
        mp_file = BACKEND_DATA_PATH / "mp_data.csv"
        if mp_file.exists():
            _mp_data_cache = pd.read_csv(mp_file)
            logger.info(f"Loaded MP data from backend: {len(_mp_data_cache)} samples")
            return _mp_data_cache

        # 备选：从 predict_block 的数据中加载
        mp_file = DATASETS_PATH / "bp_mp_fp" / "MP_clean.csv"
        if mp_file.exists():
            _mp_data_cache = pd.read_csv(mp_file)
            logger.info(f"Loaded MP data from predict_block: {len(_mp_data_cache)} samples")
            return _mp_data_cache
    except Exception as e:
        logger.error(f"Failed to load MP data: {e}")

    return None


def load_fp_data() -> Optional[pd.DataFrame]:
    """加载闪点数据集"""
    global _fp_data_cache
    if _fp_data_cache is not None:
        return _fp_data_cache

    try:
        # 优先使用小型数据集以提高速度
        fp_file = BACKEND_DATA_PATH / "fp_data_small.csv"
        if fp_file.exists():
            _fp_data_cache = pd.read_csv(fp_file)
            logger.info(f"Loaded FP data (small) from backend: {len(_fp_data_cache)} samples")
            return _fp_data_cache

        # 备选：使用完整数据集
        fp_file = BACKEND_DATA_PATH / "fp_data.csv"
        if fp_file.exists():
            _fp_data_cache = pd.read_csv(fp_file)
            logger.info(f"Loaded FP data from backend: {len(_fp_data_cache)} samples")
            return _fp_data_cache

        # 备选：从 predict_block 的数据中加载
        fp_file = DATASETS_PATH / "bp_mp_fp" / "FP_clean.csv"
        if fp_file.exists():
            _fp_data_cache = pd.read_csv(fp_file)
            logger.info(f"Loaded FP data from predict_block: {len(_fp_data_cache)} samples")
            return _fp_data_cache
    except Exception as e:
        logger.error(f"Failed to load FP data: {e}")

    return None


def load_cluster_data() -> Optional[pd.DataFrame]:
    """加载聚类数据"""
    global _cluster_data_cache
    if _cluster_data_cache is not None:
        return _cluster_data_cache

    try:
        # 优先使用小型数据集以提高速度
        cluster_file = BACKEND_DATA_PATH / "cluster_data_small.csv"
        if cluster_file.exists():
            _cluster_data_cache = pd.read_csv(cluster_file)
            logger.info(f"Loaded cluster data (small) from backend: {len(_cluster_data_cache)} samples")
            return _cluster_data_cache

        # 备选：使用完整数据集
        cluster_file = BACKEND_DATA_PATH / "cluster_data.csv"
        if cluster_file.exists():
            _cluster_data_cache = pd.read_csv(cluster_file)
            logger.info(f"Loaded cluster data from backend: {len(_cluster_data_cache)} samples")
            return _cluster_data_cache

        # 备选：从 predict_block 的数据中加载
        cluster_file = DATASETS_PATH / "QM9" / "raw" / "gdb9_clusters_with_expr_1+2+5_50lei.csv"
        if cluster_file.exists():
            _cluster_data_cache = pd.read_csv(cluster_file)
            logger.info(f"Loaded cluster data from predict_block: {len(_cluster_data_cache)} samples")
            return _cluster_data_cache
    except Exception as e:
        logger.error(f"Failed to load cluster data: {e}")

    return None


def load_precomputed_fingerprints(property_name: str) -> Optional[Dict]:
    """加载预计算的指纹"""
    fingerprint_files = {
        'bp': BACKEND_DATA_PATH / "bp_fingerprints_small.pkl",
        'mp': BACKEND_DATA_PATH / "mp_fingerprints_small.pkl",
        'fp': BACKEND_DATA_PATH / "fp_fingerprints_small.pkl",
        'cluster': BACKEND_DATA_PATH / "cluster_fingerprints_small.pkl"
    }

    file_path = fingerprint_files.get(property_name)
    if not file_path or not file_path.exists():
        return None

    try:
        with open(file_path, 'rb') as f:
            fingerprint_data = pickle.load(f)
        logger.info(f"Loaded precomputed fingerprints for {property_name}: {len(fingerprint_data['fingerprints'])} samples")
        return fingerprint_data
    except Exception as e:
        logger.error(f"Failed to load precomputed fingerprints for {property_name}: {e}")
        return None


def check_if_real_value(smiles: str, property_name: str) -> Tuple[bool, Optional[float]]:
    """
    检查分子的属性值是否为真实值（存在于数据集中）

    Args:
        smiles: 分子的 SMILES 表示
        property_name: 属性名称 ('bp', 'mp', 'fp', 'alpha', 'mu', 'gap', 'homo', 'lumo')

    Returns:
        (是否为真实值, 真实值或None)
    """
    # 处理QM9属性
    if property_name in ["alpha", "mu", "gap", "homo", "lumo"]:
        cluster_data = load_cluster_data()
        if cluster_data is None:
            return False, None

        matching_rows = cluster_data[cluster_data['SMILES'] == smiles]
        if not matching_rows.empty:
            real_value = matching_rows.iloc[0][property_name]
            if pd.notna(real_value):
                return True, float(real_value)
        return False, None

    # 处理BP/MP/FP属性
    if property_name == "bp":
        data = load_bp_data()
        col_name = "BoilingPoint"
    elif property_name == "mp":
        data = load_mp_data()
        col_name = "MeltingPoint"
    elif property_name == "fp":
        data = load_fp_data()
        col_name = "FlashPoint"
    else:
        return False, None

    if data is None:
        return False, None

    # 查找完全匹配的 SMILES
    matching_rows = data[data['SMILES'] == smiles]
    if not matching_rows.empty:
        real_value = matching_rows.iloc[0][col_name]
        if pd.notna(real_value):
            return True, float(real_value)

    return False, None


def calculate_molecular_fingerprint(smiles: str) -> Optional[np.ndarray]:
    """计算分子指纹"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # 使用Morgan指纹
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        return np.array(fp)
    except Exception as e:
        logger.error(f"Failed to calculate fingerprint for {smiles}: {e}")
        return None


def predict_property_detailed(smiles: str, property_name: str) -> Optional[Dict]:
    """
    预测分子性质（详细版本）

    返回包含预测值、是否为真实值、分子信息等的详细字典
    """
    # 检查是否为真实值
    is_real, real_value = check_if_real_value(smiles, property_name)

    if is_real:
        # 如果是真实值，直接返回
        predicted_value = real_value
    else:
        # 否则进行预测
        predicted_value = predict_property(smiles, property_name)

    if predicted_value is None:
        return None

    # 获取分子信息
    molecule_info = get_molecule_info_fast(smiles)

    return {
        'value': predicted_value,
        'is_real': is_real,
        'property_name': property_name,
        'molecule_info': molecule_info
    }


def predict_property(smiles: str, property_name: str) -> Optional[float]:
    """
    预测分子性质

    支持8种属性：bp, mp, fp, alpha, mu, gap, homo, lumo
    使用最近邻算法基于相似分子的已知性质进行预测
    """
    global _bp_fingerprints_cache, _mp_fingerprints_cache, _fp_fingerprints_cache

    # 处理QM9属性
    if property_name in ["alpha", "mu", "gap", "homo", "lumo"]:
        return predict_qm9_property(smiles, property_name)

    # 处理BP/MP/FP属性
    if property_name == "bp":
        data = load_bp_data()
        col_name = "BoilingPoint"
        cache_key = 'bp'
    elif property_name == "mp":
        data = load_mp_data()
        col_name = "MeltingPoint"
        cache_key = 'mp'
    elif property_name == "fp":
        data = load_fp_data()
        col_name = "FlashPoint"
        cache_key = 'fp'
    else:
        return None

    if data is None or data.empty:
        return None

    try:
        # 计算输入分子的指纹
        query_fp = calculate_molecular_fingerprint(smiles)
        if query_fp is None:
            return None

        # 加载预计算的指纹（只做一次）
        fingerprint_data = None
        if cache_key == 'bp' and _bp_fingerprints_cache is None:
            fingerprint_data = load_precomputed_fingerprints('bp')
            if fingerprint_data:
                _bp_fingerprints_cache = fingerprint_data
                logger.info(f"Loaded {len(fingerprint_data['fingerprints'])} precomputed BP fingerprints")
            else:
                logger.warning("No precomputed BP fingerprints found, using fallback")
                return None
        elif cache_key == 'mp' and _mp_fingerprints_cache is None:
            fingerprint_data = load_precomputed_fingerprints('mp')
            if fingerprint_data:
                _mp_fingerprints_cache = fingerprint_data
                logger.info(f"Loaded {len(fingerprint_data['fingerprints'])} precomputed MP fingerprints")
            else:
                logger.warning("No precomputed MP fingerprints found, using fallback")
                return None
        elif cache_key == 'fp' and _fp_fingerprints_cache is None:
            fingerprint_data = load_precomputed_fingerprints('fp')
            if fingerprint_data:
                _fp_fingerprints_cache = fingerprint_data
                logger.info(f"Loaded {len(fingerprint_data['fingerprints'])} precomputed FP fingerprints")
            else:
                logger.warning("No precomputed FP fingerprints found, using fallback")
                return None

        # 获取缓存的指纹数据
        if cache_key == 'bp':
            fingerprint_cache = _bp_fingerprints_cache
        elif cache_key == 'mp':
            fingerprint_cache = _mp_fingerprints_cache
        else:
            fingerprint_cache = _fp_fingerprints_cache

        if fingerprint_cache is None:
            logger.error(f"No fingerprint cache available for {cache_key}")
            return None

        # 使用预计算的指纹数据
        if isinstance(fingerprint_cache, dict) and 'fingerprints' in fingerprint_cache:
            fingerprints = fingerprint_cache['fingerprints']
            valid_indices = fingerprint_cache['valid_indices']
        else:
            logger.error(f"Invalid fingerprint cache format for {cache_key}")
            return None

        if len(fingerprints) == 0:
            return None

        # 使用KNN找到最相似的分子
        knn = NearestNeighbors(n_neighbors=min(5, len(fingerprints)), metric='euclidean')
        knn.fit(fingerprints)

        distances, indices = knn.kneighbors([query_fp])

        # 使用距离加权平均预测
        neighbor_indices = [valid_indices[i] for i in indices[0]]

        # 获取对应的属性值，处理可能的索引问题
        neighbor_values = []
        for idx in neighbor_indices:
            if idx < len(data):
                value = data.iloc[idx][col_name]
                if pd.notna(value):
                    neighbor_values.append(float(value))

        if not neighbor_values:
            logger.warning(f"No valid neighbor values found for {smiles}")
            return None

        neighbor_distances = distances[0][:len(neighbor_values)]

        # 距离越小，权重越大
        weights = 1.0 / (neighbor_distances + 1e-6)
        weights /= weights.sum()

        predicted_value = np.sum(neighbor_values * weights)
        return float(predicted_value)

    except Exception as e:
        logger.error(f"Failed to predict {property_name} for {smiles}: {e}")
        return None


def predict_qm9_property(smiles: str, property_name: str) -> Optional[float]:
    """
    预测单个QM9属性

    基于相似分子的属性值进行加权平均预测
    """
    try:
        cluster_data = load_cluster_data()
        if cluster_data is None or cluster_data.empty:
            return None

        # 查找查询分子
        query_matches = cluster_data[cluster_data['SMILES'] == smiles]
        if not query_matches.empty:
            # 如果分子在数据集中，直接返回其属性值
            query_row = query_matches.iloc[0]
            if property_name in query_row and pd.notna(query_row[property_name]):
                return float(query_row[property_name])

        # 如果分子不在数据集中，基于相似分子预测
        # 使用QM9属性距离查找最相似的5个分子
        qm9_features = ['mu', 'alpha', 'homo', 'lumo', 'gap']

        # 计算查询分子的指纹
        query_fp = calculate_molecular_fingerprint(smiles)
        if query_fp is None:
            # 如果无法计算指纹，使用全局平均值
            values = cluster_data[property_name].dropna()
            if len(values) > 0:
                return float(values.mean())
            return None

        # 加载预计算的指纹
        fingerprint_data = load_precomputed_fingerprints('bp')
        if not fingerprint_data:
            # 如果没有指纹数据，使用全局平均值
            values = cluster_data[property_name].dropna()
            if len(values) > 0:
                return float(values.mean())
            return None

        fingerprints = fingerprint_data['fingerprints']
        valid_indices = fingerprint_data['valid_indices']

        # 使用KNN查找最相似的分子
        knn = NearestNeighbors(n_neighbors=min(5, len(fingerprints)), metric='euclidean')
        knn.fit(fingerprints)
        distances, indices = knn.kneighbors([query_fp])

        # 获取相似分子的属性值
        neighbor_values = []
        neighbor_distances = []
        for i, idx in enumerate([valid_indices[j] for j in indices[0]]):
            if idx < len(cluster_data):
                row = cluster_data.iloc[idx]
                if property_name in row and pd.notna(row[property_name]):
                    neighbor_values.append(float(row[property_name]))
                    neighbor_distances.append(distances[0][i])

        if not neighbor_values:
            # 如果没有找到相似分子的属性值，使用全局平均值
            values = cluster_data[property_name].dropna()
            if len(values) > 0:
                return float(values.mean())
            return None

        # 距离加权平均
        weights = 1.0 / (np.array(neighbor_distances) + 1e-6)
        weights /= weights.sum()
        predicted_value = np.sum(np.array(neighbor_values) * weights)
        return float(predicted_value)

    except Exception as e:
        logger.error(f"Failed to predict QM9 property {property_name} for {smiles}: {e}")
        return None


def predict_qm9_properties(smiles: str) -> Optional[Dict]:
    """
    预测 QM9 量子化学性质

    返回 mu, alpha, homo, lumo, gap 等属性
    """
    try:
        qm9_features = ['mu', 'alpha', 'homo', 'lumo', 'gap']
        predicted_props = {}

        for feature in qm9_features:
            value = predict_qm9_property(smiles, feature)
            if value is not None:
                predicted_props[feature] = value

        return predicted_props if predicted_props else None

    except Exception as e:
        logger.error(f"Failed to predict QM9 properties for {smiles}: {e}")
        return None


def get_qm9_properties(mol_id: str) -> Optional[Dict]:
    """
    从 QM9 聚类数据中获取分子的量子化学性质

    返回 mu, alpha, homo, lumo, gap 等属性
    """
    try:
        data = load_cluster_data()
        if data is None or data.empty:
            return None

        # 查找匹配的分子
        mol_data = data[data['mol_expr'] == mol_id]
        if mol_data.empty:
            return None

        row = mol_data.iloc[0]
        properties = {}

        # 提取 QM9 属性
        qm9_features = ['mu', 'alpha', 'homo', 'lumo', 'gap', 'r2', 'zpve', 'u0', 'u298', 'h298', 'g298', 'cv']
        for feature in qm9_features:
            if feature in row and pd.notna(row[feature]):
                properties[feature] = float(row[feature])

        return properties if properties else None

    except Exception as e:
        logger.error(f"Failed to get QM9 properties for {mol_id}: {e}")
        return None


def find_similar_molecules(smiles: str, n_similar: int = 50) -> List[Dict]:
    """
    查找相似分子（使用预计算指纹）

    基于分子指纹相似性查找数据集中最相似的分子
    """
    global _bp_fingerprints_cache

    try:
        # 计算查询分子的指纹
        query_fp = calculate_molecular_fingerprint(smiles)
        if query_fp is None:
            return []

        # 加载预计算的指纹
        if _bp_fingerprints_cache is None:
            fingerprint_data = load_precomputed_fingerprints('bp')
            if fingerprint_data:
                _bp_fingerprints_cache = fingerprint_data
                logger.info(f"Loaded {len(fingerprint_data['fingerprints'])} precomputed fingerprints for similarity search")
            else:
                logger.warning("No precomputed fingerprints available for similarity search")
                return []

        fingerprints = _bp_fingerprints_cache['fingerprints']
        valid_indices = _bp_fingerprints_cache['valid_indices']

        # 加载数据
        bp_data = load_bp_data()
        if bp_data is None:
            return []

        # 使用KNN查找相似分子
        knn = NearestNeighbors(n_neighbors=min(n_similar, len(fingerprints)), metric='euclidean')
        knn.fit(fingerprints)

        distances, indices = knn.kneighbors([query_fp])

        results = []
        for i, idx in enumerate([valid_indices[j] for j in indices[0]]):
            if idx >= len(bp_data):
                continue

            row = bp_data.iloc[idx]
            properties = {'bp': float(row['BoilingPoint'])}

            # 尝试获取其他属性
            try:
                mp_data = load_mp_data()
                if mp_data is not None:
                    mp_match = mp_data[mp_data['SMILES'] == row['SMILES']]
                    if not mp_match.empty:
                        properties['mp'] = float(mp_match.iloc[0]['MeltingPoint'])
            except:
                pass

            try:
                fp_data = load_fp_data()
                if fp_data is not None:
                    fp_match = fp_data[fp_data['SMILES'] == row['SMILES']]
                    if not fp_match.empty:
                        properties['fp'] = float(fp_match.iloc[0]['FlashPoint'])
            except:
                pass

            results.append({
                'smiles': row['SMILES'],
                'properties': properties,
                'similarity': float(1.0 / (distances[0][i] + 1e-6))
            })

        return results

    except Exception as e:
        logger.error(f"Failed to find similar molecules for {smiles}: {e}")
        return []


def find_similar_molecules_by_qm9_distance(
    smiles: str,
    n_similar: int = 50,
    scope: str = "global"
) -> List[Dict]:
    """
    基于QM9属性距离查找相似分子

    使用alpha, mu, homo, lumo, gap属性计算欧氏距离
    支持两种scope：
    - 'cluster': 在同一聚类内查找相似分子
    - 'global': 在全局范围内查找相似分子
    """
    try:
        cluster_data = load_cluster_data()
        if cluster_data is None or cluster_data.empty:
            logger.warning("Cluster data not available for QM9 distance search")
            return []

        # QM9属性列表
        qm9_features = ['mu', 'alpha', 'homo', 'lumo', 'gap']

        # 查找查询分子在数据集中的索引
        query_matches = cluster_data[cluster_data['SMILES'] == smiles]

        if query_matches.empty:
            # 如果查询分子不在数据库中，使用预测的属性
            logger.info(f"Query SMILES {smiles} not found in cluster data, using predicted properties")

            # 预测查询分子的QM9属性
            query_props = []
            for feature in qm9_features:
                predicted_value = predict_qm9_property(smiles, feature)
                if predicted_value is not None:
                    query_props.append(float(predicted_value))
                else:
                    query_props.append(0.0)
            query_props = np.array(query_props)

            # 如果预测失败，返回空列表
            if np.all(query_props == 0.0):
                logger.warning(f"Failed to predict properties for {smiles}")
                return []
        else:
            # 使用数据库中的属性
            query_idx = query_matches.index[0]
            query_row = cluster_data.iloc[query_idx]

            # 获取查询分子的QM9属性向量
            query_props = []
            for feature in qm9_features:
                if feature in query_row and pd.notna(query_row[feature]):
                    query_props.append(float(query_row[feature]))
                else:
                    query_props.append(0.0)
            query_props = np.array(query_props)

        # 确定候选分子范围
        if scope == "cluster" and not query_matches.empty:
            # 在同一聚类内查找（仅当查询分子在数据库中时）
            query_cluster = query_row.get('cluster')
            if pd.isna(query_cluster):
                candidates = cluster_data
            else:
                candidates = cluster_data[cluster_data['cluster'] == query_cluster]
        else:
            # 全局范围
            candidates = cluster_data

        # 计算所有候选分子与查询分子的距离
        distances = []
        for idx, (_, row) in enumerate(candidates.iterrows()):
            # 如果查询分子在数据库中，跳过它本身
            if not query_matches.empty and idx == query_idx:
                continue

            mol_props = []
            for feature in qm9_features:
                if feature in row and pd.notna(row[feature]):
                    mol_props.append(float(row[feature]))
                else:
                    mol_props.append(0.0)
            mol_props = np.array(mol_props)

            # 计算欧氏距离
            dist = np.linalg.norm(query_props - mol_props)
            distances.append((dist, idx, row))

        # 按距离排序
        distances.sort(key=lambda x: x[0])

        # 返回最相似的n_similar个分子
        results = []
        all_features = ['alpha', 'mu', 'gap', 'homo', 'lumo', 'BP', 'FP', 'MP']

        for dist, idx, row in distances[:n_similar]:
            props = {}
            # 包含所有8个属性
            for feature in all_features:
                col_name = feature.lower() if feature in ['bp', 'fp', 'mp'] else feature
                if col_name in row and pd.notna(row[col_name]):
                    props[feature.lower()] = float(row[col_name])

            results.append({
                'smiles': row['SMILES'],
                'properties': props,
                'distance': float(dist),
                'similarity': float(1.0 / (dist + 1e-6))  # 距离越小，相似度越高
            })

        return results

    except Exception as e:
        logger.error(f"Failed to find similar molecules by QM9 distance for {smiles}: {e}")
        return []


def get_cluster_statistics() -> Optional[Dict]:
    """获取聚类统计信息"""
    try:
        data = load_cluster_data()
        if data is None or data.empty:
            return None

        cluster_sizes = data['cluster'].value_counts()

        return {
            'total_clusters': int(data['cluster'].max()),
            'total_molecules': len(data),
            'max_cluster_size': int(cluster_sizes.max()),
            'min_cluster_size': int(cluster_sizes.min()),
            'avg_cluster_size': float(cluster_sizes.mean()),
            'cluster_sizes': cluster_sizes.to_dict()
        }
    except Exception as e:
        logger.error(f"Failed to get cluster statistics: {e}")
        return None


def get_cluster_molecules(cluster_id: int, limit: int = 100) -> Optional[List[Dict]]:
    """
    获取聚类中的分子列表

    返回聚类中的分子及其属性信息
    """
    try:
        data = load_cluster_data()
        if data is None or data.empty:
            return None

        cluster_data = data[data['cluster'] == cluster_id]
        if cluster_data.empty:
            return None

        molecules = []
        features = ["mu", "alpha", "homo", "lumo", "gap"]

        # 限制返回的分子数量
        for idx, (_, row) in enumerate(cluster_data.head(limit).iterrows()):
            mol_id = row.get('mol_expr', f'mol_{idx}')

            # 提取属性
            properties = {}
            for feature in features:
                if feature in row and pd.notna(row[feature]):
                    properties[feature] = float(row[feature])

            molecules.append({
                'mol_id': str(mol_id),
                'properties': properties
            })

        return molecules

    except Exception as e:
        logger.error(f"Failed to get cluster molecules for cluster {cluster_id}: {e}")
        return None


def get_similar_molecules_in_cluster(cluster_id: int, center_mol_index: int, n_similar: int = 10) -> Optional[Dict]:
    """
    获取聚类中与指定分子相似的分子

    基于 QM9 属性的欧几里得距离计算相似性
    """
    try:
        data = load_cluster_data()
        if data is None or data.empty:
            return None

        cluster_data = data[data['cluster'] == cluster_id].reset_index(drop=True)
        if cluster_data.empty or center_mol_index >= len(cluster_data):
            return None

        # 获取中心分子
        center_row = cluster_data.iloc[center_mol_index]
        qm9_features = ['mu', 'alpha', 'homo', 'lumo', 'gap']

        # 获取中心分子的属性向量
        center_props = []
        for feature in qm9_features:
            if feature in center_row and pd.notna(center_row[feature]):
                center_props.append(float(center_row[feature]))
            else:
                center_props.append(0.0)
        center_props = np.array(center_props)

        # 计算所有分子与中心分子的距离
        distances = []
        for idx, (_, row) in enumerate(cluster_data.iterrows()):
            mol_props = []
            for feature in qm9_features:
                if feature in row and pd.notna(row[feature]):
                    mol_props.append(float(row[feature]))
                else:
                    mol_props.append(0.0)
            mol_props = np.array(mol_props)

            # 计算欧几里得距离
            dist = np.linalg.norm(center_props - mol_props)
            distances.append((idx, dist))

        # 按距离排序（最相似的在前）
        distances.sort(key=lambda x: x[1])

        # 获取最相似的分子（包括中心分子本身）
        similar_molecules = []
        for idx, dist in distances[:n_similar + 1]:  # +1 包括中心分子
            row = cluster_data.iloc[idx]
            properties = {}
            for feature in qm9_features:
                if feature in row and pd.notna(row[feature]):
                    properties[feature] = float(row[feature])

            similar_molecules.append({
                'index': int(idx),
                'properties': properties,
                'distance': float(dist),
                'is_center': idx == center_mol_index
            })

        return {
            'cluster_id': cluster_id,
            'center_index': center_mol_index,
            'center_properties': {
                feature: float(center_row[feature])
                for feature in qm9_features
                if feature in center_row and pd.notna(center_row[feature])
            },
            'similar_molecules': similar_molecules
        }

    except Exception as e:
        logger.error(f"Failed to get similar molecules in cluster {cluster_id}: {e}")
        return None


def get_cluster_center_molecule(cluster_id: int) -> Optional[Dict]:
    """
    获取聚类的中心分子（属性均值）

    计算聚类中所有分子属性的平均值作为中心分子
    """
    try:
        data = load_cluster_data()
        if data is None or data.empty:
            return None

        cluster_data = data[data['cluster'] == cluster_id]
        if cluster_data.empty:
            return None

        features = ["mu", "alpha", "homo", "lumo", "gap"]
        center_properties = {}

        for feature in features:
            if feature in cluster_data.columns:
                values = cluster_data[feature].dropna()
                if len(values) > 0:
                    center_properties[feature] = float(values.mean())

        return {
            'cluster_id': cluster_id,
            'size': len(cluster_data),
            'center_properties': center_properties
        }

    except Exception as e:
        logger.error(f"Failed to get cluster center molecule for cluster {cluster_id}: {e}")
        return None


def get_similar_molecules_with_weights(cluster_id: int, center_mol_index: int, n_similar: int = 10, weights: Optional[Dict[str, float]] = None) -> Optional[Dict]:
    """
    获取聚类中与指定分子相似的分子，支持自定义权重

    基于 QM9 属性的加权欧几里得距离计算相似性
    """
    try:
        data = load_cluster_data()
        if data is None or data.empty:
            return None

        cluster_data = data[data['cluster'] == cluster_id].reset_index(drop=True)
        if cluster_data.empty or center_mol_index >= len(cluster_data):
            return None

        # 获取中心分子
        center_row = cluster_data.iloc[center_mol_index]
        qm9_features = ['mu', 'alpha', 'homo', 'lumo', 'gap']

        # 设置默认权重
        if weights is None:
            weights = {feature: 1.0 for feature in qm9_features}
        else:
            # 确保所有特征都有权重
            for feature in qm9_features:
                if feature not in weights:
                    weights[feature] = 1.0

        # 获取中心分子的属性向量
        center_props = []
        for feature in qm9_features:
            if feature in center_row and pd.notna(center_row[feature]):
                center_props.append(float(center_row[feature]))
            else:
                center_props.append(0.0)
        center_props = np.array(center_props)

        # 计算所有分子与中心分子的加权距离
        distances = []
        for idx, (_, row) in enumerate(cluster_data.iterrows()):
            mol_props = []
            for feature in qm9_features:
                if feature in row and pd.notna(row[feature]):
                    mol_props.append(float(row[feature]))
                else:
                    mol_props.append(0.0)
            mol_props = np.array(mol_props)

            # 计算加权欧几里得距离
            weight_vector = np.array([weights.get(feature, 1.0) for feature in qm9_features])
            weighted_diff = (center_props - mol_props) * weight_vector
            dist = np.linalg.norm(weighted_diff)
            distances.append((idx, dist))

        # 按距离排序（最相似的在前）
        distances.sort(key=lambda x: x[1])

        # 获取最相似的分子（包括中心分子本身）
        similar_molecules = []
        for idx, dist in distances[:n_similar + 1]:  # +1 包括中心分子
            row = cluster_data.iloc[idx]
            properties = {}
            for feature in qm9_features:
                if feature in row and pd.notna(row[feature]):
                    properties[feature] = float(row[feature])

            similar_molecules.append({
                'index': int(idx),
                'properties': properties,
                'distance': float(dist),
                'is_center': idx == center_mol_index
            })

        return {
            'cluster_id': cluster_id,
            'center_index': center_mol_index,
            'center_properties': {
                feature: float(center_row[feature])
                for feature in qm9_features
                if feature in center_row and pd.notna(center_row[feature])
            },
            'similar_molecules': similar_molecules,
            'weights': weights
        }

    except Exception as e:
        logger.error(f"Failed to get similar molecules with weights for cluster {cluster_id}: {e}")
        return None


def generate_molecule_image(smiles: str, size: Tuple[int, int] = (200, 200), dark_mode: bool = False) -> Optional[str]:
    """
    生成分子结构图像并转换为base64编码的PNG字符串

    Args:
        smiles: 分子的SMILES表示
        size: 图像大小 (宽, 高)
        dark_mode: 是否为暗色模式

    Returns:
        base64编码的PNG字符串，如果失败返回None
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Failed to parse SMILES: {smiles}")
            return None

        # 生成2D坐标
        AllChem.Compute2DCoords(mol)

        # 设置绘图选项
        drawer = Draw.rdMolDraw2D.MolDraw2DCairo(size[0], size[1])

        # 配置绘图选项
        opts = drawer.drawOptions()
        if dark_mode:
            # 暗色模式：透明背景，白色原子和键
            opts.setBackgroundColour((0, 0, 0, 0))  # 透明背景
            opts.atomLabelFontSize = 14
            opts.bondLineWidth = 2
            # 设置原子颜色为白色或浅色
            opts.updateAtomPalette({6: (0.9, 0.9, 0.9)})  # 碳原子
        else:
            # 浅色模式：透明背景，深色原子和键
            opts.setBackgroundColour((0, 0, 0, 0))  # 透明背景
            opts.atomLabelFontSize = 14
            opts.bondLineWidth = 2

        # 绘制分子
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()

        # 获取PNG数据
        png_data = drawer.GetDrawingText()

        # 转换为base64字符串
        img_base64 = base64.b64encode(png_data).decode('utf-8')
        return img_base64

    except Exception as e:
        logger.error(f"Failed to generate molecule image for {smiles}: {e}")
        # 如果Cairo绘制失败，回退到PIL方法
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None

            AllChem.Compute2DCoords(mol)

            # 使用PIL方法生成图像
            img = Draw.MolToImage(mol, size=size)

            # 如果需要透明背景，转换图像
            if img.mode != 'RGBA':
                img = img.convert('RGBA')

            # 将白色背景转为透明
            data = img.getdata()
            new_data = []
            for item in data:
                # 将白色或接近白色的像素设为透明
                if item[0] > 240 and item[1] > 240 and item[2] > 240:
                    new_data.append((255, 255, 255, 0))  # 透明
                else:
                    new_data.append(item)
            img.putdata(new_data)

            # 转换为PNG字节流
            img_bytes = io.BytesIO()
            img.save(img_bytes, format='PNG')
            img_bytes.seek(0)

            # 转换为base64字符串
            img_base64 = base64.b64encode(img_bytes.getvalue()).decode('utf-8')
            return img_base64

        except Exception as e2:
            logger.error(f"Fallback image generation also failed for {smiles}: {e2}")
            return None

