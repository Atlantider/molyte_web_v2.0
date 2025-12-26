"""
Reaction Network API Routes
反应网络API路由

现代化REST API设计:
- 清晰的资源命名
- 完整的CRUD操作
- 智能的错误处理
- 统一的响应格式
"""

from fastapi import APIRouter, Depends, HTTPException, Query, status
from sqlalchemy.orm import Session
from sqlalchemy import desc, asc, or_, and_
from typing import List, Optional
import logging

from app.database import get_db
from app.dependencies import get_current_active_user
from app.models.user import User
from app.models.reaction_network import (
    ReactionNetworkJob,
    ReactionNetworkMolecule,
    ReactionNetworkReaction,
    ReactionNetworkJobStatus
)
from app.schemas.reaction_network import (
    ReactionNetworkJobCreate,
    ReactionNetworkJobUpdate,
    ReactionNetworkJobSchema,
    ReactionNetworkJobDetailSchema,
    ReactionNetworkJobListResponse,
    MoleculeSchema,
    MoleculeListResponse,
    ReactionSchema,
    ReactionListResponse,
    NetworkVisualizationData,
    NetworkNode,
    NetworkEdge,
    JobSubmissionResponse,
    ReactionNetworkJobQuery,
    MoleculeQuery,
    ReactionQuery
)

logger = logging.getLogger(__name__)
router = APIRouter()


# ============================================================================
# Job CRUD Operations
# ============================================================================

@router.post("/jobs", response_model=ReactionNetworkJobSchema)
def create_reaction_network_job(
    job_data: ReactionNetworkJobCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """创建反应网络生成任务"""
    logger.info(f"User {current_user.id} creating reaction network job: {job_data.job_name}")
    
    # 创建任务记录
    db_job = ReactionNetworkJob(
        user_id=current_user.id,
        job_name=job_data.job_name,
        description=job_data.description,
        initial_smiles=job_data.initial_smiles,
        temperature=job_data.temperature,
        electrode_type=job_data.electrode_type,
        voltage=job_data.voltage,
        max_generations=job_data.max_generations,
        max_species=job_data.max_species,
        energy_cutoff=job_data.energy_cutoff,
        slurm_partition=job_data.slurm_partition or "cpu",
        slurm_cpus=job_data.slurm_cpus or 16,
        slurm_time=job_data.slurm_time or 7200,
        config=job_data.config or {},
        status=ReactionNetworkJobStatus.CREATED
    )
    
    db.add(db_job)
    db.commit()
    db.refresh(db_job)
    
    logger.info(f"Created reaction network job {db_job.id}")
    return db_job


@router.get("/jobs", response_model=ReactionNetworkJobListResponse)
def list_reaction_network_jobs(
    status: Optional[str] = Query(None, description="按状态筛选"),
    skip: int = Query(0, ge=0),
    limit: int = Query(20, ge=1, le=100),
    sort_by: str = Query("created_at", description="排序字段"),
    sort_desc: bool = Query(True, description="降序"),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """列出反应网络任务"""
    
    # 基础查询
    query = db.query(ReactionNetworkJob).filter(
        ReactionNetworkJob.user_id == current_user.id,
        ReactionNetworkJob.is_deleted == False
    )
    
    # 状态筛选
    if status:
        query = query.filter(ReactionNetworkJob.status == status)
    
    # 统计总数
    total = query.count()
    
    # 排序
    order_col = getattr(ReactionNetworkJob, sort_by, ReactionNetworkJob.created_at)
    if sort_desc:
        query = query.order_by(desc(order_col))
    else:
        query = query.order_by(asc(order_col))
    
    # 分页
    jobs = query.offset(skip).limit(limit).all()
    
    return ReactionNetworkJobListResponse(total=total, jobs=jobs)


@router.get("/jobs/{job_id}", response_model=ReactionNetworkJobDetailSchema)
def get_reaction_network_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取反应网络任务详情"""
    
    job = db.query(ReactionNetworkJob).filter(
        ReactionNetworkJob.id == job_id,
        ReactionNetworkJob.user_id == current_user.id,
        ReactionNetworkJob.is_deleted == False
    ).first()
    
    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Job {job_id} not found"
        )
    
    return job


@router.post("/jobs/{job_id}/submit", response_model=JobSubmissionResponse)
def submit_reaction_network_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """提交任务到计算集群"""
    
    job = db.query(ReactionNetworkJob).filter(
        ReactionNetworkJob.id == job_id,
        ReactionNetworkJob.user_id == current_user.id,
        ReactionNetworkJob.is_deleted == False
    ).first()
    
    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Job {job_id} not found"
        )
    
    if job.status != ReactionNetworkJobStatus.CREATED:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Job {job_id} cannot be submitted (current status: {job.status})"
        )
    
    # 更新状态为QUEUED
    job.status = ReactionNetworkJobStatus.QUEUED
    db.commit()
    
    logger.info(f"Job {job_id} submitted to queue")
    
    return JobSubmissionResponse(
        success=True,
        message="Job submitted successfully",
        job_id=job.id,
        slurm_job_id=None  # Will be set by worker
    )


@router.delete("/jobs/{job_id}")
def delete_reaction_network_job(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """删除反应网络任务（软删除）"""
    
    job = db.query(ReactionNetworkJob).filter(
        ReactionNetworkJob.id == job_id,
        ReactionNetworkJob.user_id == current_user.id,
        ReactionNetworkJob.is_deleted == False
    ).first()
    
    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Job {job_id} not found"
        )
    
    # 软删除
    job.is_deleted = True
    from datetime import datetime
    job.deleted_at = datetime.utcnow()
    job.deleted_by = current_user.id
    
    db.commit()
    
    logger.info(f"Job {job_id} deleted by user {current_user.id}")
    
    return {"message": "Job deleted successfully", "job_id": job.id}


# ============================================================================
# Molecules
# ============================================================================

@router.get("/jobs/{job_id}/molecules", response_model=MoleculeListResponse)
def get_job_molecules(
    job_id: int,
    generation: Optional[int] = Query(None, ge=0),
    min_energy: Optional[float] = None,
    max_energy: Optional[float] = None,
    skip: int = Query(0, ge=0),
    limit: int = Query(50, ge=1, le=200),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取任务的分子列表"""
    
    # 验证任务权限
    job = db.query(ReactionNetworkJob).filter(
        ReactionNetworkJob.id == job_id,
        ReactionNetworkJob.user_id == current_user.id
    ).first()
    
    if not job:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Job not found")
    
    # 查询分子
    query = db.query(ReactionNetworkMolecule).filter(
        ReactionNetworkMolecule.job_id == job_id
    )
    
    if generation is not None:
        query = query.filter(ReactionNetworkMolecule.generation == generation)
    
    if min_energy is not None:
        query = query.filter(ReactionNetworkMolecule.energy_kcal >= min_energy)
    
    if max_energy is not None:
        query = query.filter(ReactionNetworkMolecule.energy_kcal <= max_energy)
    
    total = query.count()
    molecules = query.order_by(ReactionNetworkMolecule.generation, ReactionNetworkMolecule.id)\
                     .offset(skip).limit(limit).all()
    
    return MoleculeListResponse(total=total, molecules=molecules)


# ============================================================================
# Reactions
# ============================================================================

@router.get("/jobs/{job_id}/reactions", response_model=ReactionListResponse)
def get_job_reactions(
    job_id: int,
    operator: Optional[str] = None,
    reaction_type: Optional[str] = None,
    min_energy: Optional[float] = None,
    max_energy: Optional[float] = None,
    skip: int = Query(0, ge=0),
    limit: int = Query(50, ge=1, le=200),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取任务的反应列表"""
    
    # 验证任务权限
    job = db.query(ReactionNetworkJob).filter(
        ReactionNetworkJob.id == job_id,
        ReactionNetworkJob.user_id == current_user.id
    ).first()
    
    if not job:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Job not found")
    
    # 查询反应
    query = db.query(ReactionNetworkReaction).filter(
        ReactionNetworkReaction.job_id == job_id
    )
    
    if operator:
        query = query.filter(ReactionNetworkReaction.operator_name == operator)
    
    if reaction_type:
        query = query.filter(ReactionNetworkReaction.reaction_type == reaction_type)
    
    if min_energy is not None:
        query = query.filter(ReactionNetworkReaction.reaction_energy >= min_energy)
    
    if max_energy is not None:
        query = query.filter(ReactionNetworkReaction.reaction_energy <= max_energy)
    
    total = query.count()
    reactions = query.order_by(ReactionNetworkReaction.id)\
                     .offset(skip).limit(limit).all()
    
    return ReactionListResponse(total=total, reactions=reactions)


# ============================================================================
# Network Visualization
# ============================================================================

@router.get("/jobs/{job_id}/network", response_model=NetworkVisualizationData)
def get_network_visualization_data(
    job_id: int,
    max_generation: Optional[int] = Query(None, ge=0, description="限制最大代数"),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取网络可视化数据"""
    
    # 验证任务权限
    job = db.query(ReactionNetworkJob).filter(
        ReactionNetworkJob.id == job_id,
        ReactionNetworkJob.user_id == current_user.id
    ).first()
    
    if not job:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Job not found")
    
    # 获取分子（节点）
    mol_query = db.query(ReactionNetworkMolecule).filter(
        ReactionNetworkMolecule.job_id == job_id
    )
    
    if max_generation is not None:
        mol_query = mol_query.filter(ReactionNetworkMolecule.generation <= max_generation)
    
    molecules = mol_query.all()
    
    # 构建节点
    nodes = []
    for mol in molecules:
        nodes.append(NetworkNode(
            id=mol.smiles,
            label=mol.name,
            generation=mol.generation,
            energy=mol.energy_kcal,
            properties={
                "molecular_weight": mol.molecular_weight,
                "num_atoms": mol.num_atoms,
                "formal_charge": mol.formal_charge,
                "num_rings": mol.num_rings
            }
        ))
    
    # 获取反应（边）
    reactions = db.query(ReactionNetworkReaction).filter(
        ReactionNetworkReaction.job_id == job_id
    ).all()
    
    # 构建边
    edges = []
    for rxn in reactions:
        # 简化：假设单一反应物->单一产物
        if rxn.reactant_smiles and rxn.product_smiles:
            for reactant in rxn.reactant_smiles:
                for product in rxn.product_smiles:
                    edges.append(NetworkEdge(
                        source=reactant,
                        target=product,
                        label=rxn.operator_name,
                        operator=rxn.operator_name,
                        energy=rxn.reaction_energy,
                        properties={
                            "reaction_type": rxn.reaction_type,
                            "activation_energy": rxn.activation_energy,
                            "is_reversible": rxn.is_reversible
                        }
                    ))
    
    # 统计信息
    statistics = {
        "num_nodes": len(nodes),
        "num_edges": len(edges),
        "num_molecules": job.num_molecules,
        "num_reactions": job.num_reactions,
        "max_generation": job.max_generation_reached,
        "generation_distribution": {}
    }
    
    # 代数分布
    for mol in molecules:
        gen = mol.generation
        statistics["generation_distribution"][gen] = \
            statistics["generation_distribution"].get(gen, 0) + 1
    
    return NetworkVisualizationData(
        nodes=nodes,
        edges=edges,
        statistics=statistics
    )


# ============================================================================
# Operator Information
# ============================================================================

@router.get("/jobs/{job_id}/operators")
def get_activated_operators(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取激活的算符信息
    
    返回基于任务配置激活的反应算符列表，包括：
    - 算符名称
    - 描述
    - 适用条件
    - 激活原因
    """
    
    # 验证任务权限
    job = db.query(ReactionNetworkJob).filter(
        ReactionNetworkJob.id == job_id,
        ReactionNetworkJob.user_id == current_user.id
    ).first()
    
    if not job:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Job not found")
    
    try:
        # 导入rsnet库
        import sys
        import os
        rsnet_path = os.path.join(os.path.dirname(__file__), "../../../../rsnet-main/rsnet-main")
        if rsnet_path not in sys.path:
            sys.path.insert(0, rsnet_path)
        
        from rsnet.operators.registry import OPERATOR_REGISTRY
        from rsnet.operators.activation_rules import get_activation_rule
        from rsnet.environment import Environment
        
        # 构建环境对象
        env = Environment(
            temperature=job.temperature,
            voltage=job.voltage,
            electrode_type=job.electrode_type.value if hasattr(job.electrode_type, 'value') else job.electrode_type
        )
        
        # 获取所有算符信息
        all_operators = OPERATOR_REGISTRY.list_operators()
        
        # 构建返回结果
        activated_operators = []
        
        for op_name, op_info in all_operators.items():
            # 获取激活规则
            rule = get_activation_rule(op_name)
            
            if not rule:
                continue
            
            # 检查环境驱动力匹配
            env_drives = env.get_active_drives()
            required_drives = rule.get('required_drives', [])
            enhancing_drives = rule.get('enhancing_drives', [])
            
            # 判断是否激活
            is_activated = False
            activation_reasons = []
            
            if required_drives:
                matched_drives = [d for d in required_drives if env_drives.get(d, False)]
                if matched_drives:
                    is_activated = True
                    activation_reasons.extend([f"驱动力:{d}" for d in matched_drives])
            else:
                # 没有严格要求，检查增强驱动力
                matched_drives = [d for d in enhancing_drives if env_drives.get(d, False)]
                if matched_drives:
                    is_activated = True
                    activation_reasons.extend([f"增强驱动力:{d}" for d in matched_drives])
                else:
                    # 通用算符，总是激活
                    is_activated = True
                    activation_reasons.append("通用算符")
            
            if is_activated:
                # 构建适用条件
                conditions = []
                if required_drives:
                    conditions.append(f"需要:{', '.join(required_drives)}")
                if enhancing_drives:
                    conditions.append(f"增强:{', '.join(enhancing_drives)}")
                
                molecular_checks = rule.get('molecular_checks', {})
                if molecular_checks:
                    conditions.append(f"分子特征:{', '.join(molecular_checks.keys())}")
                
                activated_operators.append({
                    "name": op_name,
                    "description": op_info.get('description', '无描述'),
                    "weight": rule.get('weight', 0.5),
                    "conditions": conditions,
                    "activation_reasons": activation_reasons,
                    "required_drives": required_drives,
                    "enhancing_drives": enhancing_drives,
                    "molecular_checks": list(molecular_checks.keys()) if molecular_checks else []
                })
        
        # 按权重排序
        activated_operators.sort(key=lambda x: x['weight'], reverse=True)
        
        return {
            "job_id": job.id,
            "num_operators": len(activated_operators),
            "operators": activated_operators,
            "environment": {
                "temperature": job.temperature,
                "voltage": job.voltage,
                "electrode_type": job.electrode_type.value if hasattr(job.electrode_type, 'value') else job.electrode_type,
                "active_drives": list(env_drives.keys())
            }
        }
        
    except ImportError as e:
        logger.error(f"Failed to import rsnet: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to load operator information"
        )
    except Exception as e:
        logger.error(f"Error getting activated operators: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=str(e)
        )

