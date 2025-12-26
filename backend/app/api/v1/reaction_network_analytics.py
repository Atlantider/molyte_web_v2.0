"""
Reaction Network Advanced Analytics API
高级分析API - 统计、对比、导出功能
"""

from fastapi import APIRouter, Depends, HTTPException, Query, Response
from sqlalchemy.orm import Session
from sqlalchemy import func, and_, or_
from typing import List, Optional, Dict, Any
import json
import io
from datetime import datetime

from app.database import get_db
from app.dependencies import get_current_active_user
from app.models.user import User
from app.models.reaction_network import (
    ReactionNetworkJob,
    ReactionNetworkMolecule,
    ReactionNetworkReaction
)

router = APIRouter()


# ============================================================================
# 统计分析
# ============================================================================

@router.get("/jobs/{job_id}/statistics")
def get_job_statistics(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """获取任务的详细统计信息"""
    
    job = db.query(ReactionNetworkJob).filter(
        ReactionNetworkJob.id == job_id,
        ReactionNetworkJob.user_id == current_user.id
    ).first()
    
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    
    # 分子统计
    molecules = db.query(ReactionNetworkMolecule).filter(
        ReactionNetworkMolecule.job_id == job_id
    ).all()
    
    # 反应统计
    reactions = db.query(ReactionNetworkReaction).filter(
        ReactionNetworkReaction.job_id == job_id
    ).all()
    
    # 代数分布
    generation_dist = {}
    for mol in molecules:
        gen = mol.generation
        generation_dist[gen] = generation_dist.get(gen, 0) + 1
    
    # 能量分布
    energy_bins = {
        'very_negative': 0,  # < -50
        'negative': 0,       # -50 ~ 0
        'positive': 0,       # 0 ~ 50
        'very_positive': 0   # > 50
    }
    
    for rxn in reactions:
        if rxn.reaction_energy is None:
            continue
        if rxn.reaction_energy < -50:
            energy_bins['very_negative'] += 1
        elif rxn.reaction_energy < 0:
            energy_bins['negative'] += 1
        elif rxn.reaction_energy < 50:
            energy_bins['positive'] += 1
        else:
            energy_bins['very_positive'] += 1
    
    # 算符统计
    operator_stats = {}
    for rxn in reactions:
        op = rxn.operator_name or 'unknown'
        operator_stats[op] = operator_stats.get(op, 0) + 1
    
    # 分子量分布
    mw_bins = {
        'light': 0,     # < 100
        'medium': 0,    # 100-300
        'heavy': 0      # > 300
    }
    
    for mol in molecules:
        if mol.molecular_weight is None:
            continue
        if mol.molecular_weight < 100:
            mw_bins['light'] += 1
        elif mol.molecular_weight < 300:
            mw_bins['medium'] += 1
        else:
            mw_bins['heavy'] += 1
    
    return {
        'job_id': job_id,
        'job_name': job.job_name,
        'total_molecules': len(molecules),
        'total_reactions': len(reactions),
        'max_generation': job.max_generation_reached,
        'generation_distribution': generation_dist,
        'energy_distribution': energy_bins,
        'operator_distribution': operator_stats,
        'molecular_weight_distribution': mw_bins,
        'avg_reaction_energy': sum(r.reaction_energy for r in reactions if r.reaction_energy) / len([r for r in reactions if r.reaction_energy]) if reactions else None,
        'cpu_hours': job.actual_cpu_hours,
        'created_at': job.created_at.isoformat(),
        'finished_at': job.finished_at.isoformat() if job.finished_at else None
    }


@router.get("/jobs/{job_id}/energy-analysis")
def get_energy_analysis(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """能量分析 - 找出低能路径和关键反应"""
    
    job = db.query(ReactionNetworkJob).filter(
        ReactionNetworkJob.id == job_id,
        ReactionNetworkJob.user_id == current_user.id
    ).first()
    
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    
    # 获取所有反应
    reactions = db.query(ReactionNetworkReaction).filter(
        ReactionNetworkReaction.job_id == job_id,
        ReactionNetworkReaction.reaction_energy.isnot(None)
    ).order_by(ReactionNetworkReaction.reaction_energy).all()
    
    # 最低能反应（前10）
    lowest_energy_reactions = reactions[:10]
    
    # 最高能反应（后10）
    highest_energy_reactions = reactions[-10:]
    
    # 放能反应（负能量）
    exothermic_reactions = [r for r in reactions if r.reaction_energy < 0]
    
    # 吸能反应（正能量）
    endothermic_reactions = [r for r in reactions if r.reaction_energy > 0]
    
    return {
        'lowest_energy_reactions': [
            {
                'id': r.id,
                'reactants': r.reactant_smiles,
                'products': r.product_smiles,
                'operator': r.operator_name,
                'energy': r.reaction_energy
            }
            for r in lowest_energy_reactions
        ],
        'highest_energy_reactions': [
            {
                'id': r.id,
                'reactants': r.reactant_smiles,
                'products': r.product_smiles,
                'operator': r.operator_name,
                'energy': r.reaction_energy
            }
            for r in highest_energy_reactions
        ],
        'exothermic_count': len(exothermic_reactions),
        'endothermic_count': len(endothermic_reactions),
        'avg_exothermic_energy': sum(r.reaction_energy for r in exothermic_reactions) / len(exothermic_reactions) if exothermic_reactions else 0,
        'avg_endothermic_energy': sum(r.reaction_energy for r in endothermic_reactions) / len(endothermic_reactions) if endothermic_reactions else 0
    }


# ============================================================================
# 导出功能
# ============================================================================

@router.get("/jobs/{job_id}/export/graphml")
def export_graphml(
    job_id: int,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """导出为GraphML格式（用于Cytoscape等工具）"""
    
    job = db.query(ReactionNetworkJob).filter(
        ReactionNetworkJob.id == job_id,
        ReactionNetworkJob.user_id == current_user.id
    ).first()
    
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    
    molecules = db.query(ReactionNetworkMolecule).filter(
        ReactionNetworkMolecule.job_id == job_id
    ).all()
    
    reactions = db.query(ReactionNetworkReaction).filter(
        ReactionNetworkReaction.job_id == job_id
    ).all()
    
    # 生成GraphML
    graphml = f'''<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns">
  <key id="label" for="node" attr.name="label" attr.type="string"/>
  <key id="smiles" for="node" attr.name="smiles" attr.type="string"/>
  <key id="generation" for="node" attr.name="generation" attr.type="int"/>
  <key id="energy" for="node" attr.name="energy" attr.type="double"/>
  <key id="operator" for="edge" attr.name="operator" attr.type="string"/>
  <key id="reaction_energy" for="edge" attr.name="reaction_energy" attr.type="double"/>
  
  <graph id="G" edgedefault="directed">
'''
    
    # 添加节点
    for mol in molecules:
        graphml += f'''    <node id="{mol.smiles}">
      <data key="label">{mol.name}</data>
      <data key="smiles">{mol.smiles}</data>
      <data key="generation">{mol.generation}</data>
      <data key="energy">{mol.energy_kcal or 0}</data>
    </node>
'''
    
    # 添加边
    edge_id = 0
    for rxn in reactions:
        for reactant in rxn.reactant_smiles:
            for product in rxn.product_smiles:
                graphml += f'''    <edge id="e{edge_id}" source="{reactant}" target="{product}">
      <data key="operator">{rxn.operator_name or 'unknown'}</data>
      <data key="reaction_energy">{rxn.reaction_energy or 0}</data>
    </edge>
'''
                edge_id += 1
    
    graphml += '''  </graph>
</graphml>'''
    
    return Response(
        content=graphml,
        media_type="application/xml",
        headers={
            "Content-Disposition": f"attachment; filename=network_{job_id}.graphml"
        }
    )


@router.get("/jobs/{job_id}/export/csv")
def export_csv(
    job_id: int,
    export_type: str = Query("molecules", description="molecules or reactions"),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """导出为CSV格式"""
    
    job = db.query(ReactionNetworkJob).filter(
        ReactionNetworkJob.id == job_id,
        ReactionNetworkJob.user_id == current_user.id
    ).first()
    
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    
    if export_type == "molecules":
        molecules = db.query(ReactionNetworkMolecule).filter(
            ReactionNetworkMolecule.job_id == job_id
        ).all()
        
        csv = "id,name,smiles,generation,energy_kcal,molecular_weight,num_atoms,formal_charge\n"
        for mol in molecules:
            csv += f"{mol.id},{mol.name},{mol.smiles},{mol.generation},{mol.energy_kcal or ''},{mol.molecular_weight or ''},{mol.num_atoms or ''},{mol.formal_charge or ''}\n"
        
        filename = f"molecules_{job_id}.csv"
    
    else:  # reactions
        reactions = db.query(ReactionNetworkReaction).filter(
            ReactionNetworkReaction.job_id == job_id
        ).all()
        
        csv = "id,reactants,products,operator,reaction_type,reaction_energy,activation_energy\n"
        for rxn in reactions:
            reactants = ';'.join(rxn.reactant_smiles)
            products = ';'.join(rxn.product_smiles)
            csv += f"{rxn.id},{reactants},{products},{rxn.operator_name or ''},{rxn.reaction_type or ''},{rxn.reaction_energy or ''},{rxn.activation_energy or ''}\n"
        
        filename = f"reactions_{job_id}.csv"
    
    return Response(
        content=csv,
        media_type="text/csv",
        headers={
            "Content-Disposition": f"attachment; filename={filename}"
        }
    )


# ============================================================================
# 批量操作
# ============================================================================

@router.post("/jobs/batch/compare")
def batch_compare_jobs(
    job_ids: List[int],
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """批量对比多个任务"""
    
    jobs = db.query(ReactionNetworkJob).filter(
        ReactionNetworkJob.id.in_(job_ids),
        ReactionNetworkJob.user_id == current_user.id
    ).all()
    
    if len(jobs) != len(job_ids):
        raise HTTPException(status_code=404, detail="Some jobs not found")
    
    comparison_data = []
    
    for job in jobs:
        num_molecules = db.query(func.count(ReactionNetworkMolecule.id)).filter(
            ReactionNetworkMolecule.job_id == job.id
        ).scalar()
        
        num_reactions = db.query(func.count(ReactionNetworkReaction.id)).filter(
            ReactionNetworkReaction.job_id == job.id
        ).scalar()
        
        comparison_data.append({
            'job_id': job.id,
            'job_name': job.job_name,
            'temperature': job.temperature,
            'electrode_type': job.electrode_type,
            'voltage': job.voltage,
            'num_molecules': num_molecules,
            'num_reactions': num_reactions,
            'max_generation': job.max_generation_reached,
            'cpu_hours': job.actual_cpu_hours,
            'created_at': job.created_at.isoformat()
        })
    
    return {
        'jobs': comparison_data,
        'total_compared': len(comparison_data)
    }
