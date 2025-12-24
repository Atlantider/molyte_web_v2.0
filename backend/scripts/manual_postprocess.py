#!/usr/bin/env python3
"""
手动后处理脚本

用于手动触发MD任务的后处理，不依赖Celery
"""

import sys
import os
from pathlib import Path
import json
import logging

# 添加项目根目录到Python路径
sys.path.insert(0, str(Path(__file__).parent.parent))

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def manual_postprocess(job_id: int):
    """手动执行后处理"""
    
    # 避免导入Celery相关模块
    from app.database import SessionLocal
    from app.models.job import MDJob
    from app.models.result import ResultSummary, RDFResult, MSDResult
    
    db = SessionLocal()
    try:
        # 获取任务
        job = db.query(MDJob).filter(MDJob.id == job_id).first()
        if not job:
            logger.error(f'任务{job_id}不存在')
            return False
            
        logger.info(f'任务{job_id}状态: {job.status}')
        logger.info(f'工作目录: {job.work_dir}')
        
        work_dir = Path(job.work_dir)
        if not work_dir.exists():
            logger.error(f'工作目录不存在: {work_dir}')
            return False
            
        logger.info(f'工作目录存在，开始后处理...')
        
        # 检查是否已有结果摘要
        existing_summary = db.query(ResultSummary).filter(ResultSummary.md_job_id == job_id).first()
        if existing_summary:
            logger.info('已存在结果摘要，删除重新创建')
            db.delete(existing_summary)
            db.commit()
        
        # 删除已有的RDF和MSD结果
        existing_rdf = db.query(RDFResult).filter(RDFResult.md_job_id == job_id).all()
        for rdf in existing_rdf:
            db.delete(rdf)
        
        existing_msd = db.query(MSDResult).filter(MSDResult.md_job_id == job_id).all()
        for msd in existing_msd:
            db.delete(msd)
        
        db.commit()
        
        # 创建结果摘要
        result_summary = ResultSummary(
            md_job_id=job_id,
            final_structure_path=str(work_dir / 'final_structure.data'),
            trajectory_path=str(work_dir / 'trajectory.dcd'),
            log_file_path=str(work_dir / 'log.lammps'),
            atom_mapping_path=str(work_dir / 'atom_mapping.json')
        )
        db.add(result_summary)
        db.commit()
        logger.info('✅ 结果摘要已创建')
        
        # 计算RDF
        try:
            rdf_results = calculate_rdf_from_lammps(work_dir, job_id)
            logger.info(f'✅ RDF计算完成，共{len(rdf_results)}个结果')
            
            # 保存RDF结果到数据库
            for rdf_data in rdf_results:
                rdf_result = RDFResult(
                    md_job_id=job_id,
                    center_species=rdf_data['center_label'],
                    shell_species=rdf_data['target_label'],
                    r_values=rdf_data['r_values'],
                    gr_values=rdf_data['gr_values'],
                    coordination_number=rdf_data.get('coordination_number'),
                    first_peak_position=rdf_data.get('first_peak_position'),
                    first_minimum_position=rdf_data.get('first_minimum_position')
                )
                db.add(rdf_result)
            
            db.commit()
            logger.info('✅ RDF结果已保存到数据库')
            
        except Exception as e:
            logger.error(f'❌ RDF计算失败: {e}')
            import traceback
            traceback.print_exc()
        
        # 计算MSD（如果有数据）
        try:
            msd_results = calculate_msd_from_lammps(work_dir, job_id)
            if msd_results:
                logger.info(f'✅ MSD计算完成，共{len(msd_results)}个结果')
                
                for msd_data in msd_results:
                    msd_result = MSDResult(
                        md_job_id=job_id,
                        species=msd_data['species'],
                        time_values=msd_data['time_values'],
                        msd_values=msd_data['msd_values'],
                        diffusion_coefficient=msd_data.get('diffusion_coefficient')
                    )
                    db.add(msd_result)
                
                db.commit()
                logger.info('✅ MSD结果已保存到数据库')
            else:
                logger.info('未找到MSD数据')
                
        except Exception as e:
            logger.error(f'❌ MSD计算失败: {e}')
            import traceback
            traceback.print_exc()
        
        logger.info('✅ 后处理完成')
        return True
        
    except Exception as e:
        logger.error(f'后处理失败: {e}')
        import traceback
        traceback.print_exc()
        return False
        
    finally:
        db.close()


def calculate_rdf_from_lammps(work_dir: Path, job_id: int):
    """从LAMMPS输出文件读取RDF数据"""
    
    rdf_file = work_dir / "out_rdf.dat"
    if not rdf_file.exists():
        logger.warning(f"RDF文件不存在: {rdf_file}")
        return []
    
    logger.info(f"读取RDF文件: {rdf_file}")
    
    results = []
    
    try:
        with open(rdf_file, 'r') as f:
            lines = f.readlines()
        
        # 解析RDF数据
        # LAMMPS RDF输出格式: # Bin r g(r) coordination_number
        current_pair = None
        r_values = []
        gr_values = []
        
        for line in lines:
            line = line.strip()
            
            # 跳过注释行和空行
            if not line or line.startswith('#'):
                # 检查是否是新的RDF对
                if 'RDF' in line and '-' in line:
                    # 保存之前的数据
                    if current_pair and r_values and gr_values:
                        center, shell = current_pair.split('-')
                        results.append({
                            'center_label': center.strip(),
                            'target_label': shell.strip(),
                            'r_values': r_values.copy(),
                            'gr_values': gr_values.copy(),
                            'coordination_number': None,
                            'first_peak_position': None,
                            'first_minimum_position': None
                        })
                    
                    # 开始新的RDF对
                    if '-' in line:
                        parts = line.split()
                        for part in parts:
                            if '-' in part and not part.startswith('#'):
                                current_pair = part
                                r_values = []
                                gr_values = []
                                break
                continue
            
            # 解析数据行
            try:
                parts = line.split()
                if len(parts) >= 3:
                    r = float(parts[1])
                    gr = float(parts[2])
                    r_values.append(r)
                    gr_values.append(gr)
            except (ValueError, IndexError):
                continue
        
        # 保存最后一个RDF对
        if current_pair and r_values and gr_values:
            center, shell = current_pair.split('-')
            results.append({
                'center_label': center.strip(),
                'target_label': shell.strip(),
                'r_values': r_values,
                'gr_values': gr_values,
                'coordination_number': None,
                'first_peak_position': None,
                'first_minimum_position': None
            })
        
        logger.info(f"成功解析{len(results)}个RDF对")
        
    except Exception as e:
        logger.error(f"解析RDF文件失败: {e}")
        return []
    
    return results


def calculate_msd_from_lammps(work_dir: Path, job_id: int):
    """从LAMMPS输出文件读取MSD数据"""
    
    msd_file = work_dir / "out_msd.dat"
    if not msd_file.exists():
        logger.warning(f"MSD文件不存在: {msd_file}")
        return []
    
    logger.info(f"读取MSD文件: {msd_file}")
    
    # 这里可以实现MSD数据解析
    # 暂时返回空列表
    return []


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("用法: python manual_postprocess.py <job_id>")
        sys.exit(1)
    
    try:
        job_id = int(sys.argv[1])
        success = manual_postprocess(job_id)
        sys.exit(0 if success else 1)
    except ValueError:
        print("错误: job_id必须是整数")
        sys.exit(1)
