#!/usr/bin/env python3
"""
从COS下载文件并执行后处理

用于腾讯云服务器从COS下载MD结果文件，然后执行后处理分析
"""

import sys
import os
import tempfile
import shutil
from pathlib import Path
import json
import logging

# 添加项目根目录到Python路径
sys.path.insert(0, str(Path(__file__).parent.parent))

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def postprocess_from_cos(job_id: int):
    """从COS下载文件并执行后处理"""

    from app.database import SessionLocal
    from app.models.job import MDJob
    from app.models.result import ResultSummary, RDFResult, MSDResult
    import yaml

    # 初始化COS客户端
    try:
        from qcloud_cos import CosConfig, CosS3Client

        # 从polling_worker_config.yaml读取COS配置
        config_path = Path(__file__).parent.parent.parent / "deployment" / "polling_worker_config_tencent.yaml"
        if not config_path.exists():
            config_path = Path(__file__).parent.parent.parent / "deployment" / "polling_worker_config.yaml"

        if not config_path.exists():
            logger.error(f"配置文件不存在: {config_path}")
            return False

        with open(config_path, 'r', encoding='utf-8') as f:
            worker_config = yaml.safe_load(f)

        if 'cos' not in worker_config:
            logger.error("配置文件中未找到COS配置")
            return False

        cos_config = worker_config['cos']
        config = CosConfig(
            Region=cos_config['region'],
            SecretId=cos_config['secret_id'],
            SecretKey=cos_config['secret_key'],
            Scheme='https'
        )
        cos_client = CosS3Client(config)
        bucket_name = cos_config['bucket']

        logger.info(f"COS客户端已初始化 (Bucket: {bucket_name})")

    except Exception as e:
        logger.error(f"初始化COS客户端失败: {e}")
        return False
    
    db = SessionLocal()
    temp_dir = None
    
    try:
        # 获取任务
        job = db.query(MDJob).filter(MDJob.id == job_id).first()
        if not job:
            logger.error(f'任务{job_id}不存在')
            return False
            
        logger.info(f'任务{job_id}状态: {job.status}')
        
        # 创建临时目录
        temp_dir = Path(tempfile.mkdtemp(prefix=f'md_postprocess_{job_id}_'))
        logger.info(f'创建临时目录: {temp_dir}')
        
        # 从COS下载文件
        cos_prefix = f"MD_results/{job_id}/"
        
        # 列出COS中的文件
        try:
            response = cos_client.list_objects(
                Bucket=bucket_name,
                Prefix=cos_prefix
            )
            
            if 'Contents' not in response:
                logger.error(f'COS中未找到任务{job_id}的文件')
                return False
            
            files = response['Contents']
            logger.info(f'找到{len(files)}个文件')
            
            # 下载所有文件
            for file_info in files:
                cos_key = file_info['Key']
                file_name = cos_key.replace(cos_prefix, '')
                
                if not file_name:  # 跳过目录
                    continue
                
                local_path = temp_dir / file_name
                
                logger.info(f'下载文件: {cos_key} -> {local_path}')
                
                try:
                    cos_client.download_file(
                        Bucket=bucket_name,
                        Key=cos_key,
                        DestFilePath=str(local_path)
                    )
                except Exception as e:
                    logger.warning(f'下载文件失败 {cos_key}: {e}')
                    continue
            
        except Exception as e:
            logger.error(f'从COS下载文件失败: {e}')
            return False
        
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
            final_structure_path=f"{cos_prefix}final_structure.data",
            trajectory_path=f"{cos_prefix}trajectory.dcd", 
            log_file_path=f"{cos_prefix}log.lammps",
            atom_mapping_path=f"{cos_prefix}atom_mapping.json"
        )
        db.add(result_summary)
        db.commit()
        logger.info('✅ 结果摘要已创建')
        
        # 计算RDF
        try:
            rdf_results = calculate_rdf_from_files(temp_dir, job_id)
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
            msd_results = calculate_msd_from_files(temp_dir, job_id)
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
        # 清理临时目录
        if temp_dir and temp_dir.exists():
            shutil.rmtree(temp_dir)
            logger.info(f'清理临时目录: {temp_dir}')


def calculate_rdf_from_files(work_dir: Path, job_id: int):
    """从下载的文件中读取RDF数据"""
    
    # 查找RDF文件
    rdf_files = list(work_dir.glob("*rdf*.dat")) + list(work_dir.glob("*RDF*.dat"))
    
    if not rdf_files:
        logger.warning(f"未找到RDF文件")
        return []
    
    rdf_file = rdf_files[0]
    logger.info(f"读取RDF文件: {rdf_file}")
    
    results = []
    
    try:
        with open(rdf_file, 'r') as f:
            lines = f.readlines()
        
        # 解析RDF数据
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


def calculate_msd_from_files(work_dir: Path, job_id: int):
    """从下载的文件中读取MSD数据"""
    
    # 查找MSD文件
    msd_files = list(work_dir.glob("*msd*.dat")) + list(work_dir.glob("*MSD*.dat"))
    
    if not msd_files:
        logger.warning(f"未找到MSD文件")
        return []
    
    # 这里可以实现MSD数据解析
    # 暂时返回空列表
    return []


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("用法: python postprocess_from_cos.py <job_id>")
        sys.exit(1)
    
    try:
        job_id = int(sys.argv[1])
        success = postprocess_from_cos(job_id)
        sys.exit(0 if success else 1)
    except ValueError:
        print("错误: job_id必须是整数")
        sys.exit(1)
