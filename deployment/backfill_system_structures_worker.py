#!/usr/bin/env python3
"""
Worker 端脚本：回填已完成 MD 任务的系统结构

这个脚本在校园网 Worker 上运行，用于处理已完成的 MD 任务，
提取它们的系统结构并上传到腾讯云后端。

使用方法:
    python backfill_system_structures_worker.py [--job-id JOB_ID] [--limit LIMIT]
    
    --job-id: 指定单个任务 ID（可选）
    --limit: 限制处理的任务数量（默认：所有）
"""

import sys
import os
from pathlib import Path
from typing import Optional, Dict, Any
import argparse
import logging
import requests
import json
import yaml

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_config(config_file: Optional[str] = None) -> Dict[str, Any]:
    """
    加载配置文件

    Args:
        config_file: 配置文件路径（可选）

    Returns:
        配置字典
    """
    if not config_file:
        # 尝试从默认位置加载
        default_paths = [
            Path(__file__).parent / 'polling_worker_config.yaml',
            Path(__file__).parent / 'polling_worker_config_tencent.yaml',
            Path.home() / '.molyte' / 'polling_worker_config.yaml',
        ]

        for path in default_paths:
            if path.exists():
                config_file = str(path)
                logger.info(f"从 {config_file} 加载配置")
                break

    if not config_file or not Path(config_file).exists():
        logger.warning("未找到配置文件，使用默认值")
        return {}

    with open(config_file, 'r') as f:
        return yaml.safe_load(f) or {}


class SystemStructureBackfiller:
    """系统结构回填器"""
    
    def __init__(self, api_base_url: str, api_token: str):
        """
        初始化回填器
        
        Args:
            api_base_url: 后端 API 基础 URL
            api_token: API 认证令牌
        """
        self.api_base_url = api_base_url.rstrip('/')
        self.api_headers = {
            'Authorization': f'Bearer {api_token}',
            'Content-Type': 'application/json'
        }
    
    def extract_system_structure(self, work_dir: Path) -> Optional[Dict[str, Any]]:
        """
        从工作目录提取系统结构

        Args:
            work_dir: 工作目录路径

        Returns:
            包含 xyz_content 和元数据的字典，或 None
        """
        try:
            # 优先查找 after_nvt 文件
            dump_files = list(work_dir.glob('**/*after_nvt*.lammpstrj'))
            if not dump_files:
                dump_files = list(work_dir.glob('**/NVT_*.lammpstrj'))
            if not dump_files:
                dump_files = list(work_dir.glob('**/*.lammpstrj'))

            if not dump_files:
                logger.warning(f"未找到轨迹文件: {work_dir}")
                return None

            dump_file = dump_files[0]
            logger.info(f"读取轨迹文件: {dump_file}")

            # 读取最后一帧
            with open(dump_file, 'r') as f:
                lines = f.readlines()

            # 找到所有帧的开始位置
            frame_starts = []
            for i, line in enumerate(lines):
                if line.startswith('ITEM: TIMESTEP'):
                    frame_starts.append(i)

            if not frame_starts:
                logger.warning(f"未找到有效的帧数据: {dump_file}")
                return None

            # 获取最后一帧
            last_frame_start = frame_starts[-1]
            last_frame_lines = lines[last_frame_start:]

            # 解析帧数据
            timestep = int(last_frame_lines[1].strip())
            num_atoms = int(last_frame_lines[3].strip())

            # 解析盒子信息 (BOX BOUNDS pp pp pp)
            box_x_line = last_frame_lines[5].strip().split()
            box_y_line = last_frame_lines[6].strip().split()
            box_z_line = last_frame_lines[7].strip().split()

            box_x = float(box_x_line[1]) - float(box_x_line[0])
            box_y = float(box_y_line[1]) - float(box_y_line[0])
            box_z = float(box_z_line[1]) - float(box_z_line[0])

            # 提取原子坐标 (跳过 ITEM: ATOMS 行)
            atom_lines = last_frame_lines[9:9+num_atoms]

            # 转换为 XYZ 格式
            xyz_lines = [str(num_atoms), f"Frame {len(frame_starts)-1} from {dump_file.name}"]

            for atom_line in atom_lines:
                parts = atom_line.strip().split()
                if len(parts) >= 6:
                    # 格式: id element mol type x y z q
                    atom_type = parts[1]  # element
                    x, y, z = float(parts[4]), float(parts[5]), float(parts[6])
                    xyz_lines.append(f"{atom_type} {x:.6f} {y:.6f} {z:.6f}")

            xyz_content = '\n'.join(xyz_lines)

            logger.info(f"成功提取 {num_atoms} 个原子，盒子尺寸: {box_x:.2f} x {box_y:.2f} x {box_z:.2f}")

            return {
                'xyz_content': xyz_content,
                'frame_index': len(frame_starts) - 1,
                'total_frames': len(frame_starts),
                'atom_count': num_atoms,
                'box': [box_x, box_y, box_z]
            }

        except Exception as e:
            logger.error(f"提取系统结构失败: {e}", exc_info=True)
            return None
    
    def upload_system_structure(self, job_id: int, structure_data: Dict[str, Any]) -> bool:
        """
        上传系统结构到后端
        
        Args:
            job_id: MD 任务 ID
            structure_data: 系统结构数据
            
        Returns:
            True 如果成功，False 如果失败
        """
        try:
            endpoint = f"{self.api_base_url}/workers/jobs/{job_id}/system_structure"
            
            response = requests.post(
                endpoint,
                json=structure_data,
                headers=self.api_headers,
                timeout=30
            )
            
            if response.status_code == 200:
                logger.info(f"✅ 任务 {job_id} 系统结构上传成功")
                return True
            else:
                logger.error(
                    f"❌ 任务 {job_id} 上传失败: "
                    f"HTTP {response.status_code} - {response.text}"
                )
                return False
                
        except Exception as e:
            logger.error(f"❌ 任务 {job_id} 上传异常: {e}")
            return False
    
    def process_job(self, job_id: int, work_dir: Path) -> bool:
        """
        处理单个任务
        
        Args:
            job_id: MD 任务 ID
            work_dir: 工作目录
            
        Returns:
            True 如果成功，False 如果失败
        """
        logger.info(f"处理任务 {job_id}...")
        
        # 检查工作目录
        if not work_dir.exists():
            logger.warning(f"工作目录不存在: {work_dir}")
            return False
        
        # 提取系统结构
        structure_data = self.extract_system_structure(work_dir)
        if not structure_data:
            logger.warning(f"无法提取任务 {job_id} 的系统结构")
            return False
        
        # 上传到后端
        return self.upload_system_structure(job_id, structure_data)


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='Worker 端回填已完成 MD 任务的系统结构'
    )
    parser.add_argument(
        '--job-id',
        type=int,
        help='指定单个任务 ID'
    )
    parser.add_argument(
        '--work-dir',
        type=str,
        help='工作目录（用于单个任务）'
    )
    parser.add_argument(
        '--config',
        type=str,
        help='配置文件路径'
    )
    parser.add_argument(
        '--api-url',
        type=str,
        help='后端 API 基础 URL（覆盖配置文件）'
    )
    parser.add_argument(
        '--api-token',
        type=str,
        help='API 认证令牌（覆盖配置文件）'
    )

    args = parser.parse_args()

    # 加载配置文件
    config = load_config(args.config)

    # 获取 API 配置
    api_url = args.api_url or config.get('api', {}).get('base_url')
    api_token = args.api_token or config.get('api', {}).get('worker_token')

    # 检查必要参数
    if not api_url or not api_token:
        logger.error("错误: 必须提供 API URL 和 Token（通过参数或配置文件）")
        sys.exit(1)

    backfiller = SystemStructureBackfiller(api_url, api_token)

    if args.job_id and args.work_dir:
        # 处理单个任务
        work_dir = Path(args.work_dir)
        success = backfiller.process_job(args.job_id, work_dir)
        sys.exit(0 if success else 1)
    else:
        logger.error("错误: 必须提供 --job-id 和 --work-dir")
        sys.exit(1)


if __name__ == '__main__':
    main()

