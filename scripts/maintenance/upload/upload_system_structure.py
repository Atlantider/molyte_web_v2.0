#!/usr/bin/env python3
"""
手动上传系统结构到数据库
从 after_nvt 轨迹文件中提取系统结构并上传到 SystemStructure 表
"""

import sys
import yaml
import requests
from pathlib import Path

# 添加backend路径
backend_path = Path(__file__).parent / "backend"
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

# 需要处理的任务
TASKS = [
    {
        'job_id': 179,
        'work_dir': '/public/home/xiaoji/molyte_web/data/md_work/EL-20251205-0002-Li-PF6-EC-DEC-PC-MD1-298K',
        'job_name': 'EL-20251205-0002-Li-PF6-EC-DEC-PC-MD1-298K'
    },
    {
        'job_id': 178,
        'work_dir': '/public/home/xiaoji/molyte_web/data/md_work/EL-20251205-0001-Li-PF6-EC-DEC-PC-MD1-298K',
        'job_name': 'EL-20251205-0001-Li-PF6-EC-DEC-PC-MD1-298K'
    }
]

def load_config():
    """加载polling worker配置"""
    try:
        config_path = Path('deployment/polling_worker_config_tencent.yaml')
        with open(config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        return config
    except Exception as e:
        print(f"❌ 加载配置失败: {e}")
        return None

def extract_system_structure(work_dir, job_name):
    """从轨迹文件提取系统结构"""
    try:
        print(f"📊 从轨迹文件提取系统结构: {work_dir}")

        work_path = Path(work_dir)

        # 查找 after_nvt.lammpstrj 文件
        after_nvt_file = work_path / f"{job_name}_after_nvt.lammpstrj"

        if not after_nvt_file.exists():
            print(f"❌ 未找到 after_nvt 文件: {after_nvt_file}")
            return None

        print(f"📁 找到 after_nvt 文件: {after_nvt_file.name}")

        # 读取最后一帧的原子坐标
        atoms = []
        in_atoms_section = False
        current_timestep = None

        with open(after_nvt_file, 'r') as f:
            lines = f.readlines()

        # 从后往前读，找到最后一帧
        i = len(lines) - 1
        while i >= 0:
            line = lines[i].strip()
            if line.startswith('ITEM: TIMESTEP'):
                current_timestep = int(lines[i+1].strip())
                break
            i -= 1

        if current_timestep is None:
            print(f"❌ 未找到时间步信息")
            return None

        # 从找到的时间步开始，读取原子信息
        i += 2  # 跳过 TIMESTEP 和时间步数值
        while i < len(lines):
            line = lines[i].strip()
            if line.startswith('ITEM: ATOMS'):
                in_atoms_section = True
                i += 1
                continue
            elif line.startswith('ITEM:'):
                in_atoms_section = False
            elif in_atoms_section and line:
                # 解析原子行: id element mol type x y z q
                parts = line.split()
                if len(parts) >= 8:
                    atom_id = int(parts[0])
                    element = parts[1]
                    x, y, z = float(parts[4]), float(parts[5]), float(parts[6])
                    atoms.append((atom_id, element, x, y, z))
            i += 1

        if not atoms:
            print(f"❌ 未找到原子坐标数据")
            return None

        # 转换为XYZ格式
        xyz_lines = [str(len(atoms)), f"Frame from {job_name} at timestep {current_timestep}"]

        for atom_id, element, x, y, z in sorted(atoms):
            xyz_lines.append(f"{element} {x:.6f} {y:.6f} {z:.6f}")

        xyz_content = '\n'.join(xyz_lines)

        result = {
            'xyz_content': xyz_content,
            'frame_index': 0,
            'total_frames': 1,
            'atom_count': len(atoms),
            'box': [0, 0, 0]  # 简化处理
        }

        print(f"✅ 系统结构提取成功:")
        print(f"   原子数: {result.get('atom_count', 0)}")
        print(f"   时间步: {current_timestep}")
        print(f"   XYZ内容长度: {len(result['xyz_content'])}")
        return result

    except Exception as e:
        print(f"❌ 提取系统结构失败: {e}")
        import traceback
        traceback.print_exc()
        return None

def upload_system_structure(config, job_id, structure_data):
    """上传系统结构到API"""
    try:
        api_base_url = config['api']['base_url']
        api_headers = {
            'Content-Type': 'application/json',
            'Authorization': f"Bearer {config['api']['worker_token']}"
        }
        
        endpoint = f"{api_base_url}/workers/jobs/{job_id}/system_structure"
        
        # 准备上传数据
        upload_data = {
            'xyz_content': structure_data['xyz_content'],
            'frame_index': structure_data.get('frame_index', 0),
            'total_frames': structure_data.get('total_frames', 1),
            'atom_count': structure_data.get('atom_count', 0),
            'box': structure_data.get('box', [0, 0, 0])
        }
        
        print(f"📡 上传系统结构到API: {endpoint}")
        print(f"   原子数: {upload_data['atom_count']}")
        print(f"   盒子: {upload_data['box']}")
        
        response = requests.post(
            endpoint,
            headers=api_headers,
            json=upload_data,
            timeout=config['api']['timeout']
        )
        
        if response.status_code == 200:
            result = response.json()
            print(f"✅ 系统结构上传成功")
            print(f"   消息: {result.get('message', 'N/A')}")
            return True
        else:
            print(f"❌ 系统结构上传失败: HTTP {response.status_code}")
            print(f"   响应: {response.text}")
            return False
            
    except Exception as e:
        print(f"❌ 上传系统结构失败: {e}")
        return False

def process_single_task(config, task):
    """处理单个任务"""
    print(f"\n{'='*60}")
    print(f"🔧 处理任务 {task['job_id']}: {task['job_name']}")
    
    # 1. 检查工作目录
    work_path = Path(task['work_dir'])
    if not work_path.exists():
        print(f"❌ 工作目录不存在: {task['work_dir']}")
        return False
    
    # 2. 提取系统结构
    structure_data = extract_system_structure(task['work_dir'], task['job_name'])
    if not structure_data:
        return False
    
    # 3. 上传系统结构
    success = upload_system_structure(config, task['job_id'], structure_data)
    
    if success:
        print(f"✅ 任务 {task['job_id']} 系统结构处理完成!")
    else:
        print(f"❌ 任务 {task['job_id']} 系统结构处理失败!")
    
    return success

def main():
    """主函数"""
    print("🔧 开始上传系统结构到数据库...")
    print(f"📋 共有 {len(TASKS)} 个任务需要处理")
    
    # 1. 加载配置
    config = load_config()
    if not config:
        return False
    
    # 2. 处理每个任务
    success_count = 0
    failed_count = 0
    
    for task in TASKS:
        try:
            if process_single_task(config, task):
                success_count += 1
            else:
                failed_count += 1
        except Exception as e:
            print(f"❌ 任务 {task['job_id']} 处理时发生异常: {e}")
            failed_count += 1
    
    # 3. 总结
    print(f"\n{'='*60}")
    print(f"🎉 系统结构上传完成!")
    print(f"   ✅ 成功上传: {success_count} 个任务")
    print(f"   ❌ 上传失败: {failed_count} 个任务")
    print(f"   📊 总计: {len(TASKS)} 个任务")
    
    if success_count > 0:
        print(f"\n🌐 前端的整体溶液结构 (System) 现在应该可以正常显示了!")
        print(f"💡 请刷新前端页面查看效果")
    
    return success_count > 0

if __name__ == "__main__":
    try:
        success = main()
        if success:
            print(f"\n🎉 脚本执行成功!")
            sys.exit(0)
        else:
            print(f"\n💥 脚本执行失败!")
            sys.exit(1)
    except KeyboardInterrupt:
        print(f"\n⚠️ 用户中断执行")
        sys.exit(1)
    except Exception as e:
        print(f"\n💥 脚本执行时发生异常: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
