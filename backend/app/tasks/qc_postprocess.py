"""
QC 后处理 Celery Worker

负责处理 QC 计算完成后的结果提取和ESP可视化
"""

import logging
import os
import re
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from sqlalchemy.orm import Session
    from app.models.qc import QCJob, QCResult

logger = logging.getLogger(__name__)

# 延迟导入 Celery 和数据库相关模块，避免在校园网 polling_worker 中导入时出错
# 这些模块只在 Celery 任务中使用，不在可视化函数中使用
_celery_app = None
_SessionLocal = None
_QCJob = None
_QCResult = None
_MoleculeQCCache = None
_QCJobStatus = None


def get_slurm_job_cpu_hours(slurm_job_id: str) -> float:
    """
    从 Slurm 获取任务的 CPU 核时数

    使用 sacct 获取已完成任务的 CPUTimeRAW

    Args:
        slurm_job_id: Slurm 任务 ID

    Returns:
        CPU 核时数（小时）
    """
    try:
        result = subprocess.run(
            ['sacct', '-j', slurm_job_id, '-o', 'CPUTimeRAW', '-n', '-X'],
            capture_output=True,
            text=True,
            timeout=10
        )

        if result.returncode == 0 and result.stdout.strip():
            try:
                cpu_time_seconds = int(result.stdout.strip().split()[0])
                if cpu_time_seconds > 0:
                    cpu_hours = cpu_time_seconds / 3600.0
                    logger.debug(f"Job {slurm_job_id}: CPUTimeRAW={cpu_time_seconds}s, CPU hours={cpu_hours:.2f}h")
                    return cpu_hours
            except (ValueError, IndexError) as e:
                logger.warning(f"Failed to parse CPUTimeRAW for job {slurm_job_id}: {e}")

        return 0.0

    except subprocess.TimeoutExpired:
        logger.warning(f"sacct timeout for job {slurm_job_id}")
        return 0.0
    except Exception as e:
        logger.error(f"Failed to get CPU hours for job {slurm_job_id}: {e}")
        return 0.0

def _get_celery_app():
    global _celery_app
    if _celery_app is None:
        from app.celery_app import celery_app
        _celery_app = celery_app
    return _celery_app

def _get_db_models():
    global _SessionLocal, _QCJob, _QCResult, _MoleculeQCCache, _QCJobStatus
    if _SessionLocal is None:
        from app.database import SessionLocal
        from app.models.qc import QCJob, QCResult, MoleculeQCCache, QCJobStatus
        _SessionLocal = SessionLocal
        _QCJob = QCJob
        _QCResult = QCResult
        _MoleculeQCCache = MoleculeQCCache
        _QCJobStatus = QCJobStatus
    return _SessionLocal, _QCJob, _QCResult, _MoleculeQCCache, _QCJobStatus

# Hartree to eV conversion factor
HARTREE_TO_EV = 27.2114


class DatabaseTask:
    """自动管理数据库会话的任务基类"""
    _db = None

    def after_return(self, *args, **kwargs):
        if self._db is not None:
            self._db.close()
            self._db = None

    @property
    def db(self):
        if self._db is None:
            SessionLocal, _, _, _, _ = _get_db_models()
            self._db = SessionLocal()
        return self._db


def extract_gaussian_results(log_file: str) -> Dict[str, Any]:
    """
    从Gaussian输出文件提取计算结果
    
    Args:
        log_file: Gaussian输出文件路径
        
    Returns:
        Dict containing energy, HOMO, LUMO values
    """
    results = {
        "energy_au": None,
        "homo": None,
        "lumo": None,
    }
    
    if not os.path.exists(log_file):
        logger.error(f"Log file not found: {log_file}")
        return results
    
    # 正则表达式模式
    energy_pattern = re.compile(r"SCF Done:\s+E\([A-Z0-9]+\)\s+=\s+([-\d.]+)")
    alpha_occ_pattern = re.compile(r"Alpha\s+occ\.\seigenvalues\s+--\s+([-\d.]+(?:\s+[-\d.]+)*)")
    alpha_virt_pattern = re.compile(r"Alpha\s+virt\.\seigenvalues\s+--\s+([-\d.]+(?:\s+[-\d.]+)*)")
    
    last_energy = None
    last_homo = None
    last_lumo = None
    
    # 使用errors='replace'来处理可能的编码问题（Gaussian输出可能包含非UTF-8字符）
    with open(log_file, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
        for i in range(len(lines) - 1):
            # 匹配SCF能量
            match_energy = energy_pattern.search(lines[i])
            if match_energy:
                last_energy = float(match_energy.group(1))
            
            # 匹配HOMO和LUMO
            match_occ = alpha_occ_pattern.search(lines[i])
            match_virt = alpha_virt_pattern.search(lines[i + 1]) if i + 1 < len(lines) else None
            
            if match_occ and match_virt:
                occ_values = match_occ.group(1).split()
                virt_values = match_virt.group(1).split()
                
                if occ_values:
                    last_homo = float(occ_values[-1])
                if virt_values:
                    last_lumo = float(virt_values[0])
    
    results["energy_au"] = last_energy
    results["homo"] = last_homo
    results["lumo"] = last_lumo
    
    return results


def extract_esp_values(surfanalysis_file: str) -> Tuple[Optional[float], Optional[float]]:
    """
    从Multiwfn的surfanalysis.txt提取ESP极值
    
    Returns:
        Tuple of (ESP_min, ESP_max) in kcal/mol
    """
    if not os.path.exists(surfanalysis_file):
        logger.warning(f"Surface analysis file not found: {surfanalysis_file}")
        return None, None
    
    min_pattern = r"\*\s*(\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)"
    max_pattern = r"\s*(\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)"
    
    with open(surfanalysis_file, 'r') as f:
        data = f.read()
    
    min_matches = re.findall(min_pattern, data)
    max_matches = re.findall(max_pattern, data)
    
    esp_min = None
    esp_max = None
    
    if min_matches:
        min_values = [float(match[3]) for match in min_matches]
        esp_min = min(min_values)
    
    if max_matches:
        max_values = [float(match[3]) for match in max_matches]
        esp_max = max(max_values)
    
    return esp_min, esp_max


def _generate_esp_visualization_vmd(work_dir: Path, molecule_name: str, fchk_file: str) -> Optional[str]:
    """
    使用VMD生成ESP可视化图像

    应用用户提供的VMD参数配置：
    - material change mirror Opaque 0.15
    - material change outline Opaque 4.000000
    - material change outlinewidth Opaque 0.5
    - material change ambient Glossy 0.1
    - material change diffuse Glossy 0.600000
    - material change opacity Glossy 0.75
    - material change ambient Opaque 0.08
    - material change mirror Opaque 0.0
    - material change shininess Glossy 1.0
    - mol modcolor 1 top ColorID 12 (蓝色: 0.129 0.580 0.812)
    - mol modcolor 2 top ColorID 22 (红色: 0.925 0.627 0.635)
    - display distance -7.0
    - display height 10
    - light 3 on

    Args:
        work_dir: 工作目录
        molecule_name: 分子名称
        fchk_file: fchk文件路径

    Returns:
        ESP图像路径，失败返回None
    """
    try:
        # 创建简化的 Multiwfn 输入序列生成 ESP cube 文件
        # 12 = 计算分子表面上的性质
        # 1 = ESP (静电势)
        # 创建Multiwfn输入文件（使用成功例子中的完整输入序列）
        txt_content = '''12
0
1
-1
-1
5
1
3
2
0
5
12
1
2'''
        txt_path = work_dir / "ESPiso.txt"
        with open(txt_path, 'w') as f:
            f.write(txt_content)

        # 创建VMD脚本 - 使用Tachyon渲染器用于无头模式
        # 应用用户提供的参数，支持自适应缩放
        vmd_content = '''# ESP可视化脚本 - 使用Tachyon渲染器（自适应缩放）
color scale method BWR
color Display Background white
axes location Off
display depthcue off
display projection Orthographic
display nearclip set 0.01

# 应用用户提供的材质参数
material change mirror Opaque 0.15
material change outline Opaque 4.000000
material change outlinewidth Opaque 0.5
material change ambient Glossy 0.1
material change diffuse Glossy 0.600000
material change opacity Glossy 0.75
material change ambient Opaque 0.08
material change mirror Opaque 0.0
material change shininess Glossy 1.0

# 应用用户提供的颜色配置
mol modcolor 1 top ColorID 12
color change rgb "12" 0.129 0.580 0.812
mol modcolor 2 top ColorID 22
color change rgb "22" 0.925 0.627 0.635

set nsystem 1
set colorlow -0.03
set colorhigh 0.03

for {set i 1} {$i<=$nsystem} {incr i} {
set id [expr $i-1]
mol new density$i.cub
mol addfile ESP$i.cub
mol modstyle 0 $id CPK 1.000000 0.300000 22.000000 22.000000
mol addrep $id
mol modstyle 1 $id Isosurface 0.001000 0 0 0 1 1
mol modmaterial 1 $id Transparent
mol modcolor 1 $id Volume 1
mol scaleminmax $id 1 $colorlow $colorhigh
}

# 自适应视角调整 - 根据分子大小自动缩放
display resetview
mol top 0

# 获取分子的边界框
set sel [atomselect top all]
set minmax [measure minmax $sel]
set min_coords [lindex $minmax 0]
set max_coords [lindex $minmax 1]
$sel delete

# 计算分子在三个方向的尺寸
set dx [expr {[lindex $max_coords 0] - [lindex $min_coords 0]}]
set dy [expr {[lindex $max_coords 1] - [lindex $min_coords 1]}]
set dz [expr {[lindex $max_coords 2] - [lindex $min_coords 2]}]

# 获取最大维度（考虑等密度面会比原子位置大约3-4埃）
set iso_padding 6.0
set max_dim [expr {max($dx, max($dy, $dz)) + $iso_padding}]

# 计算自适应缩放比例
# 对于小分子（<20埃）：目标占画面60%
# 对于中等分子（20-30埃）：目标占画面50%
# 对于大分子（>30埃）：目标占画面40%
set target_view 15.0
if {$max_dim < 20.0} {
    set fill_ratio 0.60
} elseif {$max_dim < 30.0} {
    set fill_ratio 0.50
} else {
    set fill_ratio 0.40
}

# 检查分子是否为细长型（长宽比 > 2.5）
set aspect_ratio_xy [expr {$dx / $dy}]
set aspect_ratio_xz [expr {$dx / $dz}]
set aspect_ratio_yz [expr {$dy / $dz}]
set max_aspect_ratio [expr {max($aspect_ratio_xy, max($aspect_ratio_xz, $aspect_ratio_yz))}]

# 对于细长分子，增加填充比例以确保可见性
if {$max_aspect_ratio > 2.5} {
    set fill_ratio [expr {$fill_ratio + 0.15}]
    if {$fill_ratio > 0.85} {
        set fill_ratio 0.85
    }
    puts "ESP - Detected elongated molecule (aspect ratio: $max_aspect_ratio), increased fill_ratio to $fill_ratio"
}

set auto_scale [expr {($target_view * $fill_ratio) / $max_dim}]

# 限制缩放范围，避免极端情况
# 对于大分子，允许更小的缩放值
if {$auto_scale > 1.2} {
    set auto_scale 1.2
}
if {$auto_scale < 0.15} {
    set auto_scale 0.15
}

puts "ESP - Molecule dimensions: $dx x $dy x $dz Angstrom"
puts "ESP - Max dimension with padding: $max_dim Angstrom"
puts "ESP - Fill ratio: $fill_ratio"
puts "ESP - Auto scale factor: $auto_scale"

# 应用用户提供的显示参数
display distance -7.0
display height 10

# 旋转到较好的观察角度（来自成功例子）
rotate x by 15
rotate y by 25
rotate z by 5

# 应用自适应缩放
scale by $auto_scale

# 应用用户提供的光照设置
light 3 on

# 使用Tachyon渲染器生成图片（支持无头模式）
render Tachyon ESP.dat
exit'''
        vmd_path = work_dir / "ESPiso_vmd.vmd"
        with open(vmd_path, 'w') as f:
            f.write(vmd_content)

        # 运行Multiwfn生成cube文件
        fchk_basename = os.path.basename(fchk_file)
        esp_script = f"""#!/bin/bash
cd {work_dir}
Multiwfn {fchk_basename} -ESPrhoiso 0.001 < ESPiso.txt
mv -f density.cub density1.cub
mv -f totesp.cub ESP1.cub
"""
        script_path = work_dir / "run_esp.sh"
        with open(script_path, 'w') as f:
            f.write(esp_script)
        os.chmod(script_path, 0o755)

        # 执行Multiwfn
        result = subprocess.run(
            ["bash", str(script_path)],
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            timeout=600
        )

        if result.returncode != 0:
            logger.warning(f"Multiwfn cube generation failed: {result.stderr}")
            return None

        # 运行VMD生成Tachyon数据文件
        result = subprocess.run(
            ["vmd", "-dispdev", "text", "-e", "ESPiso_vmd.vmd"],
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            timeout=120
        )

        # 使用Tachyon渲染器将.dat转换为.tga
        tachyon_dat = work_dir / "ESP.dat"
        tachyon_tga = work_dir / "ESP.tga"
        esp_image = work_dir / "ESP.png"

        if tachyon_dat.exists():
            # 查找tachyon可执行文件路径
            tachyon_cmd = "tachyon_LINUXAMD64"  # VMD自带的tachyon
            tachyon_path = "/usr/local/lib/vmd/tachyon_LINUXAMD64"
            if os.path.exists(tachyon_path):
                tachyon_cmd = tachyon_path

            # 运行Tachyon渲染（添加完整参数）
            result = subprocess.run(
                [tachyon_cmd, "-aasamples", "12", str(tachyon_dat),
                 "-format", "TGA", "-res", "1024", "768", "-o", str(tachyon_tga)],
                cwd=str(work_dir),
                capture_output=True,
                text=True,
                timeout=120
            )

            if tachyon_tga.exists():
                try:
                    result = subprocess.run(
                        ["convert", str(tachyon_tga), str(esp_image)],
                        cwd=str(work_dir),
                        capture_output=True,
                        text=True,
                        timeout=60
                    )
                except Exception as e:
                    logger.warning(f"ImageMagick convert failed: {e}, trying Python PIL")
                    # 尝试用PIL转换
                    try:
                        from PIL import Image as PILImage
                        img = PILImage.open(str(tachyon_tga))
                        img.save(str(esp_image))
                    except Exception as e2:
                        logger.warning(f"PIL convert failed: {e2}")

        if esp_image.exists() and esp_image.stat().st_size > 0:
            return str(esp_image)

        return None

    except Exception as e:
        logger.error(f"VMD ESP visualization failed: {e}")
        return None


def generate_esp_visualization(work_dir: Path, molecule_name: str, fchk_file: str) -> Optional[str]:
    """
    生成ESP可视化图像，使用Multiwfn生成cube文件，然后用VMD渲染
    支持自适应画幅和图像大小

    Args:
        work_dir: 工作目录
        molecule_name: 分子名称
        fchk_file: fchk文件路径

    Returns:
        ESP图像路径，失败返回None
    """
    try:
        image_file = work_dir / f"ESP_{molecule_name}.png"

        # 使用Multiwfn生成ESP cube文件
        # 12: 绘制分子结构和性质
        # 0: 返回主菜单
        # 1: 绘制等值面
        # -1: 返回
        # 5: 指定等值面类型（ESP）
        # 1: 等值面
        # 3: 等值面值
        # 2: 中等质量网格
        # 0: 返回
        # 使用验证过的完整 Multiwfn 输入序列（来自成功的 QC-1485 例子）
        txt_content = '''12
0
1
-1
-1
5
1
3
2
0
5
12
1
2'''
        txt_path = work_dir / "ESPiso.txt"
        with open(txt_path, 'w') as f:
            f.write(txt_content)

        # 创建 run_esp.sh 脚本，完全照抄成功例子的方式
        # 只使用文件名，因为脚本会 cd 到工作目录
        fchk_filename = os.path.basename(fchk_file)
        script_content = f'''#!/bin/bash
cd {work_dir}
Multiwfn {fchk_filename} -ESPrhoiso 0.001 < ESPiso.txt
mv -f density.cub density1.cub
mv -f totesp.cub ESP1.cub
'''
        script_path = work_dir / "run_esp.sh"
        with open(script_path, 'w') as f:
            f.write(script_content)

        # 使脚本可执行
        os.chmod(script_path, 0o755)

        # 运行脚本生成cube文件
        result = subprocess.run(
            ["bash", str(script_path)],
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            timeout=300
        )

        if result.returncode != 0:
            logger.warning(f"Multiwfn ESP cube generation failed: {result.stderr}")
            return None

        # 检查脚本是否成功生成并重命名了cube文件
        density1_cube = work_dir / "density1.cub"
        esp1_cube = work_dir / "ESP1.cub"

        if not density1_cube.exists() or not esp1_cube.exists():
            logger.warning(f"Multiwfn script did not generate ESP cube files")
            return None

        # 使用VMD渲染ESP图像
        return _generate_esp_visualization_vmd(work_dir, molecule_name, fchk_file)

    except Exception as e:
        logger.error(f"ESP visualization failed: {e}")
        return None


def generate_orbital_visualization(work_dir: Path, fchk_file: str, orbital_type: str) -> Optional[str]:
    """
    使用VMD生成HOMO或LUMO轨道可视化图像

    应用用户提供的VMD参数配置：
    - material change mirror Opaque 0.15
    - material change outline Opaque 4.000000
    - material change outlinewidth Opaque 0.5
    - material change ambient Glossy 0.1
    - material change diffuse Glossy 0.600000
    - material change opacity Glossy 0.75
    - material change ambient Opaque 0.08
    - material change mirror Opaque 0.0
    - material change shininess Glossy 1.0
    - mol modcolor 1 top ColorID 12 (蓝色: 0.129 0.580 0.812)
    - mol modcolor 2 top ColorID 22 (红色: 0.925 0.627 0.635)
    - display distance -7.0
    - display height 10
    - light 3 on

    Args:
        work_dir: 工作目录
        fchk_file: fchk文件路径
        orbital_type: "HOMO" 或 "LUMO"

    Returns:
        图像路径，失败返回None
    """
    return _generate_orbital_visualization_vmd(work_dir, fchk_file, orbital_type)


def _generate_orbital_visualization_vmd(work_dir: Path, fchk_file: str, orbital_type: str) -> Optional[str]:
    """
    使用VMD渲染轨道图像（备用方案）

    应用用户提供的VMD参数配置：
    - material change mirror Opaque 0.15
    - material change outline Opaque 4.000000
    - material change outlinewidth Opaque 0.5
    - material change ambient Glossy 0.1
    - material change diffuse Glossy 0.600000
    - material change opacity Glossy 0.75
    - material change ambient Opaque 0.08
    - material change mirror Opaque 0.0
    - material change shininess Glossy 1.0
    - mol modcolor 1 top ColorID 12 (蓝色: 0.129 0.580 0.812)
    - mol modcolor 2 top ColorID 22 (红色: 0.925 0.627 0.635)
    - display distance -7.0
    - display height 10
    - light 3 on
    """
    try:
        orbital_lower = orbital_type.lower()
        cube_file = work_dir / f"{orbital_lower}.cub"
        image_file = work_dir / f"{orbital_type}_dimer.png"

        # 先检查是否已有cube文件
        if not cube_file.exists():
            # 使用Multiwfn生成轨道cube文件
            orbital_key = "h" if orbital_type == "HOMO" else "l"

            txt_content = f'''5
4
{orbital_key}
2
2
'''
            txt_path = work_dir / f"{orbital_lower}_gen.txt"
            with open(txt_path, 'w') as f:
                f.write(txt_content)

            # 运行Multiwfn生成cube
            multiwfn_cmd = "/home/iei/share/software/Multiwfn_3.8_dev_bin_Linux/Multiwfn"
            if os.path.exists(multiwfn_cmd):
                result = subprocess.run(
                    [multiwfn_cmd, fchk_file],
                    stdin=open(txt_path, 'r'),
                    cwd=str(work_dir),
                    capture_output=True,
                    text=True,
                    timeout=300
                )

                # Multiwfn输出为 MOvalue.cub，需要重命名
                mo_cube = work_dir / "MOvalue.cub"
                if mo_cube.exists():
                    mo_cube.rename(cube_file)
                    logger.info(f"Generated {orbital_type} cube: {cube_file}")

        if not cube_file.exists():
            logger.warning(f"{orbital_type} cube file not found")
            return None

        # 使用VMD渲染轨道图像，应用用户提供的参数
        vmd_script = f'''# {orbital_type}轨道可视化脚本 - 使用Tachyon渲染器（自适应缩放）
color Display Background white
axes location Off
display depthcue off
display projection Orthographic
display nearclip set 0.01
light 2 on
light 3 on

# 加载分子和轨道
mol new {cube_file.name}
mol modstyle 0 0 CPK 0.800000 0.300000 22.000000 22.000000

# 添加正相位等值面（蓝色 - ColorID 12: 0.129 0.580 0.812）
mol addrep 0
mol modstyle 1 0 Isosurface 0.02 0 0 0 1 1
mol modcolor 1 0 ColorID 12
color change rgb "12" 0.129 0.580 0.812
mol modmaterial 1 0 Glossy
material change ambient Glossy 0.1
material change diffuse Glossy 0.600000
material change opacity Glossy 0.75
material change shininess Glossy 1.0

# 添加负相位等值面（红色 - ColorID 22: 0.925 0.627 0.635）
mol addrep 0
mol modstyle 2 0 Isosurface -0.02 0 0 0 1 1
mol modcolor 2 0 ColorID 22
color change rgb "22" 0.925 0.627 0.635
mol modmaterial 2 0 Glossy

# 应用Opaque材质参数
material change mirror Opaque 0.15
material change outline Opaque 4.000000
material change outlinewidth Opaque 0.5
material change ambient Opaque 0.08
material change mirror Opaque 0.0

# 自适应视角调整
display resetview
mol top 0

# 获取分子的边界框
set sel [atomselect top all]
set minmax [measure minmax $sel]
set min_coords [lindex $minmax 0]
set max_coords [lindex $minmax 1]
$sel delete

# 计算分子在三个方向的尺寸
set dx [expr {{[lindex $max_coords 0] - [lindex $min_coords 0]}}]
set dy [expr {{[lindex $max_coords 1] - [lindex $min_coords 1]}}]
set dz [expr {{[lindex $max_coords 2] - [lindex $min_coords 2]}}]

# 获取最大维度（考虑轨道云会比原子位置大约5-6埃）
set orb_padding 8.0
set max_dim [expr {{max($dx, max($dy, $dz)) + $orb_padding}}]

# 计算自适应缩放比例 - HOMO/LUMO使用较激进的缩放因子
# 小分子（<20埃）：目标占画面65%
# 中等分子（20-30埃）：目标占画面60%
# 大分子（>30埃）：目标占画面50%
set target_view 15.0
if {{$max_dim < 20.0}} {{
    set fill_ratio 0.65
}} elseif {{$max_dim < 30.0}} {{
    set fill_ratio 0.60
}} else {{
    set fill_ratio 0.50
}}

# 检查分子是否为细长型（长宽比 > 2.5）
set aspect_ratio_xy [expr {{$dx / $dy}}]
set aspect_ratio_xz [expr {{$dx / $dz}}]
set aspect_ratio_yz [expr {{$dy / $dz}}]
set max_aspect_ratio [expr {{max($aspect_ratio_xy, max($aspect_ratio_xz, $aspect_ratio_yz))}}]

# 对于细长分子，增加填充比例以确保可见性
if {{$max_aspect_ratio > 2.5}} {{
    set fill_ratio [expr {{$fill_ratio + 0.15}}]
    if {{$fill_ratio > 0.90}} {{
        set fill_ratio 0.90
    }}
    puts "{orbital_type} - Detected elongated molecule (aspect ratio: $max_aspect_ratio), increased fill_ratio to $fill_ratio"
}}

set auto_scale [expr {{($target_view * $fill_ratio) / $max_dim}}]

# HOMO/LUMO 缩放范围调整（基于成功例子的改进）
# 对于大分子，允许更小的缩放值
if {{$auto_scale > 1.2}} {{
    set auto_scale 1.2
}}
if {{$auto_scale < 0.15}} {{
    set auto_scale 0.15
}}

puts "{orbital_type} - Molecule dimensions: $dx x $dy x $dz Angstrom"
puts "{orbital_type} - Max dimension with padding: $max_dim Angstrom"
puts "{orbital_type} - Fill ratio: $fill_ratio"
puts "{orbital_type} - Auto scale factor: $auto_scale"

# 应用用户提供的显示参数
display distance -7.0
display height 10

# 旋转到较好的观察角度（轨道图像使用稍微不同的角度）
rotate x by 20
rotate y by 30
rotate z by 10

# 应用自适应缩放
scale by $auto_scale

# 使用Tachyon渲染器生成图片
render Tachyon {orbital_type}.dat
exit
'''

        vmd_path = work_dir / f"{orbital_lower}_render.vmd"
        with open(vmd_path, 'w') as f:
            f.write(vmd_script)

        # 运行VMD
        result = subprocess.run(
            ["vmd", "-dispdev", "text", "-e", str(vmd_path)],
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            timeout=120
        )

        # 使用Tachyon渲染
        tachyon_dat = work_dir / f"{orbital_type}.dat"
        tachyon_tga = work_dir / f"{orbital_type}.tga"

        if tachyon_dat.exists():
            tachyon_cmd = "/usr/local/lib/vmd/tachyon_LINUXAMD64"
            result = subprocess.run(
                [tachyon_cmd, "-aasamples", "12", str(tachyon_dat),
                 "-format", "TGA", "-res", "1024", "768", "-o", str(tachyon_tga)],
                cwd=str(work_dir),
                capture_output=True,
                text=True,
                timeout=120
            )

            # 转换TGA到PNG
            if tachyon_tga.exists():
                try:
                    from PIL import Image as PILImage
                    img = PILImage.open(str(tachyon_tga))
                    img.save(str(image_file))
                except Exception as e:
                    logger.warning(f"PIL convert failed: {e}, trying ImageMagick")
                    subprocess.run(
                        ["convert", str(tachyon_tga), str(image_file)],
                        cwd=str(work_dir),
                        capture_output=True,
                        timeout=60
                    )

        if image_file.exists() and image_file.stat().st_size > 0:
            logger.info(f"Generated {orbital_type} image using VMD: {image_file}")
            return str(image_file)

        return None

    except Exception as e:
        logger.error(f"VMD {orbital_type} visualization failed: {e}")
        return None


def update_molecule_cache(db: "Session", job: "QCJob", result: "QCResult"):
    """更新分子QC缓存"""
    _, _, _, MoleculeQCCache, _ = _get_db_models()
    cache = db.query(MoleculeQCCache).filter(
        MoleculeQCCache.smiles == job.smiles
    ).first()

    homo_ev = result.homo * HARTREE_TO_EV if result.homo else None
    lumo_ev = result.lumo * HARTREE_TO_EV if result.lumo else None
    gap_ev = (lumo_ev - homo_ev) if (homo_ev and lumo_ev) else None

    if cache:
        # 更新现有缓存
        cache.molecule_name = job.molecule_name
        cache.basis_set = job.basis_set
        cache.functional = job.functional
        cache.energy_au = result.energy_au
        cache.homo_ev = homo_ev
        cache.lumo_ev = lumo_ev
        cache.homo_lumo_gap_ev = gap_ev
        cache.esp_min_kcal = result.esp_min_kcal
        cache.esp_max_kcal = result.esp_max_kcal
        cache.esp_image_path = result.esp_image_path
        cache.homo_image_path = result.homo_image_path
        cache.lumo_image_path = result.lumo_image_path
        cache.preferred_qc_result_id = result.id
        cache.calculation_count += 1
        cache.updated_at = datetime.now()
    else:
        # 创建新缓存
        cache = MoleculeQCCache(
            smiles=job.smiles,
            molecule_name=job.molecule_name,
            basis_set=job.basis_set,
            functional=job.functional,
            energy_au=result.energy_au,
            homo_ev=homo_ev,
            lumo_ev=lumo_ev,
            homo_lumo_gap_ev=gap_ev,
            esp_min_kcal=result.esp_min_kcal,
            esp_max_kcal=result.esp_max_kcal,
            esp_image_path=result.esp_image_path,
            homo_image_path=result.homo_image_path,
            lumo_image_path=result.lumo_image_path,
            preferred_qc_result_id=result.id,
            calculation_count=1
        )
        db.add(cache)


def postprocess_qc_job(self, job_id: int) -> Dict[str, Any]:
    """
    QC任务后处理

    提取计算结果、生成ESP可视化、更新缓存
    """
    db = self.db

    try:
        logger.info(f"[Task {self.request.id}] Starting postprocessing for QC job {job_id}")

        job = db.query(QCJob).filter(QCJob.id == job_id).first()
        if not job:
            return {"success": False, "error": f"QC Job {job_id} not found"}

        work_dir = Path(job.work_dir)
        if not work_dir.exists():
            job.status = QCJobStatus.FAILED
            job.error_message = "Work directory not found"
            db.commit()
            return {"success": False, "error": "Work directory not found"}

        # 1. 查找输出文件 - 智能搜索
        config = job.config or {}
        safe_name = config.get("safe_name")

        # 尝试多种可能的文件名
        possible_names = []
        if safe_name:
            possible_names.append(safe_name)
        possible_names.append(job.molecule_name)

        log_file = None
        fchk_file = None
        actual_name = None

        for name in possible_names:
            test_log = work_dir / f"{name}_out.log"
            if test_log.exists():
                log_file = test_log
                fchk_file = work_dir / f"{name}.fchk"
                actual_name = name
                break

        # 如果还没找到，在目录中搜索 *_out.log 文件
        if not log_file:
            import glob
            log_files = list(work_dir.glob("*_out.log"))
            # 排除 qc_out.log（这是Slurm的输出文件，不是Gaussian的）
            log_files = [f for f in log_files if f.name not in ("qc_out.log", "qc_err.log")]
            if log_files:
                log_file = log_files[0]
                # 从文件名提取实际名称
                actual_name = log_file.stem.replace("_out", "")
                fchk_file = work_dir / f"{actual_name}.fchk"
                logger.info(f"Found output file by glob: {log_file}")

        if not log_file or not log_file.exists():
            # 检查是否可以重试
            from app.tasks.qc_retry_helpers import can_retry_postprocess, find_output_files, extract_partial_results, auto_resubmit_failed_job
            
            # 尝试从工作目录自动查找输出文件
            engine = job.qc_engine or 'gaussian'
            found_files = find_output_files(work_dir, engine)
            
            if found_files.get('log_file'):
                log_file = found_files['log_file']
                fchk_file = found_files.get('fchk_file')
                logger.info(f"Auto-recovered output files for job {job_id}: {log_file}")
            else:
                # 尝试提取部分结果
                partial_results = extract_partial_results(work_dir, engine)
                if partial_results:
                    logger.warning(f"Extracted partial results for job {job_id}, marking as COMPLETED with warning")
                    gaussian_results = partial_results
                    # 继续后续处理,但不生成可视化
                    fchk_file = None
                elif can_retry_postprocess(job):
                    # 自动重新提交
                    if auto_resubmit_failed_job(job, db):
                        return {"success": False, "error": "Output file not found, job resubmitted", "resubmitted": True}
                    else:
                        job.status = QCJobStatus.FAILED
                        job.error_message = f"Output file not found after {job.retry_count} retries"
                        db.commit()
                        return {"success": False, "error": "Output file not found, resubmit failed"}
                else:
                    job.status = QCJobStatus.FAILED
                    job.error_message = f"Output file not found, max retries ({job.max_retries}) exceeded"
                    db.commit()
                    return {"success": False, "error": f"Output file not found after max retries"}  

        # 2. 提取结果(根据引擎类型)
        engine = job.qc_engine or 'gaussian'
        
        if engine == 'gaussian':
            gaussian_results = extract_gaussian_results(str(log_file))
            logger.info(f"Extracted Gaussian results: {gaussian_results}")
            
            # 检查结果是否有效
            if not gaussian_results.get('energy_au'):
                from app.tasks.qc_retry_helpers import can_retry_postprocess, auto_resubmit_failed_job, extract_partial_results
                
                # 尝试提取部分结果
                partial_results = extract_partial_results(work_dir, engine)
                if partial_results and partial_results.get('energy_au'):
                    logger.warning(f"Using partial results for job {job_id}")
                    gaussian_results = partial_results
                elif can_retry_postprocess(job):
                    if auto_resubmit_failed_job(job, db):
                        return {"success": False, "error": "Failed to extract results, job resubmitted", "resubmitted": True}
                
                if not gaussian_results.get('energy_au'):
                    job.status = QCJobStatus.FAILED
                    job.error_message = "Failed to extract energy from output"
                    job.retry_count += 1
                    db.commit()
                    return {"success": False, "error": "Failed to extract results"}
        
        elif engine == 'pyscf':
            from app.services.qc_engines.pyscf import PySCFEngine
            pyscf_engine = PySCFEngine()
            result = pyscf_engine.parse_output(work_dir / "pyscf_out.log")
            
            if not result.success:
                from app.tasks.qc_retry_helpers import can_retry_postprocess, auto_resubmit_failed_job
                if can_retry_postprocess(job):
                    if auto_resubmit_failed_job(job, db):
                        return {"success": False, "error": result.error_message, "resubmitted": True}
                
                job.status = QCJobStatus.FAILED
                job.error_message = result.error_message
                job.retry_count += 1
                db.commit()
                return {"success": False, "error": result.error_message}
            
            # 转换PySCF结果为Gaussian格式
            gaussian_results = {
                'energy_au': result.energy_au,
                'homo': result.homo / 27.2114 if result.homo else None,  # eV to Hartree
                'lumo': result.lumo / 27.2114 if result.lumo else None,
            }
        else:
            job.status = QCJobStatus.FAILED
            job.error_message = f"Unknown QC engine: {engine}"
            db.commit()
            return {"success": False, "error": f"Unknown engine: {engine}"}

        # 3. 生成可视化图像（如果fchk存在）
        esp_image_path = None
        homo_image_path = None
        lumo_image_path = None
        esp_min, esp_max = None, None

        if fchk_file and fchk_file.exists():
            # 使用实际找到的文件名或原始分子名称
            esp_name = actual_name or job.molecule_name

            # 3.1 生成ESP可视化
            esp_image_path = generate_esp_visualization(work_dir, esp_name, str(fchk_file))
            logger.info(f"ESP visualization: {esp_image_path}")

            # 提取ESP值
            surfanalysis_file = work_dir / "surfanalysis.txt"
            esp_min, esp_max = extract_esp_values(str(surfanalysis_file))
            logger.info(f"ESP values: min={esp_min}, max={esp_max}")

            # 3.2 生成HOMO轨道图像
            homo_image_path = generate_orbital_visualization(work_dir, str(fchk_file), "HOMO")
            logger.info(f"HOMO visualization: {homo_image_path}")

            # 3.3 生成LUMO轨道图像
            lumo_image_path = generate_orbital_visualization(work_dir, str(fchk_file), "LUMO")
            logger.info(f"LUMO visualization: {lumo_image_path}")

        # 4. 计算HOMO-LUMO gap
        homo_lumo_gap = None
        if gaussian_results["homo"] and gaussian_results["lumo"]:
            homo_lumo_gap = (gaussian_results["lumo"] - gaussian_results["homo"]) * HARTREE_TO_EV

        # 5. 创建QC结果记录
        qc_result = QCResult(
            qc_job_id=job.id,
            smiles=job.smiles,
            energy_au=gaussian_results["energy_au"],
            homo=gaussian_results["homo"],
            lumo=gaussian_results["lumo"],
            homo_lumo_gap=homo_lumo_gap,
            esp_min_kcal=esp_min,
            esp_max_kcal=esp_max,
            esp_image_path=esp_image_path,
            homo_image_path=homo_image_path,
            lumo_image_path=lumo_image_path,
            fchk_file_path=str(fchk_file) if fchk_file.exists() else None,
            log_file_path=str(log_file),
        )
        db.add(qc_result)
        db.flush()  # 获取result.id

        # 6. 更新分子缓存
        update_molecule_cache(db, job, qc_result)

        # 7. 获取 CPU 核时数据
        if job.slurm_job_id:
            cpu_hours = get_slurm_job_cpu_hours(job.slurm_job_id)
            if cpu_hours > 0:
                job.actual_cpu_hours = cpu_hours
                logger.info(f"[Task {self.request.id}] QC job {job_id} CPU hours: {cpu_hours:.2f}h")

        # 8. 更新任务状态
        job.status = QCJobStatus.COMPLETED
        job.progress = 100.0
        job.finished_at = datetime.now()
        job.log_file = str(log_file)

        # 9. 自动设置为公开（QC数据必须公开）
        job.visibility = "PUBLIC"

        db.commit()

        logger.info(f"[Task {self.request.id}] QC job {job_id} postprocessing completed")

        return {
            "success": True,
            "job_id": job_id,
            "result_id": qc_result.id,
            "energy_au": gaussian_results["energy_au"],
            "homo": gaussian_results["homo"],
            "lumo": gaussian_results["lumo"],
            "esp_min": esp_min,
            "esp_max": esp_max,
        }

    except Exception as exc:
        logger.exception(f"[Task {self.request.id}] Postprocessing failed: {exc}")

        try:
            job = db.query(QCJob).filter(QCJob.id == job_id).first()
            if job:
                job.status = QCJobStatus.FAILED
                job.error_message = f"Postprocessing failed: {str(exc)}"
                db.commit()
        except Exception as db_exc:
            logger.error(f"Failed to update job status: {db_exc}")

        raise self.retry(exc=exc, countdown=120)


# 注册 Celery 任务（延迟注册，避免在导入时出错）
try:
    celery_app = _get_celery_app()
    postprocess_qc_job = celery_app.task(
        bind=True,
        base=DatabaseTask,
        name="app.tasks.qc_postprocess.postprocess_qc_job",
        max_retries=2,
        default_retry_delay=120,
    )(postprocess_qc_job)
except Exception:
    # 如果 celery_app 不可用（例如在 polling_worker 中），跳过注册
    pass
