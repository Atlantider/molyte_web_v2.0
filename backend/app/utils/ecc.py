"""
ECC (Electronic Continuum Correction) 电荷缩放工具

用于对离子的原子电荷应用ECC缩放，以补偿经典力场中缺少的电子极化效应。
"""
import re
import logging
from pathlib import Path
from typing import Optional, Tuple


logger = logging.getLogger(__name__)


def apply_ecc_to_lt_content(lt_content: str, ecc_factor: float = 0.8) -> str:
    """
    对LT文件内容应用ECC电荷缩放
    
    LT文件中电荷有两个位置需要缩放：
    1. Data Atoms 中的电荷列（第4列，在mol:id和@atom:type之后）
    2. In Charges 中的 set type @atom:X charge Y
    
    Args:
        lt_content: LT文件内容
        ecc_factor: ECC缩放因子（默认0.8）
        
    Returns:
        缩放后的LT文件内容
    """
    if not (0 < ecc_factor <= 1.0):
        raise ValueError(f"ECC factor must be between 0 and 1, got {ecc_factor}")
    
    # 1. 缩放 "In Charges" 块中的电荷
    # 匹配: set type @atom:Li charge 1
    def scale_in_charges(match):
        prefix = match.group(1)  # "set type @atom:Li  charge "
        charge_str = match.group(2)
        try:
            charge = float(charge_str)
            scaled = charge * ecc_factor
            # 保持合理的精度
            if abs(scaled) < 0.0001:
                return f"{prefix}0"
            elif scaled == int(scaled):
                return f"{prefix}{int(scaled)}"
            else:
                return f"{prefix}{scaled:.4f}".rstrip('0').rstrip('.')
        except ValueError:
            return match.group(0)
    
    pattern_in_charges = r'(set\s+type\s+@atom:\w+\s+charge\s+)([-\d.]+)'
    scaled_content = re.sub(pattern_in_charges, scale_in_charges, lt_content, flags=re.IGNORECASE)
    
    # 2. 缩放 "Data Atoms" 块中的电荷
    # 格式: $atom:S1   $mol @atom:S   1.02   -1.491799951   -0.782199979    0.004400000
    # 电荷是第4个数值（在 @atom:type 之后的第一个数字）
    def scale_data_atoms_line(line: str) -> str:
        # 检查是否是 Data Atoms 行
        if '$atom:' not in line or '@atom:' not in line:
            return line
        
        # 分割成 tokens
        parts = line.split()
        if len(parts) < 5:
            return line
        
        # 找到 @atom:type 的位置
        atom_type_idx = -1
        for i, p in enumerate(parts):
            if p.startswith('@atom:'):
                atom_type_idx = i
                break
        
        if atom_type_idx == -1 or atom_type_idx + 1 >= len(parts):
            return line
        
        # @atom:type 后的第一个数字是电荷
        charge_idx = atom_type_idx + 1
        try:
            original_charge = float(parts[charge_idx])
            scaled_charge = original_charge * ecc_factor
            
            # 格式化缩放后的电荷
            if abs(scaled_charge) < 0.0001:
                parts[charge_idx] = "0"
            elif scaled_charge == int(scaled_charge):
                parts[charge_idx] = str(int(scaled_charge))
            else:
                parts[charge_idx] = f"{scaled_charge:.4f}".rstrip('0').rstrip('.')
            
            # 重建行，保持原始格式尽可能一致
            return '    ' + '  '.join(parts)
        except (ValueError, IndexError):
            return line
    
    # 处理 Data Atoms 块
    lines = scaled_content.split('\n')
    in_data_atoms = False
    result_lines = []
    
    for line in lines:
        if 'write("Data Atoms")' in line or "write('Data Atoms')" in line:
            in_data_atoms = True
            result_lines.append(line)
        elif in_data_atoms and line.strip() == '}':
            in_data_atoms = False
            result_lines.append(line)
        elif in_data_atoms and '$atom:' in line:
            result_lines.append(scale_data_atoms_line(line))
        else:
            result_lines.append(line)
    
    return '\n'.join(result_lines)


def apply_ecc_to_lt_file(
    input_path: Path,
    output_path: Optional[Path] = None,
    ecc_factor: float = 0.8
) -> Path:
    """
    对LT文件应用ECC缩放并保存
    
    Args:
        input_path: 输入LT文件路径
        output_path: 输出LT文件路径（如果为None，则覆盖原文件）
        ecc_factor: ECC缩放因子
        
    Returns:
        输出文件路径
    """
    with open(input_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    scaled_content = apply_ecc_to_lt_content(content, ecc_factor)
    
    if output_path is None:
        output_path = input_path
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(scaled_content)
    
    logger.info(f"Applied ECC ({ecc_factor}) to {input_path.name} -> {output_path.name}")
    
    return output_path


def get_total_charge_from_lt(lt_content: str) -> float:
    """
    从LT文件内容中计算总电荷
    
    Args:
        lt_content: LT文件内容
        
    Returns:
        总电荷
    """
    total = 0.0
    
    # 从 Data Atoms 中提取电荷
    lines = lt_content.split('\n')
    in_data_atoms = False
    
    for line in lines:
        if 'write("Data Atoms")' in line or "write('Data Atoms')" in line:
            in_data_atoms = True
            continue
        elif in_data_atoms and line.strip() == '}':
            break
        elif in_data_atoms and '$atom:' in line:
            parts = line.split()
            for i, p in enumerate(parts):
                if p.startswith('@atom:'):
                    if i + 1 < len(parts):
                        try:
                            total += float(parts[i + 1])
                        except ValueError:
                            pass
                    break
    
    return round(total, 4)


def validate_ecc_scaling(
    original_content: str,
    scaled_content: str,
    ecc_factor: float
) -> Tuple[bool, str]:
    """
    验证ECC缩放是否正确
    
    Args:
        original_content: 原始LT内容
        scaled_content: 缩放后的LT内容
        ecc_factor: 使用的缩放因子
        
    Returns:
        (是否正确, 验证信息)
    """
    original_charge = get_total_charge_from_lt(original_content)
    scaled_charge = get_total_charge_from_lt(scaled_content)
    expected_charge = round(original_charge * ecc_factor, 4)
    
    if abs(scaled_charge - expected_charge) < 0.01:
        return True, f"OK: {original_charge} * {ecc_factor} = {expected_charge} (actual: {scaled_charge})"
    else:
        return False, f"MISMATCH: {original_charge} * {ecc_factor} = {expected_charge} (actual: {scaled_charge})"
