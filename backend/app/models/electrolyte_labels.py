"""
电解液标签分类系统

用于给电解液配方添加机器学习友好的标签
"""
from enum import Enum
from typing import List, Optional
from pydantic import BaseModel


class BatteryType(str, Enum):
    """电池类型 - 单选"""
    LITHIUM_ION = "lithium_ion"          # 锂离子
    LITHIUM_METAL = "lithium_metal"      # 锂金属
    SODIUM_ION = "sodium_ion"            # 钠离子
    POTASSIUM_ION = "potassium_ion"      # 钾离子
    ZINC_ION = "zinc_ion"                # 锌离子
    MAGNESIUM_ION = "magnesium_ion"      # 镁离子
    OTHER = "other"                       # 其他


class AnodeType(str, Enum):
    """负极材料 - 多选"""
    GRAPHITE = "graphite"        # 石墨
    SILICON = "silicon"          # 硅/硅碳
    LI_METAL = "li_metal"        # 锂金属
    LTO = "lto"                  # 钛酸锂 Li4Ti5O12
    HARD_CARBON = "hard_carbon"  # 硬碳
    ZINC = "zinc"                # 锌金属


class CathodeType(str, Enum):
    """正极材料 - 多选"""
    NCM = "ncm"                  # 三元材料 NCM
    NCA = "nca"                  # 镍钴铝 NCA
    LFP = "lfp"                  # 磷酸铁锂
    LCO = "lco"                  # 钴酸锂
    LMO = "lmo"                  # 锰酸锂
    HIGH_NICKEL = "high_nickel"  # 高镍 (NCM811等)
    LI_RICH = "li_rich"          # 富锂锰基
    SULFUR = "sulfur"            # 硫 (锂硫电池)
    AIR = "air"                  # 空气 (锂空电池)


class OperatingCondition(str, Enum):
    """工作条件 - 多选"""
    STANDARD = "standard"                      # 常规条件
    HIGH_TEMP = "high_temp"                    # 高温 (>45°C)
    LOW_TEMP = "low_temp"                      # 低温 (<0°C)
    ULTRA_LOW_TEMP = "ultra_low_temp"          # 超低温 (<-20°C)
    WIDE_TEMP = "wide_temp"                    # 宽温域 (-20~60°C)
    HIGH_VOLTAGE = "high_voltage"              # 高压 (>4.3V)
    ULTRA_HIGH_VOLTAGE = "ultra_high_voltage"  # 超高压 (>4.5V)
    FAST_CHARGING = "fast_charging"            # 快充


class ElectrolyteFormType(str, Enum):
    """电解液类型/形态 - 单选"""
    ORGANIC_LIQUID = "organic_liquid"          # 有机液态
    AQUEOUS = "aqueous"                        # 水系
    IONIC_LIQUID = "ionic_liquid"              # 离子液体基
    HIGH_CONCENTRATION = "high_concentration"  # 高浓度盐 (>3M)
    LHCE = "lhce"                              # 局域化高浓电解液


class ElectrolyteLabels(BaseModel):
    """电解液标签集合"""
    battery_type: Optional[BatteryType] = None
    anode_types: List[AnodeType] = []
    cathode_types: List[CathodeType] = []
    conditions: List[OperatingCondition] = []
    electrolyte_type: Optional[ElectrolyteFormType] = None
    
    class Config:
        use_enum_values = True


# 中文显示名称映射
BATTERY_TYPE_NAMES = {
    BatteryType.LITHIUM_ION: "锂离子电池",
    BatteryType.LITHIUM_METAL: "锂金属电池",
    BatteryType.SODIUM_ION: "钠离子电池",
    BatteryType.POTASSIUM_ION: "钾离子电池",
    BatteryType.ZINC_ION: "锌离子电池",
    BatteryType.MAGNESIUM_ION: "镁离子电池",
    BatteryType.OTHER: "其他",
}

ANODE_TYPE_NAMES = {
    AnodeType.GRAPHITE: "石墨",
    AnodeType.SILICON: "硅/硅碳",
    AnodeType.LI_METAL: "锂金属",
    AnodeType.LTO: "钛酸锂",
    AnodeType.HARD_CARBON: "硬碳",
    AnodeType.ZINC: "锌金属",
}

CATHODE_TYPE_NAMES = {
    CathodeType.NCM: "NCM三元",
    CathodeType.NCA: "NCA",
    CathodeType.LFP: "磷酸铁锂",
    CathodeType.LCO: "钴酸锂",
    CathodeType.LMO: "锰酸锂",
    CathodeType.HIGH_NICKEL: "高镍",
    CathodeType.LI_RICH: "富锂锰基",
    CathodeType.SULFUR: "硫正极",
    CathodeType.AIR: "空气正极",
}

CONDITION_NAMES = {
    OperatingCondition.STANDARD: "常规条件",
    OperatingCondition.HIGH_TEMP: "高温 (>45°C)",
    OperatingCondition.LOW_TEMP: "低温 (<0°C)",
    OperatingCondition.ULTRA_LOW_TEMP: "超低温 (<-20°C)",
    OperatingCondition.WIDE_TEMP: "宽温域",
    OperatingCondition.HIGH_VOLTAGE: "高压 (>4.3V)",
    OperatingCondition.ULTRA_HIGH_VOLTAGE: "超高压 (>4.5V)",
    OperatingCondition.FAST_CHARGING: "快充",
}

ELECTROLYTE_TYPE_NAMES = {
    ElectrolyteFormType.ORGANIC_LIQUID: "有机液态",
    ElectrolyteFormType.AQUEOUS: "水系",
    ElectrolyteFormType.IONIC_LIQUID: "离子液体",
    ElectrolyteFormType.HIGH_CONCENTRATION: "高浓度电解液",
    ElectrolyteFormType.LHCE: "LHCE局域化高浓",
}


def get_label_display_name(category: str, value: str) -> str:
    """获取标签的中文显示名称"""
    mappings = {
        "battery_type": BATTERY_TYPE_NAMES,
        "anode_types": ANODE_TYPE_NAMES,
        "cathode_types": CATHODE_TYPE_NAMES,
        "conditions": CONDITION_NAMES,
        "electrolyte_type": ELECTROLYTE_TYPE_NAMES,
    }
    
    if category not in mappings:
        return value
    
    # 尝试从枚举映射获取
    for enum_val, name in mappings[category].items():
        if enum_val.value == value:
            return name
    
    return value


def labels_to_display_string(labels: dict) -> str:
    """将标签字典转换为显示字符串"""
    parts = []
    
    if labels.get("battery_type"):
        parts.append(get_label_display_name("battery_type", labels["battery_type"]))
    
    if labels.get("conditions"):
        for cond in labels["conditions"]:
            parts.append(get_label_display_name("conditions", cond))
    
    if labels.get("electrolyte_type"):
        parts.append(get_label_display_name("electrolyte_type", labels["electrolyte_type"]))
    
    return " | ".join(parts) if parts else "未分类"
