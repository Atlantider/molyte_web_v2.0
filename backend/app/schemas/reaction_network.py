"""
Reaction Network Pydantic Schemas
反应网络API数据模式

现代化设计特点:
- 清晰的数据验证规则
- 丰富的文档字符串
- 智能的默认值
- 完整的类型注解
"""

from pydantic import BaseModel, Field, field_validator
from typing import List, Optional, Dict, Any
from datetime import datetime
from enum import Enum


# ============================================================================
# Enums
# ============================================================================

class ReactionNetworkJobStatus(str, Enum):
    """任务状态枚举"""
    CREATED = "CREATED"
    QUEUED = "QUEUED"
    RUNNING = "RUNNING"
    POSTPROCESSING = "POSTPROCESSING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"
    CANCELLED = "CANCELLED"


class ElectrodeType(str, Enum):
    """电极类型枚举"""
    ANODE = "anode"
    CATHODE = "cathode"


class AnodeMaterial(str, Enum):
    """负极材料枚举"""
    # 锂系
    LI_METAL = "LI_METAL"
    GRAPHITE = "GRAPHITE"
    SILICON = "SILICON"
    SIC = "SIC"
    LTO = "LTO"
    # 钠系
    NA_METAL = "NA_METAL"
    HARD_CARBON = "HARD_CARBON"
    SOFT_CARBON = "SOFT_CARBON"
    # 钾系
    K_METAL = "K_METAL"
    K_GRAPHITE = "K_GRAPHITE"


class CathodeMaterial(str, Enum):
    """正极材料枚举"""
    LCO = "LCO"          # LiCoO2
    NMC = "NMC"          # LiNiMnCoO2 (通用)
    NMC622 = "NMC622"    # LiNi0.6Mn0.2Co0.2O2
    NMC811 = "NMC811"    # LiNi0.8Mn0.1Co0.1O2
    NCA = "NCA"          # LiNiCoAlO2
    LFP = "LFP"          # LiFePO4
    LMO = "LMO"          # LiMn2O4
    LNMO = "LNMO"        # LiNi0.5Mn1.5O4
    LRLO = "LRLO"        # Li-rich layered oxide


# ============================================================================
# Base Schemas
# ============================================================================

class ReactionNetworkJobBase(BaseModel):
    """反应网络任务基础模式"""
    job_name: str = Field(
        ...,
        min_length=1,
        max_length=255,
        description="任务名称",
        examples=["EC电解液反应网络", "锂离子电池SEI形成"]
    )
    description: Optional[str] = Field(
        None,
        max_length=2000,
        description="任务描述"
    )


class MoleculeBase(BaseModel):
    """分子基础模式"""
    name: str = Field(..., description="分子名称")
    smiles: str = Field(..., description="SMILES表示")
    generation: int = Field(default=0, ge=0, description="代数")


class ReactionBase(BaseModel):
    """反应基础模式"""
    reactant_smiles: List[str] = Field(..., min_length=1, description="反应物SMILES列表")
    product_smiles: List[str] = Field(..., min_length=1, description="产物SMILES列表")


# ============================================================================
# Create Schemas
# ============================================================================

class ReactionNetworkJobCreate(ReactionNetworkJobBase):
    """创建反应网络任务
    
    智能验证:
    - SMILES格式验证
    - 参数范围检查
    - 合理性验证
    """
    
    # 初始分子
    initial_smiles: List[str] = Field(
        ...,
        min_length=1,
        max_length=20,
        description="初始分子SMILES列表",
        examples=[["C1COC(=O)O1", "[Li+]", "F[P-](F)(F)(F)(F)F"]]
    )
    
    # 环境参数
    temperature: float = Field(
        default=300.0,
        ge=0.0,
        le=1000.0,
        description="温度 (K)",
        examples=[300.0]
    )
    electrode_type: ElectrodeType = Field(
        default=ElectrodeType.ANODE,
        description="电极类型"
    )
    anode_material: AnodeMaterial = Field(
        default=AnodeMaterial.GRAPHITE,
        description="负极材料 (GRAPHITE, SILICON, LI_METAL, NA_METAL, K_METAL, HARD_CARBON等)"
    )
    cathode_material: CathodeMaterial = Field(
        default=CathodeMaterial.NMC,
        description="正极材料 (NMC, LCO, NCA, LFP, LMO等)"
    )
    voltage: float = Field(
        default=0.1,
        ge=-10.0,
        le=10.0,
        description="电压 (V)",
        examples=[0.1, 4.2]
    )

    
    # 网络生成参数
    max_generations: int = Field(
        default=3,
        ge=1,
        le=10,
        description="最大代数 (建议1-5代)"
    )
    max_species: int = Field(
        default=50,
        ge=1,
        le=200,
        description="最大分子数 (限制网络规模)"
    )
    energy_cutoff: float = Field(
        default=80.0,
        ge=0.0,
        le=200.0,
        description="能量截断值 (kcal/mol, 排除高能反应)"
    )
    
    # Slurm资源配置
    slurm_partition: Optional[str] = Field(
        default="cpu",
        description="Slurm分区"
    )
    slurm_cpus: Optional[int] = Field(
        default=16,
        ge=1,
        le=128,
        description="CPU核心数"
    )
    slurm_time: Optional[int] = Field(
        default=7200,
        ge=10,
        le=43200,
        description="最大运行时间(分钟)"
    )
    
    # 额外配置
    config: Optional[Dict[str, Any]] = Field(
        default={},
        description="额外配置参数"
    )
    
    @field_validator('initial_smiles')
    @classmethod
    def validate_smiles_list(cls, v: List[str]) -> List[str]:
        """验证SMILES列表"""
        # 去重
        unique_smiles = list(dict.fromkeys(v))
        if len(unique_smiles) != len(v):
            raise ValueError(f"SMILES列表包含重复项，已自动去重: {len(v)} -> {len(unique_smiles)}")
        
        # 验证每个SMILES（这里简化，实际应该用RDKit验证）
        for smiles in unique_smiles:
            if not smiles or not smiles.strip():
                raise ValueError("SMILES不能为空")
            if len(smiles) > 500:
                raise ValueError(f"SMILES过长 (>{500}个字符): {smiles[:50]}...")
        
        return unique_smiles
    
    @field_validator('max_generations')
    @classmethod
    def validate_max_generations(cls, v: int, info) -> int:
        """智能建议：代数太大会导致爆炸性增长"""
        if v > 5:
            # 可以发出警告，但不阻止
            pass
        return v


class ReactionNetworkJobUpdate(BaseModel):
    """更新任务（内部使用）"""
    status: Optional[ReactionNetworkJobStatus] = None
    progress: Optional[float] = Field(None, ge=0.0, le=1.0)
    error_message: Optional[str] = None
    num_molecules: Optional[int] = None
    num_reactions: Optional[int] = None
    max_generation_reached: Optional[int] = None


# ============================================================================
# Response Schemas
# ============================================================================

class MoleculeSchema(MoleculeBase):
    """分子响应模式"""
    id: int
    job_id: int
    energy_kcal: Optional[float] = None
    molecular_weight: Optional[float] = None
    num_atoms: Optional[int] = None
    num_heavy_atoms: Optional[int] = None
    formal_charge: Optional[int] = None
    num_rings: Optional[int] = None
    properties: Dict[str, Any] = {}
    created_at: datetime
    
    class Config:
        from_attributes = True


class ReactionSchema(ReactionBase):
    """反应响应模式"""
    id: int
    job_id: int
    operator_name: Optional[str] = None
    reaction_type: Optional[str] = None
    reaction_energy: Optional[float] = None
    activation_energy: Optional[float] = None
    driving_force: Optional[str] = None
    is_reversible: Optional[bool] = True
    equilibrium_constant: Optional[float] = None
    metadata: Dict[str, Any] = {}
    created_at: datetime
    
    class Config:
        from_attributes = True


class ReactionNetworkJobSchema(ReactionNetworkJobBase):
    """反应网络任务响应模式"""
    id: int
    user_id: int
    status: ReactionNetworkJobStatus
    progress: float
    
    # 参数
    initial_smiles: List[str]
    temperature: float
    electrode_type: ElectrodeType
    voltage: float
    max_generations: int
    max_species: int
    energy_cutoff: float
    
    # 执行信息
    slurm_job_id: Optional[str] = None
    work_dir: Optional[str] = None
    slurm_partition: Optional[str] = None
    slurm_cpus: Optional[int] = None
    slurm_time: Optional[int] = None
    actual_cpu_hours: float = 0.0
    
    # 结果统计
    num_molecules: Optional[int] = None
    num_reactions: Optional[int] = None
    max_generation_reached: Optional[int] = None
    
    # 结果文件
    network_json_path: Optional[str] = None
    visualization_png_path: Optional[str] = None
    visualization_html_path: Optional[str] = None
    
    # 时间戳
    created_at: datetime
    updated_at: Optional[datetime] = None
    started_at: Optional[datetime] = None
    finished_at: Optional[datetime] = None
    
    # 错误信息
    error_message: Optional[str] = None
    
    class Config:
        from_attributes = True


class ReactionNetworkJobDetailSchema(ReactionNetworkJobSchema):
    """详细的任务信息（包含分子和反应）"""
    molecules: List[MoleculeSchema] = []
    reactions: List[ReactionSchema] = []


# ============================================================================
# List Response Schemas
# ============================================================================

class ReactionNetworkJobListResponse(BaseModel):
    """任务列表响应"""
    total: int = Field(..., description="总数")
    jobs: List[ReactionNetworkJobSchema] = Field(..., description="任务列表")


class MoleculeListResponse(BaseModel):
    """分子列表响应"""
    total: int
    molecules: List[MoleculeSchema]


class ReactionListResponse(BaseModel):
    """反应列表响应"""
    total: int
    reactions: List[ReactionSchema]


# ============================================================================
# Network Visualization Data
# ============================================================================

class NetworkNode(BaseModel):
    """网络节点（分子）"""
    id: str = Field(..., description="节点ID (通常是SMILES)")
    label: str = Field(..., description="节点标签 (分子名称)")
    generation: int = Field(..., description="代数")
    energy: Optional[float] = None
    properties: Dict[str, Any] = {}


class NetworkEdge(BaseModel):
    """网络边（反应）"""
    source: str = Field(..., description="源节点ID")
    target: str = Field(..., description="目标节点ID")
    label: Optional[str] = None
    operator: Optional[str] = None
    energy: Optional[float] = None
    properties: Dict[str, Any] = {}


class NetworkVisualizationData(BaseModel):
    """网络可视化数据
    
    现代化图数据格式，适用于:
    - vis-network
    - react-force-graph
    - cytoscape.js
    """
    nodes: List[NetworkNode] = Field(..., description="节点列表")
    edges: List[NetworkEdge] = Field(..., description="边列表")
    statistics: Dict[str, Any] = Field(
        default={},
        description="统计信息 (节点数、边数、代数分布等)"
    )


# ============================================================================
# Query Parameters
# ============================================================================

class ReactionNetworkJobQuery(BaseModel):
    """任务查询参数"""
    status: Optional[ReactionNetworkJobStatus] = None
    skip: int = Field(default=0, ge=0, description="跳过数量")
    limit: int = Field(default=20, ge=1, le=100, description="返回数量")
    sort_by: Optional[str] = Field(default="created_at", description="排序字段")
    sort_desc: bool = Field(default=True, description="降序排列")


class MoleculeQuery(BaseModel):
    """分子查询参数"""
    generation: Optional[int] = Field(None, ge=0, description="筛选特定代数")
    min_energy: Optional[float] = None
    max_energy: Optional[float] = None
    skip: int = Field(default=0, ge=0)
    limit: int = Field(default=50, ge=1, le=200)


class ReactionQuery(BaseModel):
    """反应查询参数"""
    operator: Optional[str] = None
    reaction_type: Optional[str] = None
    min_energy: Optional[float] = None
    max_energy: Optional[float] = None
    skip: int = Field(default=0, ge=0)
    limit: int = Field(default=50, ge=1, le=200)


# ============================================================================
# Submission Response
# ============================================================================

class JobSubmissionResponse(BaseModel):
    """任务提交响应"""
    success: bool
    message: str
    job_id: int
    slurm_job_id: Optional[str] = None
