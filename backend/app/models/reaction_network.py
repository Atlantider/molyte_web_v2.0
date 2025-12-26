"""
Reaction Network Job Models
反应网络任务数据模型

智能化设计特点:
- 完整的任务生命周期管理
- 丰富的元数据支持
- 与现有模块一致的设计模式
- 支持软删除和可见性控制
"""

from sqlalchemy import Column, Integer, String, Float, DateTime, Text, Boolean, ForeignKey, Enum as SQLEnum
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship
from datetime import datetime
import enum

from app.database import Base


class ReactionNetworkJobStatus(str, enum.Enum):
    """反应网络任务状态"""
    CREATED = "CREATED"        # 已创建，待提交
    QUEUED = "QUEUED"          # 已提交到队列
    RUNNING = "RUNNING"        # 正在运行
    POSTPROCESSING = "POSTPROCESSING"  # 后处理中
    COMPLETED = "COMPLETED"    # 已完成
    FAILED = "FAILED"          # 失败
    CANCELLED = "CANCELLED"    # 已取消


class ElectrodeType(str, enum.Enum):
    """电极类型"""
    ANODE = "anode"      # 阳极（负极）
    CATHODE = "cathode"  # 阴极（正极）


class ReactionNetworkJob(Base):
    """反应网络生成任务
    
    智能特性:
    - 自动进度追踪
    - 完整的执行历史
    - 丰富的统计信息
    - 结果文件管理
    """
    __tablename__ = "reaction_network_jobs"
    
    # ============================================================================
    # 基础信息
    # ============================================================================
    id = Column(Integer, primary_key=True, index=True, comment="任务ID")
    user_id = Column(Integer, ForeignKey("users.id"), nullable=False, index=True, comment="用户ID")
    
    job_name = Column(String(255), nullable=False, comment="任务名称")
    description = Column(Text, comment="任务描述")
    
    # ============================================================================
    # 状态管理
    # ============================================================================
    status = Column(
        SQLEnum(ReactionNetworkJobStatus),
        default=ReactionNetworkJobStatus.CREATED,
        nullable=False,
        index=True,
        comment="任务状态"
    )
    progress = Column(Float, default=0.0, comment="进度 (0.0-1.0)")
    error_message = Column(Text, comment="错误信息")
    
    # ============================================================================
    # 输入参数
    # ============================================================================
    initial_smiles = Column(JSONB, nullable=False, comment="初始分子SMILES列表")
    
    # 环境参数
    temperature = Column(Float, default=300.0, comment="温度 (K)")
    electrode_type = Column(
        SQLEnum(ElectrodeType),
        default=ElectrodeType.ANODE,
        comment="电极类型"
    )
    voltage = Column(Float, default=0.1, comment="电压 (V)")
    
    # 网络生成参数
    max_generations = Column(Integer, default=3, comment="最大代数")
    max_species = Column(Integer, default=50, comment="最大分子数")
    energy_cutoff = Column(Float, default=80.0, comment="能量截断值 (kcal/mol)")
    
    # 额外配置
    config = Column(JSONB, default={}, comment="额外配置参数")
    
    # ============================================================================
    # 执行信息
    # ============================================================================
    slurm_job_id = Column(String(64), index=True, comment="Slurm任务ID")
    work_dir = Column(String(512), comment="工作目录")
    
    # Slurm资源配置
    slurm_partition = Column(String(64), default="cpu", comment="Slurm分区")
    slurm_cpus = Column(Integer, default=16, comment="CPU核心数")
    slurm_time = Column(Integer, default=7200, comment="最大运行时间(分钟)")
    
    # CPU核时统计
    actual_cpu_hours = Column(Float, default=0.0, comment="实际CPU核时")
    
    # ============================================================================
    # 结果统计
    # ============================================================================
    num_molecules = Column(Integer, comment="生成的分子数")
    num_reactions = Column(Integer, comment="发现的反应数")
    max_generation_reached = Column(Integer, comment="达到的最大代数")
    
    # 结果文件路径
    network_json_path = Column(String(512), comment="网络JSON文件路径")
    visualization_png_path = Column(String(512), comment="可视化PNG文件路径")
    visualization_html_path = Column(String(512), comment="交互式HTML文件路径")
    
    # ============================================================================
    # 时间戳
    # ============================================================================
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False, comment="创建时间")
    updated_at = Column(DateTime, onupdate=datetime.utcnow, comment="更新时间")
    started_at = Column(DateTime, comment="开始时间")
    finished_at = Column(DateTime, comment="结束时间")
    
    # ============================================================================
    # 可见性管理
    # ============================================================================
    is_deleted = Column(Boolean, default=False, comment="是否已删除")
    deleted_at = Column(DateTime, comment="删除时间")
    deleted_by = Column(Integer, ForeignKey("users.id"), comment="删除操作者")
    delete_reason = Column(String(512), comment="删除原因")
    
    # ============================================================================
    # 关系
    # ============================================================================
    user = relationship("User", foreign_keys=[user_id], back_populates="reaction_network_jobs")
    molecules = relationship(
        "ReactionNetworkMolecule",
        back_populates="job",
        cascade="all, delete-orphan"
    )
    reactions = relationship(
        "ReactionNetworkReaction",
        back_populates="job",
        cascade="all, delete-orphan"
    )
    
    def __repr__(self):
        return f"<ReactionNetworkJob(id={self.id}, name='{self.job_name}', status={self.status})>"


class ReactionNetworkMolecule(Base):
    """反应网络中的分子
    
    智能特性:
    - 完整的分子属性
    - 代际追踪
    - 能量信息
    """
    __tablename__ = "reaction_network_molecules"
    
    id = Column(Integer, primary_key=True, index=True, comment="分子ID")
    job_id = Column(
        Integer,
        ForeignKey("reaction_network_jobs.id", ondelete="CASCADE"),
        nullable=False,
        index=True,
        comment="任务ID"
    )
    
    # 分子基本信息
    name = Column(String(255), nullable=False, comment="分子名称")
    smiles = Column(String(512), nullable=False, index=True, comment="SMILES表示")
    generation = Column(Integer, default=0, index=True, comment="代数（0为初始分子）")
    
    # 能量信息
    energy_kcal = Column(Float, comment="XTB计算的能量 (kcal/mol)")
    
    # 分子属性
    molecular_weight = Column(Float, comment="分子量")
    num_atoms = Column(Integer, comment="原子数")
    num_heavy_atoms = Column(Integer, comment="重原子数")
    formal_charge = Column(Integer, default=0, comment="形式电荷")
    num_rings = Column(Integer, comment="环数")
    
    # 分子指纹（用于相似度搜索）
    morgan_fingerprint = Column(String(512), comment="Morgan指纹")
    
    # 额外属性
    properties = Column(JSONB, default={}, comment="其他属性")
    
    # 时间戳
    created_at = Column(DateTime, default=datetime.utcnow, comment="创建时间")
    
    # 关系
    job = relationship("ReactionNetworkJob", back_populates="molecules")
    
    def __repr__(self):
        return f"<Molecule(id={self.id}, name='{self.name}', gen={self.generation})>"


class ReactionNetworkReaction(Base):
    """反应网络中的反应
    
    智能特性:
    - 完整的反应描述
    - 能量学信息
    - 算符追踪
    """
    __tablename__ = "reaction_network_reactions"
    
    id = Column(Integer, primary_key=True, index=True, comment="反应ID")
    job_id = Column(
        Integer,
        ForeignKey("reaction_network_jobs.id", ondelete="CASCADE"),
        nullable=False,
        index=True,
        comment="任务ID"
    )
    
    # 反应物和产物
    reactant_smiles = Column(JSONB, nullable=False, comment="反应物SMILES列表")
    product_smiles = Column(JSONB, nullable=False, comment="产物SMILES列表")
    
    # 反应算符信息
    operator_name = Column(String(128), index=True, comment="反应算符名称")
    reaction_type = Column(String(128), index=True, comment="反应类型")
    
    # 能量学信息
    reaction_energy = Column(Float, comment="反应能 (kcal/mol)")
    activation_energy = Column(Float, comment="活化能 (kcal/mol)")
    
    # 反应条件
    driving_force = Column(String(128), comment="主要驱动力")
    
    # 反应平衡
    is_reversible = Column(Boolean, default=True, comment="是否可逆")
    equilibrium_constant = Column(Float, comment="平衡常数")
    
    # 额外信息
    reaction_metadata = Column(JSONB, default={}, comment="反应元数据")

    
    # 时间戳
    created_at = Column(DateTime, default=datetime.utcnow, comment="创建时间")
    
    # 关系
    job = relationship("ReactionNetworkJob", back_populates="reactions")
    
    def __repr__(self):
        return f"<Reaction(id={self.id}, operator='{self.operator_name}')>"
