import React, { useState, useEffect } from 'react';
import {
  Spin, message, Card, Row, Col, Button, Space, Typography, theme,
  Statistic, Progress, Tag, Segmented, Empty, Checkbox, Alert,
  Slider, Switch, Tooltip, Badge, Divider, Steps, Input, Select,
  Tabs, Timeline, Avatar, Rate, Collapse
} from 'antd';
import {
  ExperimentOutlined, ClusterOutlined, RocketOutlined,
  ThunderboltOutlined, StarOutlined, SettingOutlined,
  SearchOutlined, FilterOutlined, EyeOutlined, BulbOutlined,
  TrophyOutlined, FireOutlined, DashboardOutlined, LineChartOutlined
} from '@ant-design/icons';
import client from '../api/client';
import { useThemeStore } from '../stores/themeStore';
import { useAuthStore } from '../stores/authStore';
import MultiObjectiveDashboard from '../components/MultiObjectiveDashboard';
import './AIDiscovery.css';

const DEFAULT_PROPERTIES = ['bp', 'mp', 'fp', 'alpha', 'mu', 'gap', 'homo', 'lumo'];
const QM9_PROPERTIES = ['alpha', 'mu', 'gap', 'homo', 'lumo'];

// 目标预设模板
const GOAL_PRESETS = [
  {
    id: 'high_temp_wide_esw',
    name: '高温+更宽ESW',
    description: '适用于高温环境，需要更宽的电位窗口',
    objectives: [
      { property: 'gap', direction: 'up', weight: 0.4 },
      { property: 'fp', direction: 'up', weight: 0.3 },
      { property: 'bp', direction: 'up', weight: 0.3 }
    ],
    constraints: [
      { property: 'fp', operator: '>=', value: 'baseline' },
      { property: 'bp', operator: '>=', value: 'baseline' }
    ]
  },
  {
    id: 'safety_priority',
    name: '安全性优先',
    description: '优先考虑安全性，提高闪点和沸点',
    objectives: [
      { property: 'fp', direction: 'up', weight: 0.5 },
      { property: 'bp', direction: 'up', weight: 0.3 },
      { property: 'gap', direction: 'up', weight: 0.2 }
    ],
    constraints: [
      { property: 'fp', operator: '>=', value: 'baseline+10' },
      { property: 'bp', operator: '>=', value: 'baseline' }
    ]
  },
  {
    id: 'similar_replacement',
    name: '相似替代',
    description: '找到性质相近的替代分子',
    objectives: [
      { property: 'gap', direction: 'maintain', weight: 0.3 },
      { property: 'fp', direction: 'maintain', weight: 0.3 },
      { property: 'bp', direction: 'maintain', weight: 0.4 }
    ],
    constraints: [
      { property: 'gap', operator: '>=', value: 'baseline*0.9' },
      { property: 'fp', operator: '>=', value: 'baseline*0.9' }
    ]
  }
];

// 属性配置信息
const PROPERTY_CONFIG: Record<string, {
  name: string;
  unit: string;
  icon: string;
  description: string;
  category: string;
  defaultDirection: string;
  importance: string;
}> = {
  gap: {
    name: 'HOMO-LUMO Gap',
    unit: 'eV',
    icon: '⚡',
    description: '能隙，影响电子传输性能',
    category: 'electronic',
    defaultDirection: 'higher',
    importance: 'high'
  },
  homo: {
    name: 'HOMO Energy',
    unit: 'eV',
    icon: '🔋',
    description: '最高占据分子轨道能量',
    category: 'electronic',
    defaultDirection: 'higher',
    importance: 'medium'
  },
  lumo: {
    name: 'LUMO Energy',
    unit: 'eV',
    icon: '⚡',
    description: '最低未占据分子轨道能量',
    category: 'electronic',
    defaultDirection: 'lower',
    importance: 'medium'
  },
  bp: {
    name: 'Boiling Point',
    unit: '°C',
    icon: '🌡️',
    description: '沸点，影响挥发性',
    category: 'thermal',
    defaultDirection: 'higher',
    importance: 'high'
  },
  mp: {
    name: 'Melting Point',
    unit: '°C',
    icon: '❄️',
    description: '熔点，影响相态稳定性',
    category: 'thermal',
    defaultDirection: 'lower',
    importance: 'medium'
  },
  fp: {
    name: 'Flash Point',
    unit: '°C',
    icon: '🔥',
    description: '闪点，影响安全性',
    category: 'safety',
    defaultDirection: 'higher',
    importance: 'high'
  },
  alpha: {
    name: 'Polarizability',
    unit: 'Bohr³',
    icon: '🎯',
    description: '极化率，影响溶剂化性能',
    category: 'molecular',
    defaultDirection: 'higher',
    importance: 'low'
  },
  mu: {
    name: 'Dipole Moment',
    unit: 'Debye',
    icon: '🧲',
    description: '偶极矩，影响极性',
    category: 'molecular',
    defaultDirection: 'higher',
    importance: 'medium'
  }
};

const EXAMPLE_SMILES = `C1COC(=O)O1
COC(=O)OC
CCOC(=O)OCC
CCOC(=O)OC
CC1COC(=O)O1`;

// 常见电解液分子
const COMMON_ELECTROLYTE_MOLECULES = [
  { name: 'EC (碳酸乙烯酯)', smiles: 'C1COC(=O)O1', properties: {} },
  { name: 'DMC (碳酸二甲酯)', smiles: 'COC(=O)OC', properties: {} },
  { name: 'DEC (碳酸二乙酯)', smiles: 'CCOC(=O)OCC', properties: {} },
  { name: 'EMC (碳酸甲乙酯)', smiles: 'CCOC(=O)OC', properties: {} },
  { name: 'PC (碳酸丙烯酯)', smiles: 'CC1COC(=O)O1', properties: {} },
  { name: 'FEC (氟代碳酸乙烯酯)', smiles: 'FC1COC(=O)O1', properties: {} },
  { name: 'DME (乙二醇二甲醚)', smiles: 'COCCOC', properties: {} },
  { name: 'DOL (1,3-二氧戊环)', smiles: 'C1COCO1', properties: {} },
  { name: 'THF (四氢呋喃)', smiles: 'C1CCOC1', properties: {} },
  { name: 'ACN (乙腈)', smiles: 'CC#N', properties: {} },
  { name: 'GBL (γ-丁内酯)', smiles: 'C1CCC(=O)O1', properties: {} },
  { name: 'Sulfolane (环丁砜)', smiles: 'C1CCS(=O)(=O)C1', properties: {} },
];

// 默认示例分子数据（使用数据库中实际存在的分子）
const DEFAULT_ELECTROLYTE_EXAMPLES = [
  {
    name: '丙炔醇',
    smiles: 'C#CCO',
    properties: {
      bp: 113.61,
      mp: -47.78,
      fp: 32.78,
      gap: 8.87,
      alpha: 18.78,
      mu: 1.64,
      homo: -7.36,
      lumo: 1.00
    },
    similarity: 1.0,
    molecule_info: {
      name: 'Propargyl alcohol',
      molecular_formula: 'C3H4O'
    }
  },
  {
    name: '丙炔酸',
    smiles: 'C#CC(=O)O',
    properties: {
      bp: 144.0,
      mp: 18.0,
      fp: 58.0,
      gap: 7.40,
      alpha: 19.29,
      mu: 1.56,
      homo: -7.65,
      lumo: -0.86
    },
    similarity: 0.85,
    molecule_info: {
      name: 'Propiolic acid',
      molecular_formula: 'C3H2O2'
    }
  },
  {
    name: '溴丙炔',
    smiles: 'C#CCBr',
    properties: {
      bp: 88.89,
      mp: -61.06,
      fp: 10.0,
      gap: 9.60,
      alpha: 17.25,
      mu: 1.79,
      homo: -7.69,
      lumo: 0.87
    },
    similarity: 0.78,
    molecule_info: {
      name: 'Propargyl bromide',
      molecular_formula: 'C3H3Br'
    }
  },
  {
    name: '氰基溴',
    smiles: 'C(#N)Br',
    properties: {
      bp: 61.11,
      mp: -73.55,
      fp: 22.01,
      gap: 10.54,
      alpha: 10.60,
      mu: 2.25,
      homo: -10.17,
      lumo: 0.78
    },
    similarity: 0.82,
    molecule_info: {
      name: 'Cyanogen bromide',
      molecular_formula: 'CBrN'
    }
  },
  {
    name: '丙二腈',
    smiles: 'C(#N)C#N',
    properties: {
      bp: -21.17,
      mp: -27.89,
      fp: 22.01,
      gap: 7.37,
      alpha: 15.97,
      mu: 2.65,
      homo: -9.98,
      lumo: -0.74
    },
    similarity: 0.88,
    molecule_info: {
      name: 'Malononitrile',
      molecular_formula: 'C3N2'
    }
  }
];

interface Molecule {
  smiles: string;
  name?: string;
  properties?: Record<string, any>;
  image?: string;
  index?: number;
  distance?: number;
  is_center?: boolean;
  cas_number?: string;
  molecular_formula?: string;
  molecular_weight?: number;
  is_real_data?: Record<string, boolean>;
  similarity?: number;
  cluster_id?: number;
}

interface ClusterData {
  center: Molecule;
  neighbors: Molecule[];
}

const { Title, Paragraph, Text } = Typography;

const AIDiscovery: React.FC = () => {
  const { mode } = useThemeStore();
  const { token } = theme.useToken();
  const isDark = mode === 'dark';
  const { user } = useAuthStore();

  // 检查用户是否有ai-discovery模块权限
  const hasAIDiscoveryAccess = () => {
    if (!user) return false;
    // 管理员有所有权限
    if (user.role === 'ADMIN') return true;
    // 如果没有allowed_modules，默认有权限（向后兼容）
    if (!user.allowed_modules || user.allowed_modules.length === 0) return true;
    // 检查是否有ai-discovery权限
    return user.allowed_modules.includes('ai-discovery');
  };

  const canUseOptimization = hasAIDiscoveryAccess();

  // 页面状态 - 重构为三模式
  const [currentMode, setCurrentMode] = useState<'evaluate' | 'alternatives' | 'optimize'>('evaluate');
  const [loading, setLoading] = useState(false);
  const [activeStep, setActiveStep] = useState(0);

  // 任务摘要状态
  const [taskSummary, setTaskSummary] = useState({
    baselineMolecule: null as Molecule | null,
    goals: [] as string[],
    constraints: [] as string[]
  });

  // Alternatives模式状态
  const [selectedPreset, setSelectedPreset] = useState<string>('');
  const [alternativesResults, setAlternativesResults] = useState<any[]>([]);
  const [alternativesLoading, setAlternativesLoading] = useState(false);

  // 预测相关状态
  const [smilesInput, setSmilesInput] = useState('');
  const [smilesCount, setSmilesCount] = useState(0);
  const [predictResults, setPredictResults] = useState<Molecule[]>([]);
  const [totalMolecules, setTotalMolecules] = useState(0);
  const [totalSuccess, setTotalSuccess] = useState(0);
  const [lastTime, setLastTime] = useState(0);

  // 聚类相关状态
  const [clusterData, setClusterData] = useState<ClusterData | null>(null);
  const [clusterHistory, setClusterHistory] = useState<Molecule[]>([]);
  const [centerMolecules, setCenterMolecules] = useState<Molecule[]>([]);
  const [loadingTemplates, setLoadingTemplates] = useState(false);

  // 搜索配置状态
  const [searchMode, setSearchMode] = useState<'similarity' | 'optimization'>('similarity');
  const [optimizationTargets, setOptimizationTargets] = useState<Record<string, any>>({});
  const [selectedProperties, setSelectedProperties] = useState<string[]>([]);
  const [propertyWeights, setPropertyWeights] = useState<Record<string, number>>({});
  const [searchConfig, setSearchConfig] = useState({
    similarityWeight: 0.7,
    optimizationWeight: 0.3,
    maxResults: 20,
    scope: 'global'
  });

  // UI状态
  const [hoveredMolecule, setHoveredMolecule] = useState<string | null>(null);
  const [selectedMolecule, setSelectedMolecule] = useState<string | null>(null);
  const [showAdvanced, setShowAdvanced] = useState(false);

  // 多目标优化状态
  const [objectiveWeights, setObjectiveWeights] = useState<Record<string, number>>({
    bp: 0.0, mp: 0.0, fp: 0.0, gap: 1.0, alpha: 0.0, mu: 0.0, homo: 0.0, lumo: 0.0
  });
  const [optimizationMethod, setOptimizationMethod] = useState<'pareto' | 'weighted'>('pareto');
  const [multiObjectiveResults, setMultiObjectiveResults] = useState<any[]>([]);
  const [multiObjectiveLoading, setMultiObjectiveLoading] = useState(false);
  const [showParetoChart, setShowParetoChart] = useState(false);

  useEffect(() => {
    const count = smilesInput.split(/\n+/).filter(s => s.trim()).length;
    setSmilesCount(count);
  }, [smilesInput]);

  // 初始化默认示例数据
  useEffect(() => {
    // 设置默认的预测结果示例
    setPredictResults(DEFAULT_ELECTROLYTE_EXAMPLES);
    setTotalMolecules(DEFAULT_ELECTROLYTE_EXAMPLES.length);
    setTotalSuccess(DEFAULT_ELECTROLYTE_EXAMPLES.length);
    setLastTime(0.85); // 模拟预测时间
  }, []);

  // 加载中心分子列表
  const loadCenterMolecules = async () => {
    try {
      setLoadingTemplates(true);
      const response = await client.get('/ai-discovery/templates', {
        params: { limit: 40 }
      });

      const templates = response.data.templates.map((t: any) => ({
        smiles: t.smiles,
        name: t.name || `分子 ${t.smiles.substring(0, 20)}...`,
        properties: t.properties || {},
        cluster_id: t.cluster_id
      }));

      // 合并常见电解液分子和数据库分子，去重
      const commonMolecules = COMMON_ELECTROLYTE_MOLECULES.map(mol => ({
        ...mol,
        isCommon: true
      }));

      // 过滤掉与常见分子重复的数据库分子
      const uniqueTemplates = templates.filter((t: any) =>
        !commonMolecules.some(c => c.smiles === t.smiles)
      );

      // 常见分子排在前面
      const allMolecules = [...commonMolecules, ...uniqueTemplates];

      setCenterMolecules(allMolecules);
    } catch (error) {
      message.error('加载中心分子列表失败');
      console.error(error);
    } finally {
      setLoadingTemplates(false);
    }
  };

  // 当切换到alternatives或optimize模式时加载模板
  useEffect(() => {
    if ((currentMode === 'alternatives' || currentMode === 'optimize') && centerMolecules.length === 0) {
      loadCenterMolecules();
    }
  }, [currentMode]);

  // 辅助函数：更新任务摘要
  const updateTaskSummary = (baselineMolecule: Molecule | null, goals: string[] = [], constraints: string[] = []) => {
    setTaskSummary({
      baselineMolecule,
      goals,
      constraints
    });
  };

  // 辅助函数：从评估模式跳转到其他模式
  const jumpToAlternatives = (molecule: Molecule) => {
    updateTaskSummary(molecule, ['找到相似替代'], []);
    setCurrentMode('alternatives');
  };

  const jumpToOptimize = (molecule: Molecule) => {
    updateTaskSummary(molecule, ['多目标优化'], []);
    setCurrentMode('optimize');
  };

  // 多目标优化处理函数
  const handleMultiObjectiveOptimization = async () => {
    // 检查用户是否已登录
    if (!user) {
      message.error('请先登录后使用多目标优化功能');
      return;
    }

    if (!taskSummary.baselineMolecule) {
      message.error('请先选择基准分子');
      return;
    }

    // 检查是否至少选择了一个目标
    const activeObjectives = Object.entries(objectiveWeights).filter(([_, weight]) => Math.abs(weight) > 0);
    if (activeObjectives.length === 0) {
      message.error('请至少设置一个优化目标的权重大于0');
      return;
    }

    setMultiObjectiveLoading(true);
    try {
      // 构建优化目标
      const optimizationTargets: Record<string, any> = {};
      Object.entries(objectiveWeights).forEach(([prop, weight]) => {
        if (Math.abs(weight) > 0) {
          optimizationTargets[prop] = {
            direction: weight > 0 ? 'up' : 'down',
            weight: Math.abs(weight),
            target_value: taskSummary.baselineMolecule?.properties?.[prop] || 0
          };
        }
      });

      // 调用多目标优化API
      const response = await client.post('/ai-discovery/multi-objective-optimize', {
        smiles: taskSummary.baselineMolecule.smiles,
        n_results: 50, // 减少候选分子数量以提高性能
        scope: 'global',
        optimization_targets: optimizationTargets,
        similarity_weight: 0.3,
        optimization_method: optimizationMethod
      }, {
        timeout: 120000 // 增加超时时间到2分钟
      });

      console.log('API响应:', response.data); // 调试信息

      // 转换API结果为前端格式
      const candidates = response.data.candidates || [];
      const formattedResults = candidates.map((candidate: any, idx: number) => ({
        smiles: candidate.smiles,
        name: `候选分子 #${idx + 1}`,
        properties: candidate.properties,
        similarity: candidate.similarity,
        pareto_rank: 1, // 可以根据实际算法计算
        dominance_count: 0,
        objective_values: candidate.objective_values,
        molecule_info: candidate.molecule_info || {} // 添加分子信息
      }));

      if (formattedResults.length > 0) {
        setMultiObjectiveResults(formattedResults);
        setShowParetoChart(true);
        message.success(`成功找到 ${formattedResults.length} 个候选分子`);
      } else {
        // 如果没有找到候选分子，提供提示并使用演示数据
        message.warning('数据库中未找到相似分子，显示演示数据');
        const demoResults = generateDemoMultiObjectiveResults(taskSummary.baselineMolecule);
        setMultiObjectiveResults(demoResults);
        setShowParetoChart(true);
      }
    } catch (error) {
      console.error('多目标优化失败:', error);
      // 如果API调用失败，使用演示数据
      const demoResults = generateDemoMultiObjectiveResults(taskSummary.baselineMolecule);
      setMultiObjectiveResults(demoResults);
      setShowParetoChart(true);
      message.warning(`API调用失败，使用演示数据 (${demoResults.length} 个分子)`);
    } finally {
      setMultiObjectiveLoading(false);
    }
  };

  // 生成演示多目标优化结果
  const generateDemoMultiObjectiveResults = (baseline: Molecule) => {
    const demoSmiles = [
      // 醇类
      'CCO', 'CC(C)O', 'CCCCO', 'CC(C)(C)O', 'CCCO', 'CCCCCO', 'CC(C)CO',
      // 酯类
      'CC(=O)OCC', 'COC(=O)OC', 'CCOC(=O)OCC', 'CC(=O)OCCC', 'CCC(=O)OCC',
      // 醚类
      'CCOCC', 'CCCOCC', 'CC(C)OCC', 'CCCCOCC',
      // 酸类
      'CC(=O)O', 'CCC(=O)O', 'CCCC(=O)O', 'CC(C)C(=O)O',
      // 芳香族
      'c1ccccc1', 'Cc1ccccc1', 'CCc1ccccc1', 'c1ccc(O)cc1', 'c1ccc(C)cc1',
      // 环状化合物
      'C1COC(=O)O1', 'C1CCOC1', 'C1CCOCC1', 'C1COCCO1',
      // 卤代化合物
      'CCCl', 'CCBr', 'CF', 'CCF',
      // 胺类
      'CCN', 'CCCN', 'CC(C)N', 'CCNCC',
      // 腈类
      'CC#N', 'CCC#N', 'CCCC#N',
      // 其他
      'CCS', 'CCCS', 'CC(=O)N', 'CCCN(C)C'
    ];

    // 常见分子的名称映射
    const moleculeNames: Record<string, { name: string; formula: string }> = {
      'CCO': { name: '乙醇', formula: 'C2H6O' },
      'CC(C)O': { name: '异丙醇', formula: 'C3H8O' },
      'CCCCO': { name: '1-丁醇', formula: 'C4H10O' },
      'CC(C)(C)O': { name: '叔丁醇', formula: 'C4H10O' },
      'CCCO': { name: '1-丙醇', formula: 'C3H8O' },
      'CC(=O)OCC': { name: '乙酸乙酯', formula: 'C4H8O2' },
      'COC(=O)OC': { name: '碳酸二甲酯', formula: 'C3H6O3' },
      'CCOC(=O)OCC': { name: '碳酸二乙酯', formula: 'C5H10O3' },
      'CC(=O)O': { name: '乙酸', formula: 'C2H4O2' },
      'c1ccccc1': { name: '苯', formula: 'C6H6' },
      'Cc1ccccc1': { name: '甲苯', formula: 'C7H8' },
      'C1COC(=O)O1': { name: '碳酸乙烯酯', formula: 'C3H4O3' },
      'C1CCOC1': { name: '四氢呋喃', formula: 'C4H8O' },
      'CCOCC': { name: '乙醚', formula: 'C4H10O' },
      'CCCl': { name: '氯丙烷', formula: 'C3H7Cl' },
      'CCN': { name: '乙胺', formula: 'C2H7N' },
      'CC#N': { name: '乙腈', formula: 'C2H3N' }
    };

    return demoSmiles.slice(0, 50).map((smiles, idx) => {
      const properties: Record<string, number> = {};
      DEFAULT_PROPERTIES.forEach(prop => {
        const baseValue = baseline.properties?.[prop] || 0;
        const variation = (Math.random() - 0.5) * 0.6 * Math.abs(baseValue);
        properties[prop] = Math.max(0, baseValue + variation);
      });

      // 获取分子信息
      const moleculeInfo = moleculeNames[smiles] || {
        name: `分子 #${idx + 1}`,
        formula: 'Unknown'
      };

      return {
        smiles,
        name: `候选分子 #${idx + 1}`,
        properties,
        similarity: 0.6 + Math.random() * 0.3,
        pareto_rank: Math.floor(Math.random() * 4) + 1,
        dominance_count: Math.floor(Math.random() * 8),
        objective_values: {}, // 可以根据需要计算
        molecule_info: {
          name: moleculeInfo.name,
          molecular_formula: moleculeInfo.formula,
          smiles: smiles
        }
      };
    });
  };

  // 生成演示替代分子数据
  const generateDemoAlternatives = (baseline: Molecule) => {
    const demoSmiles = [
      'CCO', 'CC(C)O', 'CCCCO', 'CC(C)(C)O', 'CCCO',
      'CC(=O)OCC', 'CC(=O)O', 'CCC(=O)O', 'CCCC(=O)O',
      'c1ccccc1', 'Cc1ccccc1', 'CCc1ccccc1', 'c1ccc(O)cc1'
    ];

    return demoSmiles.slice(0, 8).map((smiles, idx) => {
      // 基于基准分子生成变化的属性值
      const properties: Record<string, number> = {};
      DEFAULT_PROPERTIES.forEach(prop => {
        const baseValue = baseline.properties?.[prop] || 0;
        // 添加一些随机变化
        const variation = (Math.random() - 0.5) * 0.3 * Math.abs(baseValue);
        properties[prop] = baseValue + variation;
      });

      return {
        smiles,
        name: `演示分子 #${idx + 1}`,
        properties,
        similarity: 0.8 - idx * 0.05, // 递减的相似度
        distance: idx * 0.1 + 0.1
      };
    });
  };

  // 处理预设选择
  const handlePresetSelect = (presetId: string) => {
    setSelectedPreset(presetId);
    const preset = GOAL_PRESETS.find(p => p.id === presetId);
    if (preset && taskSummary.baselineMolecule) {
      const goals = preset.objectives.map(obj =>
        `${PROPERTY_CONFIG[obj.property]?.name || obj.property} ${obj.direction === 'up' ? '↑' : obj.direction === 'down' ? '↓' : '='}`
      );
      const constraints = preset.constraints.map(cons =>
        `${PROPERTY_CONFIG[cons.property]?.name || cons.property} ${cons.operator} ${cons.value}`
      );
      updateTaskSummary(taskSummary.baselineMolecule, goals, constraints);
    }
  };

  // 搜索替代分子
  const handleSearchAlternatives = async () => {
    if (!taskSummary.baselineMolecule) {
      message.error('请先选择基准分子');
      return;
    }

    setAlternativesLoading(true);
    try {
      // 使用测试端点进行相似分子搜索（无需认证）
      const response = await client.get('/ai-discovery/test-similar-molecules', {
        params: {
          smiles: taskSummary.baselineMolecule.smiles,
          n_similar: 50,
          scope: 'global'
        }
      });

      let candidates = response.data.similar_molecules || [];

      // 如果还是没有结果，使用模拟数据进行演示
      if (candidates.length === 0) {
        message.warning('未找到相似分子，使用演示数据');
        candidates = generateDemoAlternatives(taskSummary.baselineMolecule);
      }

      // 计算Δ指标和得分
      const processedCandidates = candidates.map((candidate: any) => {
        const deltas: Record<string, number> = {};
        let objectiveScore = 0;
        let constraintsMet = true;

        // 计算Δ值
        DEFAULT_PROPERTIES.forEach(prop => {
          const candidateVal = candidate.properties?.[prop];
          const baselineVal = taskSummary.baselineMolecule?.properties?.[prop];
          if (candidateVal !== undefined && baselineVal !== undefined) {
            deltas[prop] = candidateVal - baselineVal;
          }
        });

        // 计算目标得分（如果有选择预设）
        if (selectedPreset) {
          const preset = GOAL_PRESETS.find(p => p.id === selectedPreset);
          if (preset) {
            preset.objectives.forEach(obj => {
              const delta = deltas[obj.property];
              if (delta !== undefined) {
                if (obj.direction === 'up' && delta > 0) {
                  objectiveScore += delta * obj.weight;
                } else if (obj.direction === 'down' && delta < 0) {
                  objectiveScore += Math.abs(delta) * obj.weight;
                }
              }
            });

            // 检查约束
            preset.constraints.forEach(cons => {
              const candidateVal = candidate.properties?.[cons.property];
              const baselineVal = taskSummary.baselineMolecule?.properties?.[cons.property];
              if (candidateVal !== undefined && baselineVal !== undefined) {
                let threshold = baselineVal;
                if (cons.value === 'baseline+10') threshold = baselineVal + 10;
                else if (cons.value === 'baseline*0.9') threshold = baselineVal * 0.9;

                if (cons.operator === '>=' && candidateVal < threshold) {
                  constraintsMet = false;
                }
              }
            });
          }
        }

        return {
          ...candidate,
          deltas,
          objectiveScore,
          constraintsMet
        };
      });

      // 排序：达标优先，然后按目标得分
      processedCandidates.sort((a: any, b: any) => {
        if (a.constraintsMet !== b.constraintsMet) {
          return a.constraintsMet ? -1 : 1;
        }
        return b.objectiveScore - a.objectiveScore;
      });

      setAlternativesResults(processedCandidates);
      message.success(`找到 ${processedCandidates.length} 个候选分子`);
    } catch (error) {
      console.error('搜索替代分子失败:', error);
      message.error('搜索失败，请重试');
    } finally {
      setAlternativesLoading(false);
    }
  };

  const handlePredict = async () => {
    const smilesList = smilesInput.split(/\n+/).map(s => s.trim()).filter(Boolean);
    if (!smilesList.length) {
      message.error('请至少输入一个 SMILES');
      return;
    }

    try {
      setLoading(true);
      const startTime = performance.now();

      // 使用批量预测 API
      const response = await client.post('/ai-discovery/predict-batch', {
        smiles_list: smilesList,
        properties: DEFAULT_PROPERTIES,
        include_qm9: true
      });

      const results: Molecule[] = response.data.results.map((result: any, i: number) => {
        const moleculeInfo = result.molecule_info || {};
        const displayName = moleculeInfo.name || moleculeInfo.molecular_formula || `分子 #${i + 1}`;

        return {
          smiles: result.smiles,
          name: displayName,
          properties: result.predicted_properties || {},
          index: i,
          cas_number: moleculeInfo.cas_number,
          molecular_formula: moleculeInfo.molecular_formula,
          molecular_weight: moleculeInfo.molecular_weight,
          is_real_data: result.is_real_data || {},
          image: result.image  // 添加分子结构图像
        };
      });

      const elapsed = (performance.now() - startTime) / 1000;
      const successCount = response.data.success_count;

      setTotalMolecules(prev => prev + response.data.total_count);
      setTotalSuccess(prev => prev + successCount);
      setLastTime(elapsed);
      setPredictResults(results);

      if (response.data.failed_smiles.length > 0) {
        message.warning(`预测完成！成功: ${successCount}/${response.data.total_count}，失败: ${response.data.failed_smiles.length}`);
      } else {
        message.success(`预测完成！成功: ${successCount}/${response.data.total_count}`);
      }
    } catch (error) {
      message.error('预测失败');
      console.error(error);
    } finally {
      setLoading(false);
    }
  };

  const handleClusterExplore = async (molecule: Molecule, scope: string = 'global', addToHistory: boolean = true) => {
    try {
      setLoading(true);

      // 将分子加入历史记录
      if (addToHistory) {
        setClusterHistory([...clusterHistory, molecule]);
      }

      let response;

      // 根据搜索模式选择不同的端点
      if (searchMode === 'optimization') {
        // 属性优化搜索
        const optimizationPayload: any = {
          smiles: molecule.smiles,
          n_results: 50,
          scope: scope,
          similarity_weight: 0.5
        };

        // 构建优化目标
        const targets: Record<string, any> = {};
        Object.entries(optimizationTargets).forEach(([prop, config]: [string, any]) => {
          if (config && config.weight && config.weight > 0) {
            targets[prop] = {
              direction: config.direction || 'higher',
              weight: config.weight || 1
            };
          }
        });

        if (Object.keys(targets).length > 0) {
          optimizationPayload.optimization_targets = targets;
        }

        response = await client.post('/ai-discovery/optimize-molecules', optimizationPayload);

        // 转换响应格式以匹配相似分子的格式
        const optimizedMols = response.data.optimized_molecules || [];
        response.data.similar_molecules = optimizedMols.map((mol: any) => ({
          smiles: mol.smiles,
          properties: mol.properties,
          similarity: mol.combined_score,
          molecule_info: { name: mol.smiles.substring(0, 20) + '...' },
          is_real_data: {}
        }));
      } else {
        // 相似性搜索（原有逻辑）
        response = await client.get('/ai-discovery/similar-molecules', {
          params: {
            smiles: molecule.smiles,
            n_similar: 50,  // 增加到50个相似分子，提供更多选择
            scope: scope  // 支持 'cluster' 或 'global' 范围
          }
        });
      }

      // 转换相似分子数据格式
      const similarMolecules = (response.data.similar_molecules || []).map((sim: any, idx: number) => {
        // 获取分子信息
        const moleculeInfo = sim.molecule_info || {};
        const displayName = moleculeInfo.name || moleculeInfo.molecular_formula || `相似分子 #${idx + 1}`;

        return {
          smiles: sim.smiles,
          name: displayName,
          properties: sim.properties || {},
          similarity: sim.similarity,
          distance: 1 - (sim.similarity || 0), // 转换相似度为距离
          index: idx,
          cas_number: moleculeInfo.cas_number,
          molecular_formula: moleculeInfo.molecular_formula,
          molecular_weight: moleculeInfo.molecular_weight,
          is_real_data: sim.is_real_data || {}
        };
      });

      setClusterData({
        center: molecule,
        neighbors: similarMolecules
      });
      setClusterHistory(prev => [...prev, molecule]);
    } catch (error) {
      message.error('聚类探索失败');
      console.error(error);
    } finally {
      setLoading(false);
    }
  };

  const handleHistoryBack = () => {
    if (clusterHistory.length > 1) {
      const newHistory = clusterHistory.slice(0, -1);
      const prevMolecule = newHistory[newHistory.length - 1];
      setClusterHistory(newHistory);
      // 不添加到历史记录，因为我们只是返回
      handleClusterExplore(prevMolecule, 'global', false);
    }
  };

  const handleHistoryHome = () => {
    if (clusterHistory.length > 0) {
      const homeMolecule = clusterHistory[0];
      setClusterHistory([homeMolecule]);
      // 不添加到历史记录，因为我们只是返回到起点
      handleClusterExplore(homeMolecule, 'global', false);
    }
  };

  const handleBackToSelection = () => {
    // 清空聚类数据，返回到选择中心分子的界面
    setClusterData(null);
    setClusterHistory([]);
  };

  const cardStyle: React.CSSProperties = {
    borderRadius: 12,
    boxShadow: isDark ? '0 4px 12px rgba(0, 0, 0, 0.3)' : '0 4px 12px rgba(0, 0, 0, 0.06)',
    border: `1px solid ${token.colorBorder}`,
    background: token.colorBgContainer,
  };

  const statCardStyle: React.CSSProperties = {
    borderRadius: 12,
    boxShadow: isDark ? '0 4px 12px rgba(0, 0, 0, 0.3)' : '0 4px 12px rgba(0, 0, 0, 0.06)',
    border: `1px solid ${token.colorBorder}`,
    background: token.colorBgContainer,
    padding: '20px',
    height: '100%',
  };

  return (
    <div style={{ padding: '24px', background: token.colorBgLayout, minHeight: '100vh' }}>
      <Spin spinning={loading}>
        {/* 页面标题 */}
        <div style={{ marginBottom: '32px' }}>
          <Title level={2} style={{ margin: '0 0 8px 0' }}>🔬 分子属性预测平台</Title>
          <Paragraph style={{ margin: 0, color: token.colorTextSecondary }}>
            基于深度学习的分子性质预测系统，支持 alpha、mu、gap、homo、lumo、BP、FP、MP 等属性预测
          </Paragraph>
        </div>

        {/* 统计卡片 */}
        <Row gutter={[16, 16]} style={{ marginBottom: '32px' }}>
          <Col xs={24} sm={12} lg={6}>
            <Card style={statCardStyle} bordered={false}>
              <Statistic title="支持属性数" value={8} prefix="⚡" />
            </Card>
          </Col>
          <Col xs={24} sm={12} lg={6}>
            <Card style={statCardStyle} bordered={false}>
              <Statistic title="已预测分子" value={totalMolecules} prefix="🧬" />
            </Card>
          </Col>
          <Col xs={24} sm={12} lg={6}>
            <Card style={statCardStyle} bordered={false}>
              <Statistic title="成功预测" value={totalSuccess} prefix="✅" />
            </Card>
          </Col>
          <Col xs={24} sm={12} lg={6}>
            <Card style={statCardStyle} bordered={false}>
              <Statistic title="上次耗时" value={lastTime} suffix="s" precision={2} prefix="⏱️" />
            </Card>
          </Col>
        </Row>

        {/* 三模式切换 - 按钮组形式 */}
        <div style={{ marginBottom: '24px', textAlign: 'center' }}>
          <Space size="large">
            <Button
              type={currentMode === 'evaluate' ? 'primary' : 'default'}
              size="large"
              icon={<ExperimentOutlined />}
              onClick={() => setCurrentMode('evaluate')}
              style={{
                height: '48px',
                minWidth: '160px',
                borderRadius: '8px',
                fontWeight: currentMode === 'evaluate' ? 'bold' : 'normal'
              }}
            >
              评估单分子
            </Button>
            <Button
              type={currentMode === 'alternatives' ? 'primary' : 'default'}
              size="large"
              icon={<SearchOutlined />}
              onClick={() => setCurrentMode('alternatives')}
              style={{
                height: '48px',
                minWidth: '160px',
                borderRadius: '8px',
                fontWeight: currentMode === 'alternatives' ? 'bold' : 'normal'
              }}
            >
              找更好的替代
            </Button>
            <Button
              type={currentMode === 'optimize' ? 'primary' : 'default'}
              size="large"
              icon={<RocketOutlined />}
              onClick={() => setCurrentMode('optimize')}
              style={{
                height: '48px',
                minWidth: '160px',
                borderRadius: '8px',
                fontWeight: currentMode === 'optimize' ? 'bold' : 'normal'
              }}
            >
              按目标改进
            </Button>
          </Space>
        </div>

        {/* 任务摘要条 */}
        {taskSummary.baselineMolecule && (
          <Card
            size="small"
            style={{
              marginBottom: '24px',
              background: token.colorPrimaryBg,
              border: `1px solid ${token.colorPrimary}`,
            }}
          >
            <div style={{ display: 'flex', alignItems: 'center', gap: '16px', flexWrap: 'wrap' }}>
              <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
                <Text strong>基准分子:</Text>
                <Tag color="blue">
                  {taskSummary.baselineMolecule.name || taskSummary.baselineMolecule.smiles.substring(0, 20) + '...'}
                </Tag>
              </div>
              {taskSummary.goals.length > 0 && (
                <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
                  <Text strong>目标:</Text>
                  {taskSummary.goals.map((goal, idx) => (
                    <Tag key={idx} color="green">{goal}</Tag>
                  ))}
                </div>
              )}
              {taskSummary.constraints.length > 0 && (
                <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
                  <Text strong>约束:</Text>
                  {taskSummary.constraints.map((constraint, idx) => (
                    <Tag key={idx} color="orange">{constraint}</Tag>
                  ))}
                </div>
              )}
            </div>
          </Card>
        )}

        {/* 主要内容区域 */}
        {currentMode === 'evaluate' && (
          <Row gutter={24}>
            {/* 左侧：输入区域 */}
            <Col xs={24} lg={10}>
              <Card title={
                <span>
                  <ExperimentOutlined style={{ marginRight: '8px' }} />
                  分子输入
                </span>
              } style={cardStyle}>
                <div style={{ marginBottom: '20px' }}>
                  <Text strong style={{ display: 'block', marginBottom: '8px' }}>📝 SMILES 输入</Text>
                  <textarea
                    value={smilesInput}
                    onChange={(e) => setSmilesInput(e.target.value)}
                    placeholder="请输入 SMILES 字符串，每行一个&#10;例如：&#10;CCO&#10;CC(=O)O&#10;c1ccccc1"
                    style={{
                      width: '100%',
                      minHeight: '200px',
                      borderRadius: '8px',
                      border: `1px solid ${token.colorBorder}`,
                      padding: '12px',
                      background: isDark ? '#1a1a1a' : '#fafafa',
                      color: token.colorText,
                      fontFamily: 'monospace',
                      fontSize: '14px',
                      resize: 'vertical',
                    }}
                  />
                  <div style={{ display: 'flex', justifyContent: 'space-between', marginTop: '8px', fontSize: '13px', color: token.colorTextSecondary }}>
                    <span>支持批量预测，最多 100 个分子</span>
                    <span>已输入: {smilesCount} 个</span>
                  </div>
                </div>

                <Space style={{ marginBottom: '24px', width: '100%' }} direction="vertical">
                  <Button
                    type="primary"
                    size="large"
                    onClick={handlePredict}
                    loading={loading}
                    block
                    disabled={smilesCount === 0}
                    style={{
                      background: `linear-gradient(135deg, ${token.colorPrimary} 0%, ${token.colorPrimaryActive} 100%)`,
                      border: 'none',
                      borderRadius: '8px',
                      height: '48px',
                      fontSize: '16px',
                      fontWeight: 'bold',
                      boxShadow: `0 4px 12px ${token.colorPrimary}20`
                    }}
                  >
                    🚀 开始预测
                  </Button>
                  <Row gutter={8} style={{ marginTop: '12px' }}>
                    <Col span={12}>
                      <Button
                        onClick={() => setSmilesInput('')}
                        block
                        style={{
                          borderRadius: '6px',
                          height: '36px',
                          border: `1px solid ${token.colorBorder}`,
                          background: token.colorBgElevated,
                          color: token.colorTextSecondary
                        }}
                      >
                        🗑️ 清空
                      </Button>
                    </Col>
                    <Col span={12}>
                      <Button
                        onClick={() => setSmilesInput(EXAMPLE_SMILES)}
                        block
                        style={{
                          borderRadius: '6px',
                          height: '36px',
                          border: `1px solid ${token.colorPrimary}`,
                          background: `${token.colorPrimary}10`,
                          color: token.colorPrimary
                        }}
                      >
                        📋 电解液示例
                      </Button>
                    </Col>
                  </Row>
                </Space>

                {/* 预测进度 */}
                {loading && (
                  <div style={{ textAlign: 'center', padding: '20px' }}>
                    <Spin size="large" />
                    <Text style={{ display: 'block', marginTop: '12px', color: token.colorTextSecondary }}>
                      正在预测分子属性...
                    </Text>
                  </div>
                )}
              </Card>
            </Col>

            {/* 右侧：结果区域 */}
            <Col xs={24} lg={14}>
              <Card title={
                <span>
                  <DashboardOutlined style={{ marginRight: '8px' }} />
                  预测结果
                </span>
              } style={cardStyle}>

                {predictResults.length === 0 ? (
                  <Empty
                    description="请在左侧输入SMILES并点击预测"
                    style={{ marginTop: '40px', marginBottom: '40px' }}
                    image={Empty.PRESENTED_IMAGE_SIMPLE}
                  />
                ) : (
                  <div>
                    {/* 紧凑的表格形式显示 */}
                    <div style={{ marginBottom: '16px' }}>
                      <Text style={{ fontSize: '12px', color: token.colorTextSecondary }}>
                        共预测 {predictResults.length} 个分子
                      </Text>
                    </div>

                    <div style={{ maxHeight: '600px', overflowY: 'auto' }}>
                      {predictResults.map((mol, idx) => (
                        <div
                          key={idx}
                          style={{
                            marginBottom: '12px',
                            border: `1px solid ${token.colorBorder}`,
                            borderRadius: '6px',
                            overflow: 'hidden',
                            background: isDark ? '#0a0a0a' : '#fafafa',
                          }}
                        >
                          <div style={{ display: 'flex', gap: '0', alignItems: 'stretch' }}>
                            {/* 左侧：分子结构图 + 名称 + SMILES + 分子式 */}
                            <div style={{
                              width: '140px',
                              flexShrink: 0,
                              background: isDark ? '#1a1a1a' : '#f5f5f5',
                              borderRight: `1px solid ${token.colorBorder}`,
                              display: 'flex',
                              flexDirection: 'column',
                              alignItems: 'center',
                              justifyContent: 'flex-start',
                              padding: '8px',
                            }}>
                              {/* 分子结构图 */}
                              <div style={{
                                width: '80px',
                                height: '80px',
                                background: isDark ? '#0a0a0a' : '#ffffff',
                                borderRadius: '4px',
                                display: 'flex',
                                alignItems: 'center',
                                justifyContent: 'center',
                                border: `1px solid ${token.colorBorder}`,
                                overflow: 'hidden',
                                marginBottom: '6px',
                              }}>
                                {mol.image ? (
                                  <img src={`data:image/png;base64,${mol.image}`} alt="分子结构" style={{ maxWidth: '100%', maxHeight: '100%', objectFit: 'contain' }} />
                                ) : (
                                  <div style={{ fontSize: '28px' }}>🧪</div>
                                )}
                              </div>

                              {/* 分子名称 */}
                              <Text style={{ fontSize: '10px', textAlign: 'center', lineHeight: '1.2', color: token.colorTextSecondary, marginBottom: '4px', fontWeight: 'bold' }}>
                                {mol.name ? mol.name.substring(0, 14) : `分子#${idx + 1}`}
                              </Text>

                              {/* SMILES */}
                              <Text code style={{ fontSize: '8px', textAlign: 'center', lineHeight: '1.2', color: token.colorTextSecondary, marginBottom: '4px', wordBreak: 'break-all' }}>
                                {mol.smiles.length > 20 ? `${mol.smiles.substring(0, 20)}...` : mol.smiles}
                              </Text>

                              {/* 分子式 */}
                              {mol.molecular_formula && (
                                <Tag color="blue" style={{ fontSize: '9px', marginBottom: '0' }}>
                                  {mol.molecular_formula}
                                </Tag>
                              )}
                            </div>

                            {/* 中间：属性网格 */}
                            <div style={{ flex: 1, padding: '8px 12px', display: 'flex', flexDirection: 'column', justifyContent: 'center' }}>
                              {/* 属性网格 - 2行4列 */}
                              <div style={{
                                display: 'grid',
                                gridTemplateColumns: 'repeat(4, 1fr)',
                                gap: '6px',
                                fontSize: '11px'
                              }}>
                                {mol.properties && Object.entries(mol.properties).slice(0, 8).map(([key, value]: [string, any]) => {
                                  const isQM9 = QM9_PROPERTIES.includes(key);
                                  const isRealValue = mol.is_real_data?.[key] || false;

                                  return (
                                    <div key={key} style={{
                                      padding: '4px 6px',
                                      background: isDark ? '#1a1a1a' : '#ffffff',
                                      borderRadius: '3px',
                                      border: `1px solid ${token.colorBorder}`,
                                      textAlign: 'center',
                                    }}>
                                      <div style={{ fontSize: '8px', color: token.colorTextSecondary, marginBottom: '2px' }}>
                                        {key.toUpperCase()}
                                      </div>
                                      <div style={{ fontWeight: 'bold', fontSize: '10px', color: isQM9 ? '#1890ff' : '#52c41a' }}>
                                        {typeof value === 'number' ? value.toFixed(2) : value}
                                      </div>
                                      <div style={{ fontSize: '7px', color: isRealValue ? '#10b981' : '#999' }}>
                                        {isRealValue ? '实测' : '预测'}
                                      </div>
                                    </div>
                                  );
                                })}
                              </div>
                            </div>

                            {/* 右侧：操作按钮 */}
                            <div style={{
                              width: '90px',
                              flexShrink: 0,
                              background: isDark ? '#1a1a1a' : '#f5f5f5',
                              borderLeft: `1px solid ${token.colorBorder}`,
                              display: 'flex',
                              flexDirection: 'column',
                              alignItems: 'center',
                              justifyContent: 'center',
                              gap: '6px',
                              padding: '8px',
                            }}>
                              <Button
                                size="small"
                                type="primary"
                                style={{
                                  fontSize: '11px',
                                  padding: '4px 8px',
                                  height: '28px',
                                  borderRadius: '4px',
                                  width: '100%',
                                  background: `linear-gradient(135deg, ${token.colorPrimary} 0%, ${token.colorPrimaryActive} 100%)`,
                                  border: 'none',
                                }}
                                onClick={() => jumpToOptimize(mol)}
                              >
                                优化
                              </Button>
                              <Button
                                size="small"
                                style={{
                                  fontSize: '11px',
                                  padding: '4px 8px',
                                  height: '28px',
                                  borderRadius: '4px',
                                  width: '100%',
                                  background: token.colorBgElevated,
                                  border: `1px solid ${token.colorBorder}`,
                                  color: token.colorText
                                }}
                                onClick={() => jumpToAlternatives(mol)}
                              >
                                找替代
                              </Button>
                            </div>
                          </div>
                        </div>
                      ))}
                    </div>
                  </div>
                )}
              </Card>
            </Col>
          </Row>
        )}

        {currentMode === 'alternatives' && (
          <Row gutter={24}>
            {/* 左侧：基准分子选择和目标设定 */}
            <Col xs={24} lg={10}>
              <Card title={
                <span>
                  <SearchOutlined style={{ marginRight: '8px' }} />
                  基准分子与目标
                </span>
              } style={cardStyle}>
                {/* 基准分子选择 */}
                <div style={{ marginBottom: '24px' }}>
                  <Text strong style={{ display: 'block', marginBottom: '12px' }}>
                    1. 选择基准分子
                  </Text>

                  {/* 从预测结果选择 */}
                  {predictResults.length > 0 && (
                    <div style={{ marginBottom: '16px' }}>
                      <Text style={{ fontSize: '12px', color: token.colorTextSecondary, display: 'block', marginBottom: '8px' }}>
                        📊 从预测结果中选择：
                      </Text>
                      <Row gutter={[8, 8]}>
                        {predictResults.slice(0, 4).map((mol, idx) => (
                          <Col key={idx} xs={12} sm={6}>
                            <Card
                              size="small"
                              hoverable
                              style={{
                                cursor: 'pointer',
                                border: taskSummary.baselineMolecule?.smiles === mol.smiles
                                  ? `2px solid ${token.colorPrimary}`
                                  : `1px solid ${token.colorBorder}`,
                              }}
                              onClick={() => updateTaskSummary(mol, [], [])}
                            >
                              <Text strong style={{ fontSize: '11px', display: 'block' }}>
                                {mol.name || `#${idx + 1}`}
                              </Text>
                              <Text style={{ fontSize: '10px', color: token.colorTextSecondary }}>
                                {mol.smiles.substring(0, 15)}...
                              </Text>
                            </Card>
                          </Col>
                        ))}
                      </Row>
                    </div>
                  )}

                  {/* 从模板选择 */}
                  {centerMolecules.length > 0 && (
                    <div>
                      <Text style={{ fontSize: '12px', color: token.colorTextSecondary, display: 'block', marginBottom: '8px' }}>
                        🔬 从数据库中选择：
                      </Text>
                      <Select
                        placeholder="选择基准分子"
                        style={{ width: '100%' }}
                        showSearch
                        optionFilterProp="children"
                        value={taskSummary.baselineMolecule?.smiles}
                        onChange={(value) => {
                          const mol = centerMolecules.find(m => m.smiles === value);
                          if (mol) updateTaskSummary(mol, [], []);
                        }}
                      >
                        {centerMolecules.slice(0, 40).map((mol: any, idx) => (
                          <Select.Option key={idx} value={mol.smiles}>
                            <span style={{ fontWeight: mol.isCommon ? 'bold' : 'normal' }}>
                              {mol.name || mol.smiles.substring(0, 30)}
                              {mol.isCommon && <Tag color="blue" style={{ marginLeft: '4px', fontSize: '10px' }}>常用</Tag>}
                            </span>
                          </Select.Option>
                        ))}
                      </Select>
                    </div>
                  )}
                </div>

                {/* 目标预设选择 */}
                {taskSummary.baselineMolecule && (
                  <div style={{ marginBottom: '24px' }}>
                    <Text strong style={{ display: 'block', marginBottom: '12px' }}>
                      2. 选择目标预设
                    </Text>
                    <Row gutter={[8, 8]}>
                      {GOAL_PRESETS.map((preset) => (
                        <Col key={preset.id} xs={24}>
                          <Card
                            size="small"
                            hoverable
                            style={{
                              cursor: 'pointer',
                              border: selectedPreset === preset.id
                                ? `2px solid ${token.colorPrimary}`
                                : `1px solid ${token.colorBorder}`,
                            }}
                            onClick={() => handlePresetSelect(preset.id)}
                          >
                            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
                              <div>
                                <Text strong style={{ display: 'block', marginBottom: '4px' }}>
                                  {preset.name}
                                </Text>
                                <Text style={{ fontSize: '12px', color: token.colorTextSecondary }}>
                                  {preset.description}
                                </Text>
                              </div>
                              {selectedPreset === preset.id && (
                                <StarOutlined style={{ color: token.colorPrimary }} />
                              )}
                            </div>
                          </Card>
                        </Col>
                      ))}
                    </Row>
                  </div>
                )}

                {/* 搜索按钮 */}
                {taskSummary.baselineMolecule && selectedPreset && (
                  <Button
                    type="primary"
                    size="large"
                    block
                    loading={alternativesLoading}
                    onClick={handleSearchAlternatives}
                  >
                    🔍 搜索替代分子
                  </Button>
                )}
              </Card>
            </Col>

            {/* 右侧：候选结果 */}
            <Col xs={24} lg={14}>
              <Card title={
                <span>
                  <TrophyOutlined style={{ marginRight: '8px' }} />
                  候选结果
                </span>
              } style={cardStyle}>
                {alternativesResults.length === 0 ? (
                  <Empty
                    description="请在左侧选择基准分子和目标，然后搜索"
                    style={{ marginTop: '40px', marginBottom: '40px' }}
                    image={Empty.PRESENTED_IMAGE_SIMPLE}
                  />
                ) : (
                  <div>
                    <div style={{ marginBottom: '16px', display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                      <Text>
                        找到 <Text strong>{alternativesResults.length}</Text> 个候选分子
                        （<Text style={{ color: '#52c41a' }}>{alternativesResults.filter(r => r.constraintsMet).length}</Text> 个达标）
                      </Text>
                    </div>

                    {/* 候选表格 */}
                    <div style={{ overflowX: 'auto' }}>
                      <table style={{
                        width: '100%',
                        borderCollapse: 'collapse',
                        fontSize: '12px',
                        background: token.colorBgContainer,
                        borderRadius: '8px',
                        overflow: 'hidden'
                      }}>
                        <thead>
                          <tr style={{ background: token.colorPrimaryBg, borderBottom: `1px solid ${token.colorBorder}` }}>
                            <th style={{ padding: '12px', textAlign: 'left', fontWeight: 600 }}>分子</th>
                            <th style={{ padding: '12px', textAlign: 'center', fontWeight: 600 }}>达标</th>
                            <th style={{ padding: '12px', textAlign: 'center', fontWeight: 600 }}>相似度</th>
                            <th style={{ padding: '12px', textAlign: 'center', fontWeight: 600 }}>Δgap</th>
                            <th style={{ padding: '12px', textAlign: 'center', fontWeight: 600 }}>Δfp</th>
                            <th style={{ padding: '12px', textAlign: 'center', fontWeight: 600 }}>Δbp</th>
                            <th style={{ padding: '12px', textAlign: 'center', fontWeight: 600 }}>得分</th>
                          </tr>
                        </thead>
                        <tbody>
                          {alternativesResults.slice(0, 20).map((candidate, idx) => (
                            <tr
                              key={idx}
                              style={{
                                borderBottom: `1px solid ${token.colorBorder}`,
                                background: idx % 2 === 0 ? token.colorBgContainer : token.colorBgElevated,
                              }}
                            >
                              <td style={{ padding: '12px', maxWidth: '150px' }}>
                                <Text strong style={{ display: 'block', fontSize: '11px' }}>
                                  {candidate.name || `分子 #${idx + 1}`}
                                </Text>
                                <Text style={{ fontSize: '10px', color: token.colorTextSecondary }}>
                                  {candidate.smiles.substring(0, 20)}...
                                </Text>
                              </td>
                              <td style={{ padding: '12px', textAlign: 'center' }}>
                                {candidate.constraintsMet ? (
                                  <Badge status="success" text="达标" />
                                ) : (
                                  <Badge status="default" text="未达标" />
                                )}
                              </td>
                              <td style={{ padding: '12px', textAlign: 'center' }}>
                                {candidate.distance ? (1 - candidate.distance).toFixed(3) : '-'}
                              </td>
                              <td style={{ padding: '12px', textAlign: 'center' }}>
                                {candidate.deltas?.gap !== undefined ? (
                                  <span style={{ color: candidate.deltas.gap > 0 ? '#52c41a' : '#ff4d4f' }}>
                                    {candidate.deltas.gap > 0 ? '↑' : '↓'} {Math.abs(candidate.deltas.gap).toFixed(2)}
                                  </span>
                                ) : '-'}
                              </td>
                              <td style={{ padding: '12px', textAlign: 'center' }}>
                                {candidate.deltas?.fp !== undefined ? (
                                  <span style={{ color: candidate.deltas.fp > 0 ? '#52c41a' : '#ff4d4f' }}>
                                    {candidate.deltas.fp > 0 ? '↑' : '↓'} {Math.abs(candidate.deltas.fp).toFixed(2)}
                                  </span>
                                ) : '-'}
                              </td>
                              <td style={{ padding: '12px', textAlign: 'center' }}>
                                {candidate.deltas?.bp !== undefined ? (
                                  <span style={{ color: candidate.deltas.bp > 0 ? '#52c41a' : '#ff4d4f' }}>
                                    {candidate.deltas.bp > 0 ? '↑' : '↓'} {Math.abs(candidate.deltas.bp).toFixed(2)}
                                  </span>
                                ) : '-'}
                              </td>
                              <td style={{ padding: '12px', textAlign: 'center' }}>
                                <Text strong style={{ color: token.colorPrimary }}>
                                  {candidate.objectiveScore.toFixed(2)}
                                </Text>
                              </td>
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    </div>
                  </div>
                )}
              </Card>
            </Col>
          </Row>
        )}

        {currentMode === 'optimize' && (
          <Row gutter={24}>
            {/* 左侧：优化目标设定 */}
            <Col xs={24} lg={10}>
              <Card title={
                <span>
                  <RocketOutlined style={{ marginRight: '8px' }} />
                  优化目标设定
                </span>
              } style={cardStyle}>
                <Space direction="vertical" style={{ width: '100%' }} size="large">
                  {/* 基准分子选择 */}
                  <div>
                    <Text strong style={{ display: 'block', marginBottom: '8px' }}>
                      基准分子
                    </Text>
                    {taskSummary.baselineMolecule ? (
                      <Card size="small" style={{ background: token.colorBgElevated }}>
                        <Text strong>{taskSummary.baselineMolecule.name || '未命名分子'}</Text>
                        <br />
                        <Text code style={{ fontSize: '11px' }}>{taskSummary.baselineMolecule.smiles}</Text>
                      </Card>
                    ) : (
                      <Alert message="请先在评估模式中预测分子，或在替代模式中选择基准分子" type="warning" />
                    )}
                  </div>

                  {/* 优化方法选择 */}
                  {taskSummary.baselineMolecule && (
                    <div>
                      <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', marginBottom: '16px' }}>
                        <Text strong>优化方法</Text>
                        <Segmented
                          value={optimizationMethod}
                          onChange={setOptimizationMethod}
                          options={[
                            { label: '📊 Pareto前沿', value: 'pareto' },
                            { label: '⚖️ 加权优化', value: 'weighted' }
                          ]}
                          size="small"
                        />
                      </div>

                      {optimizationMethod === 'pareto' && (
                        <Alert
                          message="📊 Pareto前沿: 找到所有不被支配的最优解集合，提供多种权衡选择，无需预设权重，适合探索性分析和决策支持"
                          type="info"
                          style={{ marginBottom: '16px', fontSize: '12px' }}
                        />
                      )}

                      {optimizationMethod === 'weighted' && (
                        <Alert
                          message="⚖️ 加权优化: 基于用户设定的权重找到单一最优解，适合有明确偏好的优化场景"
                          type="info"
                          style={{ marginBottom: '16px', fontSize: '12px' }}
                        />
                      )}
                    </div>
                  )}

                  {/* 多目标权重配置 */}
                  {taskSummary.baselineMolecule && (
                    <div>
                      <Text strong style={{ display: 'block', marginBottom: '12px' }}>
                        多目标权重配置
                      </Text>
                      <Text style={{ display: 'block', marginBottom: '12px', fontSize: '12px', color: token.colorTextSecondary }}>
                        与【找替代】不同，这里可以精确调节每个属性的重要性权重
                      </Text>

                      {/* 8种属性的4x2网格布局 */}
                      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(4, 1fr)', gap: '12px' }}>
                        {DEFAULT_PROPERTIES.map(prop => {
                          const isQM9 = QM9_PROPERTIES.includes(prop);
                          return (
                            <div key={prop} style={{
                              padding: '12px',
                              border: `1px solid ${token.colorBorder}`,
                              borderRadius: '8px',
                              background: token.colorBgElevated,
                              transition: 'all 0.2s'
                            }}>
                              <div style={{ display: 'flex', alignItems: 'center', marginBottom: '8px' }}>
                                <Text strong style={{ fontSize: '13px' }}>
                                  {prop.toUpperCase()}
                                </Text>
                                <Tag
                                  color={isQM9 ? 'blue' : 'green'}
                                  style={{ marginLeft: '4px', fontSize: '9px' }}
                                >
                                  {isQM9 ? 'QM9' : 'EXP'}
                                </Tag>
                              </div>

                              <Select
                                size="small"
                                value={objectiveWeights[prop] > 0 ? 'up' : objectiveWeights[prop] < 0 ? 'down' : 'none'}
                                onChange={(value) => {
                                  const newWeights = { ...objectiveWeights };
                                  if (value === 'none') {
                                    newWeights[prop] = 0;
                                  } else if (value === 'up') {
                                    newWeights[prop] = Math.abs(newWeights[prop]) || 0.5;
                                  } else if (value === 'down') {
                                    newWeights[prop] = -(Math.abs(newWeights[prop]) || 0.5);
                                  }
                                  setObjectiveWeights(newWeights);
                                }}
                                style={{ width: '100%', marginBottom: '8px' }}
                                options={[
                                  { value: 'none', label: '不优化' },
                                  { value: 'up', label: '提高 ↑' },
                                  { value: 'down', label: '降低 ↓' }
                                ]}
                              />

                              <Slider
                                min={0}
                                max={1}
                                step={0.05}
                                value={Math.abs(objectiveWeights[prop])}
                                onChange={(value) => {
                                  const newWeights = { ...objectiveWeights };
                                  const currentDirection = newWeights[prop] >= 0 ? 1 : -1;
                                  newWeights[prop] = value * currentDirection;
                                  setObjectiveWeights(newWeights);
                                }}
                                marks={{ 0: '0', 1: '1' }}
                                disabled={objectiveWeights[prop] === 0}
                                style={{ margin: '8px 0 4px 0' }}
                              />

                              <Text style={{ fontSize: '11px', color: token.colorTextTertiary, textAlign: 'center', display: 'block' }}>
                                权重: {Math.abs(objectiveWeights[prop]).toFixed(2)}
                              </Text>
                            </div>
                          );
                        })}
                      </div>
                    </div>
                  )}

                  {/* 开始优化按钮 */}
                  {taskSummary.baselineMolecule && (
                    <div>
                      <Button
                        type="primary"
                        size="large"
                        block
                        loading={multiObjectiveLoading}
                        onClick={handleMultiObjectiveOptimization}
                        style={{
                          background: `linear-gradient(135deg, ${token.colorSuccess} 0%, ${token.colorSuccessActive} 100%)`,
                          border: 'none',
                          borderRadius: '8px',
                          height: '48px',
                          fontSize: '16px',
                          fontWeight: 'bold',
                          boxShadow: `0 4px 12px ${token.colorSuccess}20`
                        }}
                      >
                        🎯 运行多目标优化
                      </Button>
                      <Text style={{ display: 'block', textAlign: 'center', marginTop: '8px', fontSize: '12px', color: token.colorTextTertiary }}>
                        将基于{optimizationMethod === 'pareto' ? 'Pareto前沿' : '权重配置'}寻找最优解
                      </Text>
                    </div>
                  )}
                </Space>
              </Card>
            </Col>

            {/* 右侧：优化结果 */}
            <Col xs={24} lg={14}>
              <Card title={
                <span>
                  <LineChartOutlined style={{ marginRight: '8px' }} />
                  优化结果
                </span>
              } style={cardStyle}>
                {!taskSummary.baselineMolecule ? (
                  <Empty
                    description="请先选择基准分子并设置优化目标"
                    style={{ marginTop: '40px', marginBottom: '40px' }}
                    image={Empty.PRESENTED_IMAGE_SIMPLE}
                  />
                ) : multiObjectiveResults.length === 0 ? (
                  <div style={{ textAlign: 'center', padding: '40px 20px' }}>
                    <div style={{ fontSize: '48px', marginBottom: '16px' }}>📊</div>
                    <Text style={{ fontSize: '16px', color: token.colorTextSecondary, display: 'block', marginBottom: '8px' }}>
                      多目标优化结果
                    </Text>
                    <Text style={{ fontSize: '14px', color: token.colorTextTertiary, display: 'block', marginBottom: '16px' }}>
                      将显示Pareto前沿、权衡分析、最优解集合
                    </Text>
                    <div style={{
                      background: token.colorBgElevated,
                      padding: '16px',
                      borderRadius: '8px',
                      border: `1px dashed ${token.colorBorder}`
                    }}>
                      <Text style={{ fontSize: '12px', color: token.colorTextTertiary }}>
                        💡 <strong>与"找替代"的区别：</strong><br />
                        • 找替代：基于相似度 + 简单目标筛选<br />
                        • 多目标优化：精确权重配置 + Pareto最优解
                      </Text>
                    </div>
                  </div>
                ) : (
                  <MultiObjectiveDashboard
                    optimizeResults={multiObjectiveResults}
                    paretoFrontier={multiObjectiveResults.filter(r => r.pareto_rank === 1)}
                    objectiveWeights={Object.fromEntries(
                      Object.entries(objectiveWeights).map(([key, weight]) => [
                        key,
                        { direction: weight > 0 ? 'higher' : 'lower', weight }
                      ])
                    )}
                    baselineMolecule={taskSummary.baselineMolecule}
                  />
                )}
              </Card>
            </Col>
          </Row>
        )}
      </Spin>
    </div>
  );
};

export default AIDiscovery;

