/**
 * 任务配置和创建页面
 */
import { useState, useEffect } from 'react';
import { useNavigate, useParams, useLocation } from 'react-router-dom';
import {
  Card,
  Form,
  Input,
  InputNumber,
  Button,
  Space,
  Typography,
  Descriptions,
  message,
  Alert,
  Divider,
  Row,
  Col,
  Select,
  Checkbox,
  Switch,
  theme,
  Tooltip,
  Tag,
  Affix,
  Collapse,
  Progress,
} from 'antd';
import {
  ArrowLeftOutlined,
  ThunderboltOutlined,
  SaveOutlined,
  WarningOutlined,
  ExperimentOutlined,
  QuestionCircleOutlined,
  RocketOutlined,
  StarOutlined,
  FireOutlined,
  CloudServerOutlined,
  SettingOutlined,
  InfoCircleOutlined,
} from '@ant-design/icons';
import { UserRole, type ElectrolyteSystem, type MDJobCreate } from '../types';
import { getElectrolyte } from '../api/electrolytes';
import { createMDJob, checkJobQuota } from '../api/jobs';
import { getPartitions, type PartitionInfo } from '../api/slurm';
import AccuracyLevelSelector from '../components/AccuracyLevelSelector';
import { useThemeStore } from '../stores/themeStore';
import { useAuthStore } from '../stores/authStore';

const { Title, Text } = Typography;

export default function JobCreate() {
  const navigate = useNavigate();
  const location = useLocation();
  const { systemId } = useParams<{ systemId: string }>();
  const { mode } = useThemeStore();
  const { user } = useAuthStore();
  const { token } = theme.useToken();
  const isDark = mode === 'dark';
  const [form] = Form.useForm();
  const [loading, setLoading] = useState(false);
  const [electrolyte, setElectrolyte] = useState<ElectrolyteSystem | null>(null);
  const [submitting, setSubmitting] = useState(false);
  const [quota, setQuota] = useState<{
    can_create: boolean;
    current_count: number;
    limit: number;
    remaining: number;
  } | null>(null);
  const [accuracyDefaults, setAccuracyDefaults] = useState<any>(null);
  const [selectedAccuracyLevel, setSelectedAccuracyLevel] = useState<string>('standard');  // QC 智能推荐开关
  const [useSmartRecommend, setUseSmartRecommend] = useState(true);

  // 集群分区列表
  const [partitions, setPartitions] = useState<PartitionInfo[]>([]);

  // QCConfig 类型定义
  interface QCConfig {
    qc_engine?: string;
    use_xtb_preopt?: boolean;
    qc_solvent_models?: string[];
    qc_solvents?: string[];
    qc_functionals?: string[];
    qc_basis_sets?: string[];
  }

  // QC 配置
  const [qcConfig, setQcConfig] = useState<Partial<QCConfig>>({});




  // 监听配方变化，如果开启了智能推荐，则自动设置参数
  useEffect(() => {
    if (useSmartRecommend && electrolyte) {
      const name = electrolyte.name.toLowerCase();
      let sol = 'Acetone'; // general organic
      if (name.includes('water') || name.includes('aq') || name.includes('h2o')) sol = 'Water';
      else if (name.includes('dmc') || name.includes('emc') || name.includes('dec')) sol = 'Chloroform';

      // 避免重复设置导致表单脏检查问题，这里简单直接设置
      form.setFieldsValue({
        qc_solvent_models: ['pcm'],
        qc_solvents: [sol],
        qc_functionals: ['B3LYP'],
        qc_basis_sets: ['6-31++g(d,p)']
      });
    }
  }, [electrolyte, useSmartRecommend, form]);

  // 从 location.state 或 API 加载电解质系统信息
  useEffect(() => {
    const loadElectrolyte = async () => {
      setLoading(true);
      try {
        // 优先使用 location.state 传递的数据
        if (location.state?.electrolyte) {
          setElectrolyte(location.state.electrolyte);
        } else if (systemId) {
          const data = await getElectrolyte(Number(systemId));
          setElectrolyte(data);
        }
      } catch (error: any) {
        message.error('加载配方信息失败: ' + (error.response?.data?.detail || error.message));
        navigate('/workspace/liquid-electrolyte/electrolytes');
      } finally {
        setLoading(false);
      }
    };

    loadElectrolyte();
  }, [systemId, location.state]);

  // 加载配额信息
  useEffect(() => {
    const loadQuota = async () => {
      try {
        const quotaData = await checkJobQuota();
        setQuota(quotaData);
      } catch (error: any) {
        console.error('加载配额信息失败:', error);
      }
    };

    loadQuota();
  }, []);

  // 加载集群分区信息
  useEffect(() => {
    const loadPartitions = async () => {
      try {
        const data = await getPartitions();
        setPartitions(data);
      } catch (error) {
        console.error('Failed to load partitions:', error);
        // 如果加载失败，使用默认值
        setPartitions([]);
      }
    };

    loadPartitions();
  }, []);

  // 加载精度等级配置
  useEffect(() => {
    const loadAccuracyLevels = async () => {
      try {
        const token = localStorage.getItem('access_token');
        const response = await fetch('/api/v1/jobs/accuracy-levels', {
          headers: { Authorization: `Bearer ${token}` }
        });
        const data = await response.json();
        setAccuracyDefaults(data);
      } catch (error) {
        console.error('加载精度等级配置失败:', error);
      }
    };
    loadAccuracyLevels();
  }, []);

  // 获取默认分区
  const getDefaultPartition = () => {
    if (partitions.length > 0) {
      const upPartition = partitions.find(p => p.state === 'up');
      return upPartition?.name || partitions[0].name;
    }
    return 'cpu';
  };

  // 设置默认值
  useEffect(() => {
    if (electrolyte) {
      form.setFieldsValue({
        job_name: '', // 留空，让后端自动生成
        accuracy_level: 'standard', // 默认标准模式
        // 所有参数留空，使用精度等级的默认值
        nsteps_npt: undefined,
        nsteps_nvt: undefined,
        timestep: undefined,
        temperature: 298.15, // 默认常温
        pressure: 1.0,       // 默认常压
        freq_trj_npt: undefined,
        freq_trj_nvt: undefined,
        thermo_freq: undefined,
        // Slurm 资源配置默认值
        slurm_partition: getDefaultPartition(),
        slurm_nodes: 1,
        slurm_ntasks: 8,
        slurm_cpus_per_task: 8,
        slurm_time: 7200,
      });
    }
  }, [electrolyte, form, partitions]);

  // 获取当前精度等级的默认值
  const getCurrentDefaults = () => {
    return accuracyDefaults?.[selectedAccuracyLevel] || {};
  };

  // 检查是否为自定义模式
  const isCustomMode = () => {
    return selectedAccuracyLevel === 'custom';
  };

  // 将步数转换为时间，自动选择合适的单位（ns、ps或fs）
  const stepsToTime = (steps: number | null | undefined, timestep: number = 1.0): string => {
    if (!steps) return '-';

    const timeFs = steps * timestep; // 总时间（飞秒）

    // 根据时间大小选择合适的单位
    if (timeFs >= 1_000_000) {
      // >= 1 ns，使用 ns
      const timeNs = timeFs / 1_000_000;
      return `${timeNs.toFixed(1)} ns`;
    } else if (timeFs >= 1_000) {
      // >= 1 ps，使用 ps
      const timePs = timeFs / 1_000;
      return `${timePs.toFixed(1)} ps`;
    } else {
      // < 1 ps，使用 fs
      return `${timeFs.toFixed(0)} fs`;
    }
  };

  // 格式化步数显示（步数 + 时间）
  const formatStepsWithTime = (steps: number | null | undefined, timestep: number = 1.0): string => {
    if (!steps) return '-';
    return `${steps.toLocaleString()} 步 (${stepsToTime(steps, timestep)})`;
  };

  // 生成字段配置的辅助函数
  const getFieldConfig = (fieldName: string, label: string, unit: string = '') => {
    const customMode = isCustomMode();
    const defaultValue = getCurrentDefaults()[fieldName];

    return {
      label: `${label}${customMode ? '（可修改）' : '（可选）'}`,
      rules: customMode ? [{ required: true, message: `请输入${label}` }] : [],
      tooltip: customMode
        ? '已填充参考值，您可以根据需要修改'
        : '留空使用精度等级的默认值',
      extra: customMode && defaultValue
        ? `参考值（标准模式）: ${typeof defaultValue === 'number' ? defaultValue.toLocaleString() : defaultValue}${unit}`
        : (!customMode && accuracyDefaults && defaultValue
          ? `默认值: ${typeof defaultValue === 'number' ? defaultValue.toLocaleString() : defaultValue}${unit}`
          : undefined),
      placeholder: customMode
        ? `参考值: ${defaultValue || ''}${unit}`
        : (defaultValue ? `默认: ${typeof defaultValue === 'number' ? defaultValue.toLocaleString() : defaultValue}${unit}` : '请选择精度等级'),
    };
  };

  // 构建任务数据
  const buildJobData = (values: any, submitToCluster: boolean = false): MDJobCreate => {
    const jobData: MDJobCreate = {
      system_id: electrolyte!.id,
      job_name: values.job_name || undefined,
      accuracy_level: values.accuracy_level || 'standard',
      charge_method: values.charge_method || undefined,  // 电荷计算方法（仅自定义模式有效）
      nsteps_npt: values.nsteps_npt || undefined,
      nsteps_nvt: values.nsteps_nvt || undefined,
      timestep: values.timestep,
      temperature: values.temperature,
      pressure: values.pressure,
      freq_trj_npt: values.freq_trj_npt || undefined,
      freq_trj_nvt: values.freq_trj_nvt || undefined,
      thermo_freq: values.thermo_freq || undefined,
      submit_to_cluster: submitToCluster,
      // Slurm 资源配置 - 从简化的slurm_cpus转换
      slurm_partition: values.slurm_partition || 'cpu',
      slurm_nodes: 1, // 固定1个节点
      slurm_ntasks: 8, // 固定8个任务  
      slurm_cpus_per_task: values.slurm_cpus ? Math.ceil(values.slurm_cpus / 8) : 2, // 根据总核心数计算
      slurm_time: values.slurm_time || 7200,
      // ECC 配置
      use_ecc: values.use_ecc,
      ecc_factor: values.ecc_factor,
    };

    // QC计算选项 - 支持多选
    if (values.qc_enabled) {
      jobData.qc_options = {
        enabled: true,
        // 使用复数形式的数组字段（后端已支持）
        functionals: values.qc_functionals || ['B3LYP'],
        basis_sets: values.qc_basis_sets || ['6-31++g(d,p)'],
        solvent_models: values.qc_solvent_models || ['pcm'],
        solvents: values.qc_solvents || ['Water'],
        molecules: [], // 将由后端从电解质配方中提取
        // 兼容旧版字段（取第一个值）
        functional: values.qc_functionals?.[0] || 'B3LYP',
        basis_set: values.qc_basis_sets?.[0] || '6-31++g(d,p)',
        solvent_model: values.qc_solvent_models?.[0] || 'pcm',
        solvent_name: values.qc_solvents?.[0] || 'Water',
        // QC引擎
        qc_engine: values.qc_engine || 'pyscf',
        // 自定义溶剂参数（如果选择了custom溶剂模型）
        custom_solvent: values.custom_solvent || undefined,
      } as any;
    }

    return jobData;
  };

  // 通用提交处理
  const handleJobSubmit = async (submitToCluster: boolean) => {
    try {
      const values = await form.validateFields();
      setSubmitting(true);

      const jobData = buildJobData(values, submitToCluster);
      const job = await createMDJob(jobData);

      const successMsg = submitToCluster ? '计算任务已提交到集群！' : '计算任务已保存！';
      const targetPath = submitToCluster ? `/workspace/liquid-electrolyte/md/${job.id}` : `/workspace/liquid-electrolyte/md/${job.id}/submit`;

      message.success(successMsg);
      navigate(targetPath);
    } catch (error: any) {
      if (error.errorFields) {
        message.error('请检查表单填写');
      } else {
        const action = submitToCluster ? '提交' : '保存';
        message.error(`${action}任务失败: ` + (error.response?.data?.detail || error.message));
      }
    } finally {
      setSubmitting(false);
    }
  };

  // 保存任务（不提交到集群）
  const handleSave = () => handleJobSubmit(false);

  // 提交任务到集群
  const handleSubmit = () => handleJobSubmit(true);

  if (loading || !electrolyte) {
    return (
      <div style={{
        display: 'flex',
        justifyContent: 'center',
        alignItems: 'center',
        height: 'calc(100vh - 64px)',
        background: token.colorBgLayout,
      }}>
        加载中...
      </div>
    );
  }

  return (
    <Form form={form} layout="vertical" onFinish={handleSubmit}>
      <div style={{ padding: 24, background: token.colorBgLayout, minHeight: 'calc(100vh - 64px)', transition: 'background 0.3s' }}>
        {/* 页面头部 */}
        <div style={{ marginBottom: 24 }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', marginBottom: 16 }}>
            <div>
              <Button
                icon={<ArrowLeftOutlined />}
                onClick={() => navigate('/workspace/electrolytes')}
                style={{ marginBottom: 12 }}
              >
                返回
              </Button>
              <Title level={2} style={{ margin: 0, marginBottom: 4 }}>
                <ThunderboltOutlined style={{ marginRight: 12, color: '#1677ff' }} />
                创建计算任务
              </Title>
              <Text type="secondary">配置分子动力学模拟参数并提交计算任务</Text>
            </div>
            {quota && (
              <Tag color={quota.can_create ? 'blue' : 'red'} style={{ fontSize: 14, padding: '6px 16px', height: 'auto' }}>
                今日配额: {quota.current_count} / {quota.limit}
              </Tag>
            )}
          </div>

          {/* 配额警告 */}
          {quota && !quota.can_create && (
            <Alert
              message="配额已用完"
              description="您今日的任务创建配额已用完，请明天再试或联系管理员。"
              type="warning"
              showIcon
              style={{ borderRadius: 8 }}
            />
          )}
        </div>

        <Row gutter={24}>
          {/* 左侧主要配置区 */}
          <Col xs={24} lg={14}>
            {/* 基础信息卡片 */}

            {/* 核心参数配置 */}
            <Card
              title="核心参数"
              style={{ marginBottom: 24 }}
            >
              <Form.Item name="accuracy_level">
                <AccuracyLevelSelector
                  value={selectedAccuracyLevel}
                  onChange={(value) => {
                    setSelectedAccuracyLevel(value);
                    form.setFieldsValue({ accuracy_level: value });

                    const defaults = accuracyDefaults?.[value];
                    if (defaults) {
                      // 无论什么模式，都填充默认值（解决用户反馈选了模式不填充的问题）
                      // 对于自定义模式，这提供了初始编辑值
                      // 对于标准模式，这确保了显示的数值正确
                      form.setFieldsValue({
                        nsteps_npt: defaults.nsteps_npt,
                        nsteps_nvt: defaults.nsteps_nvt,
                        timestep: defaults.timestep || 1.0,
                        // 温度和压力：如果配置中有默认值则填充，否则保持 298.15 / 1.0
                        temperature: defaults.temperature || 298.15,
                        pressure: defaults.pressure || 1.0,
                        freq_trj_npt: defaults.freq_trj_npt,
                        freq_trj_nvt: defaults.freq_trj_nvt,
                        thermo_freq: defaults.thermo_freq,
                      });
                    }
                  }}
                />
              </Form.Item>

              {isCustomMode() && (
                <Alert
                  message="自定义模式"
                  description="下方已填充参考值，您可根据需要修改。"
                  type="warning"
                  showIcon
                  style={{ marginBottom: 16 }}
                />
              )}

              {/* 模拟参数配置 - 统一放在左侧主栏 */}
              {/* 模拟参数配置 - 统一放在左侧主栏 */}
              <Text type="secondary" style={{ fontSize: 12, marginBottom: 8, marginTop: 24, display: 'block' }}>电荷与结构</Text>
              <Row gutter={16}>
                <Col span={8}>
                  <Form.Item
                    name="charge_method"
                    label="电荷产生方式"
                    initialValue={getCurrentDefaults().charge_method || 'ligpargen'}
                    rules={[{ required: true, message: '请选择电荷产生方式' }]}
                  >
                    {isCustomMode() ? (
                      <Select>
                        <Select.Option value="ligpargen">🚀 LigParGen</Select.Option>
                        <Select.Option value="resp">🎯 RESP</Select.Option>
                      </Select>
                    ) : (
                      <Text>{getCurrentDefaults().charge_method === 'resp' ? '🎯 RESP (高精度)' : '🚀 LigParGen (快速)'}</Text>
                    )}
                  </Form.Item>
                </Col>
                <Col span={8}>
                  <Form.Item
                    name="nsteps_npt"
                    label="NPT 平衡"
                    tooltip="恒温恒压系综，用于调整密度"
                    initialValue={getCurrentDefaults().nsteps_npt}
                  >
                    {isCustomMode() ? (
                      <InputNumber min={1000} step={100000} style={{ width: '100%' }} addonAfter="步" />
                    ) : (
                      <Text strong>{formatStepsWithTime(getCurrentDefaults().nsteps_npt, getCurrentDefaults().timestep)}</Text>
                    )}
                  </Form.Item>
                </Col>
                <Col span={8}>
                  <Form.Item
                    name="nsteps_nvt"
                    label="NVT 采样"
                    tooltip="恒温恒体积系综，用于生产采样"
                    initialValue={getCurrentDefaults().nsteps_nvt}
                  >
                    {isCustomMode() ? (
                      <InputNumber min={1000} step={100000} style={{ width: '100%' }} addonAfter="步" />
                    ) : (
                      <Text strong>{formatStepsWithTime(getCurrentDefaults().nsteps_nvt, getCurrentDefaults().timestep)}</Text>
                    )}
                  </Form.Item>
                </Col>
              </Row>

              <Text type="secondary" style={{ fontSize: 12, marginBottom: 8, marginTop: 16, display: 'block' }}>物理条件</Text>
              <Row gutter={16}>
                <Col span={8}>
                  <Form.Item
                    name="temperature"
                    label="温度"
                    rules={[{ required: true, message: '请输入温度' }]}
                    initialValue={getCurrentDefaults().temperature || 298.15}
                  >
                    <InputNumber min={200} max={500} step={1} style={{ width: '100%' }} addonAfter="K" />
                  </Form.Item>
                </Col>
                <Col span={8}>
                  <Form.Item
                    name="pressure"
                    label="压力"
                    rules={[{ required: true, message: '请输入压力' }]}
                    initialValue={getCurrentDefaults().pressure || 1.0}
                  >
                    <InputNumber min={0.1} max={100} step={0.1} style={{ width: '100%' }} addonAfter="atm" />
                  </Form.Item>
                </Col>
                <Col span={8}>
                  <Form.Item
                    name="timestep"
                    label="步长"
                    initialValue={getCurrentDefaults().timestep || 1.0}
                  >
                    {isCustomMode() ? (
                      <InputNumber min={0.1} max={10} step={0.1} style={{ width: '100%' }} addonAfter="fs" />
                    ) : (
                      <Text>{getCurrentDefaults().timestep} fs</Text>
                    )}
                  </Form.Item>
                </Col>
              </Row>

              <Text type="secondary" style={{ fontSize: 12, marginBottom: 8, marginTop: 16, display: 'block' }}>输出频率 (每多少步)</Text>
              <Row gutter={16}>
                <Col span={8}>
                  <Form.Item
                    name="freq_trj_npt"
                    label="NPT 轨迹"
                    initialValue={getCurrentDefaults().freq_trj_npt}
                  >
                    {isCustomMode() ? (
                      <InputNumber min={100} step={1000} style={{ width: '100%' }} />
                    ) : (
                      <Text>{getCurrentDefaults().freq_trj_npt}</Text>
                    )}
                  </Form.Item>
                </Col>
                <Col span={8}>
                  <Form.Item
                    name="freq_trj_nvt"
                    label="NVT 轨迹"
                    initialValue={getCurrentDefaults().freq_trj_nvt}
                  >
                    {isCustomMode() ? (
                      <InputNumber min={100} step={1000} style={{ width: '100%' }} />
                    ) : (
                      <Text>{getCurrentDefaults().freq_trj_nvt}</Text>
                    )}
                  </Form.Item>
                </Col>
                <Col span={8}>
                  <Form.Item
                    name="thermo_freq"
                    label="热力学"
                    initialValue={getCurrentDefaults().thermo_freq}
                  >
                    {isCustomMode() ? (
                      <InputNumber min={100} step={1000} style={{ width: '100%' }} />
                    ) : (
                      <Text>{getCurrentDefaults().thermo_freq}</Text>
                    )}
                  </Form.Item>
                </Col>
              </Row>


              {/* 高级设置放在右侧栏 */}
            </Card>
          </Col>

          {/* 右侧辅助配置区 */}
          <Col xs={24} lg={10}>
            {/* 基础信息卡片 - 移至右侧以平衡高度 */}
            <Card
              title={<Space><InfoCircleOutlined /> 基础信息</Space>}
              style={{ marginBottom: 24 }}
            >
              <Descriptions column={1} size="small" style={{ marginBottom: 24 }}>
                <Descriptions.Item label="配方名称">{electrolyte.name}</Descriptions.Item>
                <Descriptions.Item label="配方 ID">#{electrolyte.id}</Descriptions.Item>
                <Descriptions.Item label="盒子大小">
                  {electrolyte.box_size ? Number(electrolyte.box_size).toFixed(1) : '-'} Å
                </Descriptions.Item>
                <Descriptions.Item label="力场">{electrolyte.force_field || 'OPLS-AA'}</Descriptions.Item>
              </Descriptions>

              <Form.Item
                label="任务备注"
                name="job_name"
                tooltip="用于记录任务目的或特殊说明"
                style={{ marginBottom: 0 }}
              >
                <Input placeholder="可选备注" allowClear />
              </Form.Item>
            </Card>

            <Card title="高级设置" style={{ marginBottom: 24 }}>
              <Text strong style={{ marginBottom: 12, display: 'block' }}>ECC 电荷缩放</Text>

              <Row gutter={24}>
                <Col span={6}>
                  <Form.Item
                    label="启用 ECC"
                    name="use_ecc"
                    valuePropName="checked"
                    initialValue={false}
                    tooltip="Electronic Continuum Correction - 补偿电子极化效应"
                  >
                    <Switch />
                  </Form.Item>
                </Col>

                <Col span={18}>
                  <Form.Item
                    noStyle
                    shouldUpdate={(prevValues, currentValues) => prevValues.use_ecc !== currentValues.use_ecc}
                  >
                    {({ getFieldValue }) =>
                      getFieldValue('use_ecc') && (
                        <Form.Item
                          label="缩放系数"
                          name="ecc_factor"
                          initialValue={0.8}
                        >
                          <Select style={{ width: '100%' }}>
                            <Select.Option value={0.75}>0.75 - 强极化</Select.Option>
                            <Select.Option value={0.8}>0.80 - 推荐</Select.Option>
                            <Select.Option value={0.85}>0.85 - 弱极化</Select.Option>
                          </Select>
                        </Form.Item>
                      )
                    }
                  </Form.Item>
                </Col>
              </Row>



              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: 12 }}>
                <Text strong>QC 量子化学计算</Text>
                <Space>
                  <Text type="secondary" style={{ fontSize: 12 }}>智能参数推荐</Text>
                  <Switch
                    checked={useSmartRecommend}
                    onChange={(checked) => {
                      setUseSmartRecommend(checked);
                      if (checked && electrolyte) {
                        const name = electrolyte.name.toLowerCase();
                        let sol = 'Acetone'; // general organic
                        if (name.includes('water') || name.includes('aq') || name.includes('h2o')) sol = 'Water';
                        else if (name.includes('dmc') || name.includes('emc') || name.includes('dec')) sol = 'Chloroform'; // low dielectric
                        form.setFieldsValue({
                          qc_solvent_models: ['pcm'],
                          qc_solvents: [sol],
                          qc_functionals: ['B3LYP'],
                          qc_basis_sets: ['6-31++g(d,p)']
                        });
                        message.success(`已启用智能推荐: ${sol} / B3LYP / 6-31++G**`);
                      }
                    }}
                  />
                </Space>
              </div>

              <Form.Item
                label="启用 QC 计算"
                name="qc_enabled"
                valuePropName="checked"
                initialValue={false}
                tooltip="获取 HOMO/LUMO 等电子结构性质"
              >
                <Switch />
              </Form.Item>

              <Form.Item noStyle shouldUpdate={(prevValues, currentValues) =>
                prevValues.qc_enabled !== currentValues.qc_enabled
              }>
                {({ getFieldValue }) => {
                  const qcEnabled = getFieldValue('qc_enabled');
                  if (!qcEnabled) return null;

                  // 收集将要计算的分子列表
                  const moleculesToCalc: Array<{ name: string, smiles: string, type: string, charge: number }> = [];
                  if (electrolyte) {
                    // 溶剂分子
                    electrolyte.solvents?.forEach(sol => {
                      if (sol.smiles && !moleculesToCalc.find(m => m.smiles === sol.smiles)) {
                        moleculesToCalc.push({ name: sol.name, smiles: sol.smiles, type: 'solvent', charge: 0 });
                      }
                    });
                    // 阳离子
                    electrolyte.cations?.forEach(cat => {
                      if (cat.smiles && !moleculesToCalc.find(m => m.smiles === cat.smiles)) {
                        moleculesToCalc.push({ name: cat.name, smiles: cat.smiles, type: 'cation', charge: 1 });
                      }
                    });
                    // 阴离子
                    electrolyte.anions?.forEach(an => {
                      if (an.smiles && !moleculesToCalc.find(m => m.smiles === an.smiles)) {
                        moleculesToCalc.push({ name: an.name, smiles: an.smiles, type: 'anion', charge: -1 });
                      }
                    });
                  }

                  return (
                    <>
                      {/* QC 引擎选择 */}
                      <Card style={{ marginBottom: 16 }}>
                        <Form.Item
                          name="qc_engine"
                          label={<Text strong>计算引擎</Text>}
                          initialValue="pyscf"
                          tooltip={
                            <div>
                              <p><strong>PySCF</strong>: 开源软件，适合大多数计算</p>
                              <p><strong>Gaussian</strong>: 商业软件，需要License授权</p>
                            </div>
                          }
                          style={{ marginBottom: 0 }}
                        >
                          <Select
                            style={{ width: '100%' }}
                            placeholder="选择计算引擎"
                          >
                            <Select.Option value="pyscf">PySCF (开源软件)</Select.Option>
                            {(user?.role === UserRole.ADMIN || user?.can_use_gaussian) && (
                              <Select.Option value="gaussian">Gaussian (商业软件)</Select.Option>
                            )}
                          </Select>
                        </Form.Item>
                      </Card>

                      {/* 计算参数配置 */}
                      <Card style={{ marginBottom: 16 }}>
                        <Row gutter={16}>
                          <Col span={12}>
                            <Form.Item
                              name="qc_functionals"
                              label="泛函"
                              initialValue={['B3LYP']}
                              style={{ marginBottom: 8 }}
                              tooltip="可选择多个泛函进行对比计算"
                            >
                              <Select mode="multiple" placeholder="选择泛函（可多选）">
                                <Select.Option value="B3LYP">B3LYP (混合泛函)</Select.Option>
                                <Select.Option value="M062X">M06-2X (Minnesota泛函)</Select.Option>
                                <Select.Option value="wB97XD">ωB97X-D (长程校正)</Select.Option>
                                <Select.Option value="PBE0">PBE0 (混合GGA)</Select.Option>
                                <Select.Option value="CAM-B3LYP">CAM-B3LYP (长程校正)</Select.Option>
                              </Select>
                            </Form.Item>
                          </Col>
                          <Col span={12}>
                            <Form.Item
                              name="qc_basis_sets"
                              label="基组"
                              initialValue={['6-31++g(d,p)']}
                              style={{ marginBottom: 8 }}
                              tooltip="可选择多个基组进行对比计算"
                            >
                              <Select mode="multiple" placeholder="选择基组（可多选）">
                                <Select.Option value="6-31g(d,p)">6-31G(d,p) (标准)</Select.Option>
                                <Select.Option value="6-31++g(d,p)">6-31++G(d,p) (含弥散)</Select.Option>
                                <Select.Option value="6-311g(d,p)">6-311G(d,p) (三重劈裂)</Select.Option>
                                <Select.Option value="6-311++g(d,p)">6-311++G(d,p) (三重劈裂+弥散)</Select.Option>
                                <Select.Option value="Def2TZVP">Def2-TZVP (高精度)</Select.Option>
                              </Select>
                            </Form.Item>
                          </Col>
                        </Row>

                        <Row gutter={16}>
                          <Col span={12}>
                            <Form.Item
                              name="qc_solvent_models"
                              label="溶剂环境"
                              initialValue={['pcm']}
                              style={{ marginBottom: 8 }}
                              tooltip={
                                <div>
                                  <p><strong>气相 (Gas)</strong>: 真空环境，无溶剂效应</p>
                                  <p><strong>PCM</strong>: 极化连续介质模型，使用介电常数描述溶剂</p>
                                  <p><strong>SMD</strong>: 溶剂密度模型，更精确但计算量更大</p>
                                  <p><strong>自定义</strong>: 手动设置介电常数等参数</p>
                                  <p>可多选进行对比计算</p>
                                </div>
                              }
                            >
                              <Select mode="multiple" placeholder="选择溶剂环境（可多选）">
                                <Select.Option value="gas">气相 (Gas Phase) - 无溶剂效应</Select.Option>
                                <Select.Option value="pcm">PCM - 极化连续介质模型</Select.Option>
                                <Select.Option value="smd">SMD - 溶剂密度模型（更精确）</Select.Option>
                                <Select.Option value="custom">自定义溶剂参数</Select.Option>
                              </Select>
                            </Form.Item>
                          </Col>
                          <Col span={12}>
                            <Form.Item noStyle shouldUpdate={(prevValues, currentValues) =>
                              prevValues.qc_solvent_models !== currentValues.qc_solvent_models
                            }>
                              {({ getFieldValue }) => {
                                const solventModels = getFieldValue('qc_solvent_models') || [];
                                const hasNonCustomModel = solventModels.some((m: string) => m !== 'gas' && m !== 'custom');
                                if (!hasNonCustomModel) return null;
                                return (
                                  <Form.Item
                                    name="qc_solvents"
                                    label="隐式溶剂"
                                    initialValue={['Water']}
                                    style={{ marginBottom: 8 }}
                                    tooltip={
                                      <div>
                                        <p><strong>选择原则</strong>：选择介电常数(ε)接近您电解液的溶剂</p>
                                        <hr style={{ margin: '4px 0', borderColor: 'rgba(255,255,255,0.3)' }} />
                                        <p>• <strong>水系电解液</strong>: 选择 Water (ε=78.4)</p>
                                        <p>• <strong>高浓电解液</strong>: 选择 Acetone (ε=20.5)</p>
                                        <p>• <strong>EC基电解液</strong>: 选择 Water 或 PC (ε≈65-90)</p>
                                        <p>• <strong>DMC/EMC/DEC电解液</strong>: 选择 Chloroform (ε≈3-5)</p>
                                        <p>• <strong>离子液体</strong>: 选择 DMSO (ε=46.8)</p>
                                      </div>
                                    }
                                  >
                                    <Select mode="multiple" placeholder="选择隐式溶剂（可多选）" showSearch>
                                      <Select.OptGroup label="📌 水系电解液 (ε>50)">
                                        <Select.Option value="Water">水 (Water) ε=78.4</Select.Option>
                                      </Select.OptGroup>
                                      <Select.OptGroup label="📌 高介电常数碳酸酯 (ε=40-90)">
                                        <Select.Option value="DiMethylSulfoxide">DMSO ε=46.8 (离子液体参考)</Select.Option>
                                        <Select.Option value="1,2-EthaneDiol">乙二醇 ε=40.2</Select.Option>
                                      </Select.OptGroup>
                                      <Select.OptGroup label="📌 中等介电常数 (ε=15-40)">
                                        <Select.Option value="Acetonitrile">乙腈 ε=35.7</Select.Option>
                                        <Select.Option value="Methanol">甲醇 ε=32.6</Select.Option>
                                        <Select.Option value="Ethanol">乙醇 ε=24.9</Select.Option>
                                        <Select.Option value="Acetone">丙酮 ε=20.5 (高浓电解液)</Select.Option>
                                        <Select.Option value="1-Propanol">正丙醇 ε=20.5</Select.Option>
                                      </Select.OptGroup>
                                      <Select.OptGroup label="📌 低介电常数 (ε<15) - DMC/EMC/DEC体系">
                                        <Select.Option value="DiChloroEthane">二氯乙烷 ε=10.1</Select.Option>
                                        <Select.Option value="Dichloromethane">二氯甲烷 ε=8.9</Select.Option>
                                        <Select.Option value="TetraHydroFuran">四氢呋喃 (THF) ε=7.4</Select.Option>
                                        <Select.Option value="Chloroform">氯仿 ε=4.7 (线性碳酸酯参考)</Select.Option>
                                        <Select.Option value="DiethylEther">乙醚 ε=4.2</Select.Option>
                                        <Select.Option value="CarbonTetraChloride">四氯化碳 ε=2.2</Select.Option>
                                        <Select.Option value="Toluene">甲苯 ε=2.4</Select.Option>
                                        <Select.Option value="Benzene">苯 ε=2.3</Select.Option>
                                      </Select.OptGroup>
                                    </Select>
                                  </Form.Item>
                                );
                              }}
                            </Form.Item>
                          </Col>
                        </Row>

                        {/* 自定义溶剂参数输入 */}
                        <Form.Item noStyle shouldUpdate={(prevValues, currentValues) =>
                          prevValues.qc_solvent_models !== currentValues.qc_solvent_models
                        }>
                          {({ getFieldValue }) => {
                            const solventModels = getFieldValue('qc_solvent_models') || [];
                            if (!solventModels.includes('custom')) return null;
                            return (
                              <Card size="small" style={{ marginBottom: 12, background: isDark ? 'rgba(250, 173, 20, 0.15)' : '#fffbe6', borderColor: token.colorWarning }}>
                                <Text strong style={{ display: 'block', marginBottom: 8 }}>
                                  🔧 自定义溶剂参数（SMD模型）
                                </Text>
                                <Row gutter={[8, 8]}>
                                  <Col span={6}>
                                    <Form.Item name={['custom_solvent', 'eps']} label="介电常数 ε" style={{ marginBottom: 4 }} rules={[{ required: true, message: '请输入介电常数' }]}>
                                      <InputNumber style={{ width: '100%' }} placeholder="如: 89.6 (EC)" step={0.1} min={1} />
                                    </Form.Item>
                                  </Col>
                                  <Col span={6}>
                                    <Form.Item name={['custom_solvent', 'eps_inf']} label="光学介电常数 n²" style={{ marginBottom: 4 }}>
                                      <InputNumber style={{ width: '100%' }} placeholder="如: 2.2" step={0.01} min={1} />
                                    </Form.Item>
                                  </Col>
                                  <Col span={6}>
                                    <Form.Item name={['custom_solvent', 'hbond_acidity']} label="氢键酸度 α" style={{ marginBottom: 4 }}>
                                      <InputNumber style={{ width: '100%' }} placeholder="0.00-1.00" min={0} max={1} step={0.01} />
                                    </Form.Item>
                                  </Col>
                                  <Col span={6}>
                                    <Form.Item name={['custom_solvent', 'hbond_basicity']} label="氢键碱度 β" style={{ marginBottom: 4 }}>
                                      <InputNumber style={{ width: '100%' }} placeholder="0.00-1.00" min={0} max={1} step={0.01} />
                                    </Form.Item>
                                  </Col>
                                  <Col span={8}>
                                    <Form.Item name={['custom_solvent', 'surface_tension']} label="表面张力 γ" style={{ marginBottom: 4 }}>
                                      <InputNumber style={{ width: '100%' }} placeholder="cal/mol·Å²" step={0.1} />
                                    </Form.Item>
                                  </Col>
                                  <Col span={8}>
                                    <Form.Item name={['custom_solvent', 'carbon_aromaticity']} label="芳香碳比例 φ" style={{ marginBottom: 4 }}>
                                      <InputNumber style={{ width: '100%' }} placeholder="0.00-1.00" min={0} max={1} step={0.01} />
                                    </Form.Item>
                                  </Col>
                                  <Col span={8}>
                                    <Form.Item name={['custom_solvent', 'halogenicity']} label="卤素比例 ψ" style={{ marginBottom: 4 }}>
                                      <InputNumber style={{ width: '100%' }} placeholder="0.00-1.00" min={0} max={1} step={0.01} />
                                    </Form.Item>
                                  </Col>
                                </Row>
                                <Alert
                                  type="info"
                                  showIcon
                                  style={{ marginTop: 8 }}
                                  message={
                                    <Text style={{ fontSize: 11 }}>
                                      常用电解液介电常数参考：EC(ε≈89.6), PC(ε≈64.9), DMC(ε≈3.1), EMC(ε≈2.9), DEC(ε≈2.8)
                                    </Text>
                                  }
                                />
                              </Card>
                            );
                          }}
                        </Form.Item>

                        {/* 溶剂选择提示 */}
                        <Alert
                          type="info"
                          showIcon
                          style={{ marginBottom: 12 }}
                          message={
                            <Text style={{ fontSize: 12 }}>
                              <strong>隐式溶剂选择提示：</strong>选择介电常数(ε)接近您电解液的溶剂，或使用"自定义溶剂参数"输入精确值。
                              例如：EC体系选Water(ε≈78)或自定义(ε=89.6)，DMC/EMC体系选Chloroform(ε≈4.7)或自定义。
                            </Text>
                          }
                        />

                        {/* 显示将要计算的分子列表和任务数量 */}
                        <Form.Item noStyle shouldUpdate={(prevValues, currentValues) =>
                          prevValues.qc_functionals !== currentValues.qc_functionals ||
                          prevValues.qc_basis_sets !== currentValues.qc_basis_sets ||
                          prevValues.qc_solvent_models !== currentValues.qc_solvent_models ||
                          prevValues.qc_solvents !== currentValues.qc_solvents
                        }>
                          {({ getFieldValue }) => {
                            const functionals = getFieldValue('qc_functionals') || ['B3LYP'];
                            const basisSets = getFieldValue('qc_basis_sets') || ['6-31++g(d,p)'];
                            const solventModels = getFieldValue('qc_solvent_models') || ['pcm'];
                            const solvents = getFieldValue('qc_solvents') || ['Water'];

                            // 计算溶剂组合数
                            let solventCombinations = 0;
                            if (solventModels.includes('gas')) {
                              solventCombinations += 1;
                            }
                            if (solventModels.includes('custom')) {
                              solventCombinations += 1; // 自定义溶剂只有一个组合
                            }
                            const standardModels = solventModels.filter((m: string) => m !== 'gas' && m !== 'custom');
                            solventCombinations += standardModels.length * solvents.length;

                            const totalJobs = moleculesToCalc.length * functionals.length * basisSets.length * solventCombinations;

                            return moleculesToCalc.length > 0 ? (
                              <Alert
                                type="info"
                                showIcon
                                style={{ marginTop: 8 }}
                                message={
                                  <div>
                                    <strong>将创建 {totalJobs} 个 QC 任务</strong>
                                    <Text type="secondary" style={{ fontSize: 12, marginLeft: 8 }}>
                                      ({moleculesToCalc.length} 分子 × {functionals.length} 泛函 × {basisSets.length} 基组 × {solventCombinations} 溶剂组合)
                                    </Text>
                                  </div>
                                }
                                description={
                                  <div style={{ marginTop: 8 }}>
                                    <div style={{ marginBottom: 8 }}>
                                      <Text strong style={{ fontSize: 12 }}>分子列表：</Text>
                                    </div>
                                    {moleculesToCalc.map((mol, index) => (
                                      <div key={index} style={{
                                        display: 'inline-block',
                                        marginRight: 8,
                                        marginBottom: 4,
                                        padding: '2px 8px',
                                        background: mol.type === 'solvent' ? token.colorSuccessBg :
                                          mol.type === 'cation' ? token.colorErrorBg : token.colorInfoBg,
                                        borderRadius: 4,
                                        border: `1px solid ${mol.type === 'solvent' ? token.colorSuccessBorder :
                                          mol.type === 'cation' ? token.colorErrorBorder : token.colorInfoBorder}`
                                      }}>
                                        <Text style={{ fontSize: 12 }}>
                                          {mol.name}
                                          <Text type="secondary" style={{ fontSize: 11, marginLeft: 4 }}>
                                            ({mol.type === 'solvent' ? '溶剂' :
                                              mol.type === 'cation' ? '阳离子' : '阴离子'})
                                          </Text>
                                        </Text>
                                      </div>
                                    ))}
                                  </div>
                                }
                              />
                            ) : null;
                          }}
                        </Form.Item>
                      </Card>
                    </>
                  );
                }}
              </Form.Item>

            </Card>

            <Card title="计算资源" style={{ marginBottom: 24 }}>
              <Row gutter={16}>
                <Col span={12}>
                  <Form.Item
                    name="slurm_partition"
                    label={
                      <Space>
                        计算分区
                        <Tooltip title="选择计算资源所在的分区">
                          <QuestionCircleOutlined style={{ color: '#8c8c8c' }} />
                        </Tooltip>
                      </Space>
                    }
                    rules={[{ required: true, message: '请选择队列' }]}
                    extra="💡 不同分区有不同的资源限制"
                  >
                    <Select
                      placeholder={partitions.length > 0 ? "选择队列" : "暂无可用队列"}
                      disabled={partitions.length === 0}
                    >
                      {partitions.map(p => (
                        <Select.Option
                          key={p.name}
                          value={p.name}
                          disabled={p.state !== 'up'}
                        >
                          <span style={{ color: p.state === 'up' ? 'inherit' : '#999' }}>
                            {p.name} {p.state === 'up'
                              ? `✅ ${p.available_cpus}/${p.total_cpus} CPUs`
                              : '(不可用)'}
                          </span>
                        </Select.Option>
                      ))}
                    </Select>
                  </Form.Item>
                  {partitions.length === 0 && (
                    <Alert
                      message="暂无可用队列"
                      description="请联系管理员分配队列权限"
                      type="warning"
                      showIcon
                      style={{ marginBottom: 16 }}
                    />
                  )}
                </Col>

                <Col span={12}>
                  <Form.Item
                    name="slurm_cpus"
                    label={
                      <Space>
                        CPU核心数
                        <Tooltip title={
                          <div>
                            <p><strong>选择建议:</strong></p>
                            <p>• 16核: 通用任务（默认）</p>
                            <p>• 32核: 大体系/高精度</p>
                            <p>• 64核: 超大体系</p>
                          </div>
                        }>
                          <QuestionCircleOutlined style={{ color: '#8c8c8c' }} />
                        </Tooltip>
                      </Space>
                    }
                    initialValue={16}
                    extra="💡 默认 16 核适合大多数任务"
                  >
                    <Select>
                      <Select.Option value={16}>
                        <Space>
                          <span>16 核</span>
                          <Tag color="blue" style={{ margin: 0 }}>推荐</Tag>
                        </Space>
                      </Select.Option>
                      <Select.Option value={32}>
                        <Space>
                          <span>32 核</span>
                          <Tag color="orange" style={{ margin: 0 }}>高级</Tag>
                        </Space>
                      </Select.Option>
                      <Select.Option value={64}>
                        <Space>
                          <span>64 核</span>
                          <Tag color="red" style={{ margin: 0 }}>高级</Tag>
                        </Space>
                      </Select.Option>
                    </Select>
                  </Form.Item>
                </Col>
              </Row>

              <Form.Item
                name="slurm_time"
                initialValue={7200}
                hidden
              >
                <InputNumber />
              </Form.Item>


            </Card>
            <Card style={{ boxShadow: '0 -2px 10px rgba(0,0,0,0.05)' }} bodyStyle={{ padding: '12px 24px' }}>
              <Space style={{ width: '100%', justifyContent: 'flex-end' }}>
                <Button onClick={() => navigate('/workspace/electrolytes')} size="large">
                  取消
                </Button>
                <Button
                  icon={<SaveOutlined />}
                  onClick={handleSave}
                  loading={submitting}
                  size="large"
                  disabled={quota ? !quota.can_create : false}
                >
                  保存草稿
                </Button>
                <Button
                  type="primary"
                  icon={<RocketOutlined />}
                  onClick={handleSubmit}
                  loading={submitting}
                  size="large"
                  disabled={quota ? !quota.can_create : false}
                >
                  提交任务
                </Button>
              </Space>
            </Card>
          </Col>
        </Row>
      </div>
    </Form >
  );
}
