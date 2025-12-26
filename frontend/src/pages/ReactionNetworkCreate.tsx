/**
 * Create Reaction Network Job Page
 * 创建反应网络任务 - 统一使用离子+溶剂选择器
 */

import React, { useState, useEffect } from 'react';
import {
    Card,
    Form,
    Input,
    InputNumber,
    Select,
    Button,
    Row,
    Col,
    message,
    Space,
    Tag,
    Tooltip,
    Alert,
    Divider,
    Result
} from 'antd';

import {
    ExperimentOutlined,
    FireOutlined,
    ThunderboltOutlined,
    SettingOutlined,
    RocketOutlined,
    InfoCircleOutlined,
    DeleteOutlined
} from '@ant-design/icons';
import { useNavigate } from 'react-router-dom';
import { createReactionNetworkJob, submitReactionNetworkJob } from '../api/reactionNetwork';
import { getAvailableIons } from '../api/electrolytes';
import AnionSelectorWithGeneration from '../components/AnionSelectorWithGeneration';

const { TextArea } = Input;
const { Option } = Select;

// Define IonInfo locally since it's not exported from api
interface IonInfo {
    name: string;
    charge: number;
    smiles: string;
}

// 溶剂列表（复用配方页面的常用溶剂）
const SOLVENTS = [
    { name: 'EC', fullName: '碳酸乙烯酯', smiles: 'C1COC(=O)O1' },
    { name: 'PC', fullName: '碳酸丙烯酯', smiles: 'CC1COC(=O)O1' },
    { name: 'DMC', fullName: '碳酸二甲酯', smiles: 'COC(=O)OC' },
    { name: 'DEC', fullName: '碳酸二乙酯', smiles: 'CCOC(=O)OCC' },
    { name: 'EMC', fullName: '碳酸甲乙酯', smiles: 'CCOC(=O)OC' },
    { name: 'FEC', fullName: '氟代碳酸乙烯酯', smiles: 'C1C(OC(=O)O1)F' },
    { name: 'DME', fullName: '乙二醇二甲醚', smiles: 'COCCOC' },
    { name: 'DOL', fullName: '1,3-二氧戊环', smiles: 'C1COCO1' },
    { name: 'Sulfolane', fullName: '环丁砜', smiles: 'C1CCS(=O)(=O)C1' },
];

interface SelectedIon {
    name: string;
    charge: number;
    smiles: string;
}

interface SelectedSolvent {
    name: string;
    smiles: string;
}

const ReactionNetworkCreate: React.FC = () => {
    const navigate = useNavigate();
    const [form] = Form.useForm();
    const [submitting, setSubmitting] = useState(false);
    const [currentStep, setCurrentStep] = useState(0);
    const [createdJobId, setCreatedJobId] = useState<number | null>(null);

    // 离子和溶剂状态
    const [availableCations, setAvailableCations] = useState<IonInfo[]>([]);
    const [availableAnions, setAvailableAnions] = useState<IonInfo[]>([]);
    const [selectedCations, setSelectedCations] = useState<SelectedIon[]>([]);
    const [selectedAnions, setSelectedAnions] = useState<SelectedIon[]>([]);
    const [selectedSolvents, setSelectedSolvents] = useState<SelectedSolvent[]>([]);
    const [loading, setLoading] = useState(false);
    const [electrodeType, setElectrodeType] = useState<'anode' | 'cathode'>('anode');

    useEffect(() => {
        loadAvailableIons();
    }, []);

    const loadAvailableIons = async () => {
        setLoading(true);
        try {
            const data = await getAvailableIons();
            setAvailableCations((data.cations as any[]).map(c => ({ ...c, smiles: c.smiles || '' })));
            setAvailableAnions((data.anions as any[]).map(a => ({ ...a, smiles: a.smiles || '' })));
        } catch (error: any) {
            message.error('加载可用离子列表失败');
        } finally {
            setLoading(false);
        }
    };

    const addCation = (ionName: string) => {
        const ion = availableCations.find(i => i.name === ionName);
        if (ion && !selectedCations.find(c => c.name === ionName)) {
            setSelectedCations([...selectedCations, {
                name: ion.name,
                charge: ion.charge,
                smiles: ion.smiles
            }]);
        }
    };

    const addAnion = (ionName: string) => {
        const ion = availableAnions.find(i => i.name === ionName);
        if (ion && !selectedAnions.find(a => a.name === ionName)) {
            setSelectedAnions([...selectedAnions, {
                name: ion.name,
                charge: ion.charge,
                smiles: ion.smiles
            }]);
        }
    };

    const removeCation = (ionName: string) => {
        setSelectedCations(selectedCations.filter(c => c.name !== ionName));
    };

    const removeAnion = (ionName: string) => {
        setSelectedAnions(selectedAnions.filter(a => a.name !== ionName));
    };

    const addSolvent = (solventName: string) => {
        const solvent = SOLVENTS.find(s => s.name === solventName);
        if (solvent && !selectedSolvents.find(s => s.name === solventName)) {
            setSelectedSolvents([...selectedSolvents, {
                name: solvent.name,
                smiles: solvent.smiles
            }]);
        }
    };

    const removeSolvent = (solventName: string) => {
        setSelectedSolvents(selectedSolvents.filter(s => s.name !== solventName));
    };

    const onFinish = async (values: any) => {
        setSubmitting(true);
        try {
            // 验证至少有一个组分
            if (selectedCations.length === 0 && selectedAnions.length === 0 && selectedSolvents.length === 0) {
                message.error('请至少选择一个阳离子、阴离子或溶剂');
                setSubmitting(false);
                return;
            }

            // 将选中的离子和溶剂转换为SMILES列表
            const smilesList: string[] = [];
            selectedCations.forEach(ion => smilesList.push(ion.smiles));
            selectedAnions.forEach(ion => smilesList.push(ion.smiles));
            selectedSolvents.forEach(solvent => smilesList.push(solvent.smiles));

            const jobData = {
                job_name: values.job_name,
                description: values.description,
                initial_smiles: smilesList,
                temperature: values.temperature,
                electrode_type: values.electrode_type,
                anode_material: values.anode_material,
                cathode_material: values.cathode_material,
                voltage: values.voltage,
                max_generations: values.max_generations,
                max_species: values.max_species,
                energy_cutoff: values.energy_cutoff,
                slurm_partition: values.slurm_partition,
                slurm_cpus: values.slurm_cpus,
                slurm_time: values.slurm_time
            };

            const job = await createReactionNetworkJob(jobData);
            setCreatedJobId(job.id);
            message.success('任务创建成功！');

            if (values.submit_immediately) {
                await submitReactionNetworkJob(job.id);
                message.success('任务已提交到计算队列');
                navigate(`/workspace/liquid-electrolyte/reaction-network/${job.id}`);
            } else {
                setCurrentStep(1);
            }
        } catch (error: any) {
            message.error(error.response?.data?.detail || '任务创建失败');
        } finally {
            setSubmitting(false);
        }
    };

    const handleSubmitNow = async () => {
        if (!createdJobId) return;
        try {
            await submitReactionNetworkJob(createdJobId);
            message.success('任务已提交！');
            navigate(`/workspace/liquid-electrolyte/reaction-network/${createdJobId}`);
        } catch (error) {
            message.error('提交失败');
        }
    };

    // 预设模板
    const applyTemplate = (templateName: string) => {
        if (templateName === 'ec_lipf6') {
            setSelectedCations([{ name: 'Li+', charge: 1, smiles: '[Li+]' }]);
            setSelectedAnions([{ name: 'PF6-', charge: -1, smiles: 'F[P-](F)(F)(F)(F)F' }]);
            setSelectedSolvents([{ name: 'EC', smiles: 'C1COC(=O)O1' }]);
            form.setFieldsValue({
                job_name: 'EC + LiPF6 电解液反应网络',
                electrode_type: 'anode',
                voltage: 0.1
            });
            message.success('已应用 EC + LiPF6 模板');
        } else if (templateName === 'dmc_litfsi') {
            setSelectedCations([{ name: 'Li+', charge: 1, smiles: '[Li+]' }]);
            setSelectedAnions([{ name: 'TFSI-', charge: -1, smiles: 'O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F' }]);
            setSelectedSolvents([{ name: 'DMC', smiles: 'COC(=O)OC' }]);
            form.setFieldsValue({
                job_name: 'DMC + LiTFSI 电解液反应网络',
                electrode_type: 'anode',
                voltage: 0.1
            });
            message.success('已应用 DMC + LiTFSI 模板');
        }
    };

    if (currentStep === 1 && createdJobId) {
        return (
            <div style={{ padding: '24px', maxWidth: '800px', margin: '0 auto' }}>
                <Card>
                    <Result
                        status="success"
                        title="任务创建成功！"
                        subTitle={`任务ID: ${createdJobId}`}
                        extra={[
                            <Button
                                type="primary"
                                key="submit"
                                icon={<RocketOutlined />}
                                onClick={handleSubmitNow}
                            >
                                立即提交到计算队列
                            </Button>,
                            <Button
                                key="view"
                                onClick={() => navigate(`/workspace/liquid-electrolyte/reaction-network/${createdJobId}`)}
                            >
                                查看任务详情
                            </Button>,
                            <Button key="list" onClick={() => navigate('/workspace/liquid-electrolyte/reaction-network')}>
                                返回任务列表
                            </Button>
                        ]}
                    />
                </Card>
            </div>
        );
    }

    return (
        <div style={{ padding: '24px', maxWidth: '1200px', margin: '0 auto' }}>
            <Card
                title={
                    <Space size="large">
                        <ExperimentOutlined style={{ fontSize: '24px', color: '#1890ff' }} />
                        <span style={{ fontSize: '20px' }}>创建反应网络生成任务</span>
                    </Space>
                }
                extra={
                    <Space>
                        <Tooltip title="使用预设模板快速开始">
                            <Select
                                placeholder="选择模板"
                                style={{ width: 200 }}
                                onChange={applyTemplate}
                                allowClear
                            >
                                <Option value="ec_lipf6">EC + LiPF6</Option>
                                <Option value="dmc_litfsi">DMC + LiTFSI</Option>
                            </Select>
                        </Tooltip>
                    </Space>
                }
            >
                <Alert
                    message="什么是反应网络生成？"
                    description="基于初始分子（阳离子+阴离子+溶剂），通过智能算符系统自动发现可能的化学反应，使用XTB半经验方法计算能量，构建完整的反应网络。适用于电池电解液SEI/CEI形成、催化反应筛选、降解机制研究等场景。"
                    type="info"
                    showIcon
                    icon={<InfoCircleOutlined />}
                    style={{ marginBottom: '24px' }}
                />

                <Form
                    form={form}
                    layout="vertical"
                    onFinish={onFinish}
                    initialValues={{
                        temperature: 300.0,
                        electrode_type: 'anode',
                        anode_material: 'GRAPHITE',
                        cathode_material: 'NMC',
                        voltage: 0.1,
                        max_generations: 3,
                        max_species: 50,
                        energy_cutoff: 80.0,
                        slurm_partition: 'cpu',
                        slurm_cpus: 16,
                        slurm_time: 7200,
                        submit_immediately: true
                    }}

                >
                    {/* 基本信息 */}
                    <Card type="inner" title="📝 基本信息" style={{ marginBottom: '16px' }}>
                        <Row gutter={16}>
                            <Col span={16}>
                                <Form.Item
                                    name="job_name"
                                    label="任务名称"
                                    rules={[{ required: true, message: '请输入任务名称' }]}
                                >
                                    <Input
                                        placeholder="例如: EC电解液反应网络分析"
                                        prefix={<ExperimentOutlined />}
                                        size="large"
                                    />
                                </Form.Item>
                            </Col>
                            <Col span={8}>
                                <Form.Item name="submit_immediately" label="创建后" valuePropName="checked">
                                    <Select size="large">
                                        <Option value={true}>立即提交运行</Option>
                                        <Option value={false}>暂不提交</Option>
                                    </Select>
                                </Form.Item>
                            </Col>
                        </Row>

                        <Form.Item name="description" label="任务描述（可选）">
                            <TextArea rows={2} placeholder="描述这个反应网络任务的目的和预期..." />
                        </Form.Item>
                    </Card>

                    {/* 初始分子配置 - 离子和溶剂选择 */}
                    <Card type="inner" title="🧪 初始分子配置" style={{ marginBottom: '16px' }}>
                        <Alert
                            message="选择电解液组分"
                            description="请选择阳离子、阴离子和溶剂，系统将自动使用这些分子作为初始分子进行反应网络生成"
                            type="info"
                            showIcon
                            style={{ marginBottom: '16px' }}
                        />

                        <Row gutter={16}>
                            {/* 阳离子 */}
                            <Col span={8}>
                                <div>
                                    <div style={{ marginBottom: 8 }}>
                                        <Tag color="blue">阳离子</Tag>
                                    </div>
                                    <Select
                                        placeholder="选择阳离子"
                                        style={{ width: '100%', marginBottom: 12 }}
                                        onChange={addCation}
                                        value={undefined}
                                        showSearch
                                        filterOption={(input, option) =>
                                            String(option?.children ?? '').toLowerCase().includes(input.toLowerCase())
                                        }
                                    >
                                        {availableCations
                                            .filter(ion => !selectedCations.find(c => c.name === ion.name))
                                            .map(ion => (
                                                <Option key={ion.name} value={ion.name}>
                                                    {ion.name} (+{ion.charge})
                                                </Option>
                                            ))}
                                    </Select>
                                    <div>
                                        {selectedCations.map(ion => (
                                            <Tag
                                                key={ion.name}
                                                color="blue"
                                                closable
                                                onClose={() => removeCation(ion.name)}
                                                style={{ marginBottom: 8 }}
                                            >
                                                {ion.name}
                                            </Tag>
                                        ))}
                                        {selectedCations.length === 0 && (
                                            <Alert message="请选择阳离子" type="warning" showIcon style={{ padding: '4px 8px' }} />
                                        )}
                                    </div>
                                </div>
                            </Col>

                            {/* 阴离子 */}
                            <Col span={8}>
                                <div>
                                    <div style={{ marginBottom: 8 }}>
                                        <Tag color="red">阴离子</Tag>
                                    </div>
                                    <AnionSelectorWithGeneration
                                        availableAnions={availableAnions}
                                        selectedAnions={selectedAnions.map(a => ({ ...a, concentration: 1.0 }))}
                                        onAddAnion={addAnion}
                                        onRefresh={loadAvailableIons}
                                    />
                                    <div style={{ marginTop: 12 }}>
                                        {selectedAnions.map(ion => (
                                            <Tag
                                                key={ion.name}
                                                color="red"
                                                closable
                                                onClose={() => removeAnion(ion.name)}
                                                style={{ marginBottom: 8 }}
                                            >
                                                {ion.name}
                                            </Tag>
                                        ))}
                                        {selectedAnions.length === 0 && (
                                            <Alert message="请选择阴离子" type="warning" showIcon style={{ padding: '4px 8px' }} />
                                        )}
                                    </div>
                                </div>
                            </Col>

                            {/* 溶剂 */}
                            <Col span={8}>
                                <div>
                                    <div style={{ marginBottom: 8 }}>
                                        <Tag color="green">溶剂</Tag>
                                    </div>
                                    <Select
                                        placeholder="选择溶剂"
                                        style={{ width: '100%', marginBottom: 12 }}
                                        onChange={addSolvent}
                                        value={undefined}
                                        showSearch
                                        filterOption={(input, option) =>
                                            String(option?.children ?? '').toLowerCase().includes(input.toLowerCase())
                                        }
                                    >
                                        {SOLVENTS
                                            .filter(s => !selectedSolvents.find(sel => sel.name === s.name))
                                            .map(solvent => (
                                                <Option key={solvent.name} value={solvent.name}>
                                                    {solvent.name} - {solvent.fullName}
                                                </Option>
                                            ))}
                                    </Select>
                                    <div>
                                        {selectedSolvents.map(solvent => (
                                            <Tag
                                                key={solvent.name}
                                                color="green"
                                                closable
                                                onClose={() => removeSolvent(solvent.name)}
                                                style={{ marginBottom: 8 }}
                                            >
                                                {solvent.name}
                                            </Tag>
                                        ))}
                                        {selectedSolvents.length === 0 && (
                                            <Alert message="可选：添加溶剂分子" type="info" showIcon style={{ padding: '4px 8px' }} />
                                        )}
                                    </div>
                                </div>
                            </Col>
                        </Row>
                    </Card>

                    {/* 环境参数 */}
                    <Card type="inner" title={<><ThunderboltOutlined /> 环境参数</>} style={{ marginBottom: '16px' }}>
                        <Row gutter={16}>
                            <Col span={8}>
                                <Form.Item
                                    name="temperature"
                                    label="温度 (K)"
                                    rules={[{ required: true }]}
                                >
                                    <InputNumber
                                        min={0}
                                        max={1000}
                                        step={10}
                                        style={{ width: '100%' }}
                                        addonAfter="K"
                                    />
                                </Form.Item>
                            </Col>
                            <Col span={8}>
                                <Form.Item
                                    name="electrode_type"
                                    label="电极类型"
                                    rules={[{ required: true }]}
                                    tooltip="阳极研究SEI，阴极研究CEI"
                                >
                                    <Select onChange={(value) => setElectrodeType(value as 'anode' | 'cathode')}>
                                        <Option value="anode">阳极 (负极) - SEI形成</Option>
                                        <Option value="cathode">阴极 (正极) - CEI形成</Option>
                                    </Select>
                                </Form.Item>
                            </Col>
                            <Col span={8}>
                                <Form.Item
                                    name="voltage"
                                    label="电压 (V)"
                                    rules={[{ required: true }]}
                                >
                                    <InputNumber
                                        min={0}
                                        max={5}
                                        step={0.1}
                                        style={{ width: '100%' }}
                                        addonAfter="V"
                                    />
                                </Form.Item>
                            </Col>
                        </Row>

                        <Divider orientation="left" plain>
                            电极材料 (影响自动注入物种)
                        </Divider>

                        <Alert
                            message={`已选择${electrodeType === 'anode' ? '阳极(负极)' : '阴极(正极)'}，仅需配置对应材料`}
                            type="info"
                            showIcon
                            style={{ marginBottom: 16 }}
                        />

                        <Row gutter={16}>
                            {electrodeType === 'anode' ? (
                                <Col span={24}>
                                    <Form.Item
                                        name="anode_material"
                                        label="负极材料"
                                        tooltip="决定载流子类型(Li/Na/K)和SEI化学"
                                        extra="系统将自动注入相应的载流子离子和金属原子"
                                        rules={[{ required: true, message: '请选择负极材料' }]}
                                    >
                                        <Select size="large">
                                            <Option value="GRAPHITE">石墨 (Li+) - 成熟SEI</Option>
                                            <Option value="LI_METAL">锂金属 (Li+) - 金属SEI</Option>
                                            <Option value="SILICON">硅负极 (Li+) - Si-O-Li SEI</Option>
                                            <Option value="SIC">硅碳复合 (Li+)</Option>
                                            <Option value="LTO">钛酸锂 LTO (Li+) - 无SEI</Option>
                                            <Option value="NA_METAL">钠金属 (Na+)</Option>
                                            <Option value="HARD_CARBON">硬碳 (Na+)</Option>
                                            <Option value="SOFT_CARBON">软碳 (Na+)</Option>
                                            <Option value="K_METAL">钾金属 (K+)</Option>
                                            <Option value="K_GRAPHITE">钾石墨 KC8 (K+)</Option>
                                        </Select>
                                    </Form.Item>
                                </Col>
                            ) : (
                                <Col span={24}>
                                    <Form.Item
                                        name="cathode_material"
                                        label="正极材料"
                                        tooltip="决定氧释放行为和CEI化学"
                                        extra="系统将根据电压自动注入氧物种和自由基"
                                        rules={[{ required: true, message: '请选择正极材料' }]}
                                    >
                                        <Select size="large">
                                            <Option value="NMC">NMC三元 (通用) - 中等氧释放</Option>
                                            <Option value="NMC811">NMC811 高镍 - 高氧释放</Option>
                                            <Option value="NMC622">NMC622 - 适度氧释放</Option>
                                            <Option value="LCO">钴酸锂 LCO - 氧释放</Option>
                                            <Option value="NCA">NCA镍钴铝 - 高氧释放</Option>
                                            <Option value="LFP">磷酸铁锂 LFP - 无氧释放</Option>
                                            <Option value="LMO">锰酸锂 LMO</Option>
                                            <Option value="LNMO">高电压尖晶石 LNMO - 强氧释放</Option>
                                            <Option value="LRLO">富锂层状氧化物 - 氧损失</Option>
                                        </Select>
                                    </Form.Item>
                                </Col>
                            )}
                        </Row>
                    </Card>

                    {/* 网络生成参数 */}
                    <Card type="inner" title={<><SettingOutlined /> 网络生成参数</>} style={{ marginBottom: '16px' }}>
                        <Row gutter={16}>
                            <Col span={8}>
                                <Form.Item
                                    name="max_generations"
                                    label="最大代数"
                                    tooltip="控制反应网络的深度"
                                    rules={[{ required: true }]}
                                >
                                    <InputNumber min={1} max={10} style={{ width: '100%' }} />
                                </Form.Item>
                            </Col>
                            <Col span={8}>
                                <Form.Item
                                    name="max_species"
                                    label="最大分子数"
                                    tooltip="限制网络中的分子总数"
                                    rules={[{ required: true }]}
                                >
                                    <InputNumber min={10} max={500} step={10} style={{ width: '100%' }} />
                                </Form.Item>
                            </Col>
                            <Col span={8}>
                                <Form.Item
                                    name="energy_cutoff"
                                    label="能量截断 (kcal/mol)"
                                    tooltip="过滤高能量反应"
                                    rules={[{ required: true }]}
                                >
                                    <InputNumber min={0} max={200} step={10} style={{ width: '100%' }} addonAfter="kcal/mol" />
                                </Form.Item>
                            </Col>
                        </Row>
                    </Card>

                    {/* 计算资源 */}
                    <Card type="inner" title={<><FireOutlined /> 计算资源</>} style={{ marginBottom: '16px' }}>
                        <Row gutter={16}>
                            <Col span={8}>
                                <Form.Item name="slurm_partition" label="计算分区">
                                    <Select>
                                        <Option value="cpu">CPU分区</Option>
                                        <Option value="gpu">GPU分区</Option>
                                    </Select>
                                </Form.Item>
                            </Col>
                            <Col span={8}>
                                <Form.Item name="slurm_cpus" label="CPU核数">
                                    <InputNumber min={1} max={128} style={{ width: '100%' }} />
                                </Form.Item>
                            </Col>
                            <Col span={8}>
                                <Form.Item name="slurm_time" label="最大运行时间 (分钟)">
                                    <InputNumber min={60} max={43200} step={60} style={{ width: '100%' }} addonAfter="分钟" />
                                </Form.Item>
                            </Col>
                        </Row>
                    </Card>

                    {/* 提交按钮 */}
                    <Form.Item>
                        <Space size="large" style={{ width: '100%', justifyContent: 'center' }}>
                            <Button onClick={() => navigate('/workspace/liquid-electrolyte/reaction-network')} size="large">
                                取消
                            </Button>
                            <Button
                                type="primary"
                                htmlType="submit"
                                size="large"
                                icon={<RocketOutlined />}
                                loading={submitting}
                            >
                                创建任务
                            </Button>
                        </Space>
                    </Form.Item>
                </Form>
            </Card>
        </div>
    );
};

export default ReactionNetworkCreate;
