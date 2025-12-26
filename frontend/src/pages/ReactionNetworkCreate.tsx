/**
 * Create Reaction Network Job Page
 * 创建反应网络任务 - 现代化表单设计
 */

import React, { useState } from 'react';
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
    Steps,
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
    InfoCircleOutlined
} from '@ant-design/icons';
import { useNavigate } from 'react-router-dom';
import { createReactionNetworkJob, submitReactionNetworkJob } from '../api/reactionNetwork';

const { TextArea } = Input;
const { Option } = Select;

const ReactionNetworkCreate: React.FC = () => {
    const navigate = useNavigate();
    const [form] = Form.useForm();
    const [submitting, setSubmitting] = useState(false);
    const [currentStep, setCurrentStep] = useState(0);
    const [createdJobId, setCreatedJobId] = useState<number | null>(null);

    const onFinish = async (values: any) => {
        setSubmitting(true);
        try {
            // 解析SMILES（每行一个）
            const smilesText = values.initial_smiles_text || '';
            const smilesList = smilesText
                .split('\n')
                .map((s: string) => s.trim())
                .filter((s: string) => s.length > 0);

            if (smilesList.length === 0) {
                message.error('请至少输入一个SMILES');
                setSubmitting(false);
                return;
            }

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
                navigate(`/reaction-network/${job.id}`);
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
            navigate(`/reaction-network/${createdJobId}`);
        } catch (error) {
            message.error('提交失败');
        }
    };

    // 预设模板
    const templates = {
        ec_lipf6: {
            name: 'EC + LiPF6 电解液',
            smiles: 'C1COC(=O)O1\n[Li+]\nF[P-](F)(F)(F)(F)F',
            temperature: 300,
            electrode_type: 'anode',
            voltage: 0.1,
            max_generations: 3
        },
        dmc_litfsi: {
            name: 'DMC + LiTFSI 电解液',
            smiles: 'COC(=O)OC\n[Li+]\nO=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F',
            temperature: 300,
            electrode_type: 'anode',
            voltage: 0.1,
            max_generations: 3
        }
    };

    const applyTemplate = (templateKey: keyof typeof templates) => {
        const template = templates[templateKey];
        form.setFieldsValue({
            job_name: template.name,
            initial_smiles_text: template.smiles,
            temperature: template.temperature,
            electrode_type: template.electrode_type,
            voltage: template.voltage,
            max_generations: template.max_generations
        });
        message.success('模板已应用');
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
                                onClick={() => navigate(`/reaction-network/${createdJobId}`)}
                            >
                                查看任务详情
                            </Button>,
                            <Button key="list" onClick={() => navigate('/reaction-network')}>
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
                    description="基于初始分子，通过智能算符系统自动发现可能的化学反应，使用XTB半经验方法计算能量，构建完整的反应网络。适用于电池电解液SEI形成、催化反应筛选、降解机制研究等场景。"
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

                    {/* 初始分子 */}
                    <Card type="inner" title="🧪 初始分子 (SMILES)" style={{ marginBottom: '16px' }}>
                        <Form.Item
                            name="initial_smiles_text"
                            label={
                                <span>
                                    SMILES表达式
                                    <Tooltip title="每行输入一个SMILES，支持离子（如[Li+]）、中性分子等">
                                        <InfoCircleOutlined style={{ marginLeft: 8, color: '#1890ff' }} />
                                    </Tooltip>
                                </span>
                            }
                            rules={[{ required: true, message: '请输入至少一个SMILES' }]}
                            extra={
                                <Space>
                                    <Tag color="blue">示例: C1COC(=O)O1 (EC)</Tag>
                                    <Tag color="green">示例: [Li+] (锂离子)</Tag>
                                    <Tag color="orange">示例: F[P-](F)(F)(F)(F)F (PF6-)</Tag>
                                </Space>
                            }
                        >
                            <TextArea
                                rows={8}
                                placeholder={'每行一个SMILES，例如:\nC1COC(=O)O1\n[Li+]\nF[P-](F)(F)(F)(F)F'}
                                style={{ fontFamily: 'monospace' }}
                            />
                        </Form.Item>
                    </Card>

                    {/* 环境参数 */}
                    <Card
                        type="inner"
                        title={
                            <Space>
                                <FireOutlined />
                                <span>环境参数</span>
                            </Space>
                        }
                        style={{ marginBottom: '16px' }}
                    >
                        <Row gutter={16}>
                            <Col span={6}>
                                <Form.Item
                                    name="temperature"
                                    label="温度 (K)"
                                    tooltip="影响反应的热激活过程"
                                >
                                    <InputNumber
                                        min={0}
                                        max={1000}
                                        style={{ width: '100%' }}
                                        addonAfter="K"
                                    />
                                </Form.Item>
                            </Col>
                            <Col span={6}>
                                <Form.Item
                                    name="electrode_type"
                                    label="电极类型"
                                    tooltip="影响驱动力的方向和物种注入"
                                >
                                    <Select>
                                        <Option value="anode">阳极 (负极)</Option>
                                        <Option value="cathode">阴极 (正极)</Option>
                                    </Select>
                                </Form.Item>
                            </Col>
                            <Col span={6}>
                                <Form.Item
                                    name="voltage"
                                    label="电压 (V)"
                                    tooltip="电极电势，影响氧化/还原反应"
                                >
                                    <InputNumber
                                        min={-10}
                                        max={10}
                                        step={0.1}
                                        style={{ width: '100%' }}
                                        addonAfter="V"
                                    />
                                </Form.Item>
                            </Col>
                        </Row>

                        <Divider orientation="left" plain>电极材料 (影响自动注入物种)</Divider>

                        <Row gutter={16}>
                            <Col span={12}>
                                <Form.Item
                                    name="anode_material"
                                    label="负极材料"
                                    tooltip="决定载流子类型(Li/Na/K)和SEI化学"
                                    extra="自动注入相应的离子和原子"
                                >
                                    <Select>
                                        <Option value="GRAPHITE">石墨 (Li+)</Option>
                                        <Option value="LI_METAL">锂金属 (Li+)</Option>
                                        <Option value="SILICON">硅负极 (Li+)</Option>
                                        <Option value="SIC">硅碳复合 (Li+)</Option>
                                        <Option value="LTO">钛酸锂 LTO (Li+)</Option>
                                        <Option value="NA_METAL">钠金属 (Na+)</Option>
                                        <Option value="HARD_CARBON">硬碳 (Na+)</Option>
                                        <Option value="SOFT_CARBON">软碳 (Na+)</Option>
                                        <Option value="K_METAL">钾金属 (K+)</Option>
                                        <Option value="K_GRAPHITE">钾石墨 KC8 (K+)</Option>
                                    </Select>
                                </Form.Item>
                            </Col>
                            <Col span={12}>
                                <Form.Item
                                    name="cathode_material"
                                    label="正极材料"
                                    tooltip="决定氧释放行为和CEI化学"
                                    extra="高电压时自动注入氧物种"
                                >
                                    <Select>
                                        <Option value="NMC">NMC三元 (通用)</Option>
                                        <Option value="NMC811">NMC811 高镍</Option>
                                        <Option value="NMC622">NMC622</Option>
                                        <Option value="LCO">钴酸锂 LCO</Option>
                                        <Option value="NCA">NCA镍钴铝</Option>
                                        <Option value="LFP">磷酸铁锂 LFP (稳定)</Option>
                                        <Option value="LMO">锰酸锂 LMO</Option>
                                        <Option value="LNMO">高电压尖晶石 LNMO</Option>
                                        <Option value="LRLO">富锂层状氧化物</Option>
                                    </Select>
                                </Form.Item>
                            </Col>
                        </Row>
                    </Card>


                    {/* 网络生成参数 */}
                    <Card
                        type="inner"
                        title={
                            <Space>
                                <ThunderboltOutlined />
                                <span>网络生成参数</span>
                            </Space>
                        }
                        style={{ marginBottom: '16px' }}
                    >
                        <Row gutter={16}>
                            <Col span={8}>
                                <Form.Item
                                    name="max_generations"
                                    label="最大代数"
                                    tooltip="反应迭代的最大轮数，建议1-5代"
                                    extra="代数越大，网络越复杂"
                                >
                                    <InputNumber min={1} max={10} style={{ width: '100%' }} />
                                </Form.Item>
                            </Col>
                            <Col span={8}>
                                <Form.Item
                                    name="max_species"
                                    label="最大分子数"
                                    tooltip="限制网络规模，避免爆炸性增长"
                                >
                                    <InputNumber min={1} max={200} style={{ width: '100%' }} />
                                </Form.Item>
                            </Col>
                            <Col span={8}>
                                <Form.Item
                                    name="energy_cutoff"
                                    label="能量截断 (kcal/mol)"
                                    tooltip="排除高能反应，只保留能量低于此值的反应"
                                >
                                    <InputNumber
                                        min={0}
                                        max={200}
                                        style={{ width: '100%' }}
                                        addonAfter="kcal/mol"
                                    />
                                </Form.Item>
                            </Col>
                        </Row>
                    </Card>

                    {/* Slurm资源配置 */}
                    <Card
                        type="inner"
                        title={
                            <Space>
                                <SettingOutlined />
                                <span>计算资源配置</span>
                            </Space>
                        }
                        style={{ marginBottom: '16px' }}
                    >
                        <Row gutter={16}>
                            <Col span={8}>
                                <Form.Item name="slurm_partition" label="Slurm队列">
                                    <Select>
                                        <Option value="cpu">CPU队列</Option>
                                        <Option value="gpu">GPU队列</Option>
                                        <Option value="fat">大内存队列</Option>
                                    </Select>
                                </Form.Item>
                            </Col>
                            <Col span={8}>
                                <Form.Item name="slurm_cpus" label="CPU核心数">
                                    <InputNumber min={1} max={128} style={{ width: '100%' }} />
                                </Form.Item>
                            </Col>
                            <Col span={8}>
                                <Form.Item
                                    name="slurm_time"
                                    label="最大运行时间"
                                    extra="单位: 分钟"
                                >
                                    <InputNumber min={10} max={43200} style={{ width: '100%' }} />
                                </Form.Item>
                            </Col>
                        </Row>
                    </Card>

                    <Divider />

                    {/* 提交按钮 */}
                    <Form.Item>
                        <Space size="large" style={{ width: '100%', justifyContent: 'center' }}>
                            <Button onClick={() => navigate('/reaction-network')} size="large">
                                取消
                            </Button>
                            <Button
                                type="primary"
                                htmlType="submit"
                                loading={submitting}
                                icon={<RocketOutlined />}
                                size="large"
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
