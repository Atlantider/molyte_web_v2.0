/**
 * Reaction Network Detail Page
 * 任务详情页 - 现代化设计
 */

import React, { useState, useEffect } from 'react';
import {
    Card,
    Descriptions,
    Tag,
    Space,
    Button,
    Row,
    Col,
    Statistic,
    Progress,
    Tabs,
    Table,
    Alert,
    Spin,
    message,
    Modal
} from 'antd';
import {
    ArrowLeftOutlined,
    RocketOutlined,
    ReloadOutlined,
    DownloadOutlined,
    EyeOutlined,
    ExperimentOutlined
} from '@ant-design/icons';
import { useParams, useNavigate } from 'react-router-dom';
import {
    getReactionNetworkJob,
    submitReactionNetworkJob,
    type ReactionNetworkJob,
    type Molecule,
    type Reaction
} from '../api/reactionNetwork';

const { TabPane } = Tabs;

const ReactionNetworkDetail: React.FC = () => {
    const { id } = useParams<{ id: string }>();
    const navigate = useNavigate();
    const [job, setJob] = useState<ReactionNetworkJob | null>(null);
    const [molecules, setMolecules] = useState<Molecule[]>([]);
    const [reactions, setReactions] = useState<Reaction[]>([]);
    const [loading, setLoading] = useState(true);

    useEffect(() => {
        loadJobDetail();
        // 自动刷新（如果任务在运行）
        const interval = setInterval(() => {
            if (job?.status === 'RUNNING' || job?.status === 'QUEUED') {
                loadJobDetail();
            }
        }, 10000);
        return () => clearInterval(interval);
    }, [id]);

    const loadJobDetail = async () => {
        try {
            const data = await getReactionNetworkJob(Number(id));
            setJob(data);
            setMolecules(data.molecules || []);
            setReactions(data.reactions || []);
        } catch (error) {
            message.error('加载任务详情失败');
        } finally {
            setLoading(false);
        }
    };

    const handleSubmit = async () => {
        try {
            await submitReactionNetworkJob(Number(id));
            message.success('任务已提交');
            loadJobDetail();
        } catch (error) {
            message.error('提交失败');
        }
    };

    const getStatusTag = (status: string) => {
        const config: Record<string, { color: string; text: string }> = {
            CREATED: { color: 'default', text: '已创建' },
            QUEUED: { color: 'processing', text: '排队中' },
            RUNNING: { color: 'processing', text: '运行中' },
            COMPLETED: { color: 'success', text: '已完成' },
            FAILED: { color: 'error', text: '失败' },
        };
        const cfg = config[status] || { color: 'default', text: status };
        return <Tag color={cfg.color}>{cfg.text}</Tag>;
    };

    const moleculeColumns = [
        { title: 'ID', dataIndex: 'id', key: 'id', width: 80 },
        { title: '名称', dataIndex: 'name', key: 'name' },
        { title: 'SMILES', dataIndex: 'smiles', key: 'smiles', render: (text: string) => <code>{text}</code> },
        { title: '代数', dataIndex: 'generation', key: 'generation', render: (gen: number) => <Tag color="blue">Gen {gen}</Tag> },
        {
            title: '能量 (kcal/mol)',
            dataIndex: 'energy_kcal',
            key: 'energy',
            render: (e: number) => e ? e.toFixed(2) : '-'
        },
        { title: '原子数', dataIndex: 'num_atoms', key: 'num_atoms' },
    ];

    const reactionColumns = [
        { title: 'ID', dataIndex: 'id', key: 'id', width: 80 },
        {
            title: '反应物',
            dataIndex: 'reactant_smiles',
            key: 'reactants',
            render: (smiles: string[]) => smiles.map((s, i) => <Tag key={i}>{s}</Tag>)
        },
        {
            title: '产物',
            dataIndex: 'product_smiles',
            key: 'products',
            render: (smiles: string[]) => smiles.map((s, i) => <Tag key={i} color="green">{s}</Tag>)
        },
        {
            title: '算符',
            dataIndex: 'operator_name',
            key: 'operator',
            render: (op: string) => <Tag color="orange">{op}</Tag>
        },
        {
            title: '反应能 (kcal/mol)',
            dataIndex: 'reaction_energy',
            key: 'energy',
            render: (e: number) => e ? e.toFixed(2) : '-'
        },
    ];

    if (loading) {
        return <div style={{ textAlign: 'center', padding: '100px' }}><Spin size="large" /></div>;
    }

    if (!job) {
        return <Alert message="任务不存在" type="error" />;
    }

    return (
        <div style={{ padding: '24px' }}>
            {/* 头部操作栏 */}
            <Card style={{ marginBottom: '16px' }}>
                <Row justify="space-between" align="middle">
                    <Col>
                        <Space size="large">
                            <Button icon={<ArrowLeftOutlined />} onClick={() => navigate('/reaction-network')}>
                                返回列表
                            </Button>
                            <div>
                                <Space>
                                    <ExperimentOutlined style={{ fontSize: '24px', color: '#1890ff' }} />
                                    <span style={{ fontSize: '20px', fontWeight: 'bold' }}>{job.job_name}</span>
                                    {getStatusTag(job.status)}
                                </Space>
                            </div>
                        </Space>
                    </Col>
                    <Col>
                        <Space>
                            <Button icon={<ReloadOutlined />} onClick={loadJobDetail}>刷新</Button>
                            {job.status === 'CREATED' && (
                                <Button type="primary" icon={<RocketOutlined />} onClick={handleSubmit}>
                                    提交运行
                                </Button>
                            )}
                            {job.status === 'COMPLETED' && job.network_json_path && (
                                <Button icon={<DownloadOutlined />}>下载结果</Button>
                            )}
                            {job.status === 'COMPLETED' && (
                                <Button
                                    type="primary"
                                    icon={<EyeOutlined />}
                                    onClick={() => navigate(`/reaction-network/${id}/visualize`)}
                                >
                                    查看网络图
                                </Button>
                            )}
                        </Space>
                    </Col>
                </Row>
            </Card>

            {/* 统计卡片 */}
            {job.status === 'COMPLETED' && (
                <Row gutter={16} style={{ marginBottom: '16px' }}>
                    <Col span={6}>
                        <Card>
                            <Statistic
                                title="生成分子数"
                                value={job.num_molecules || 0}
                                prefix={<ExperimentOutlined />}
                                valueStyle={{ color: '#1890ff' }}
                            />
                        </Card>
                    </Col>
                    <Col span={6}>
                        <Card>
                            <Statistic
                                title="发现反应数"
                                value={job.num_reactions || 0}
                                valueStyle={{ color: '#52c41a' }}
                            />
                        </Card>
                    </Col>
                    <Col span={6}>
                        <Card>
                            <Statistic
                                title="最大代数"
                                value={job.max_generation_reached || 0}
                                valueStyle={{ color: '#faad14' }}
                            />
                        </Card>
                    </Col>
                    <Col span={6}>
                        <Card>
                            <Statistic
                                title="CPU核时"
                                value={job.actual_cpu_hours?.toFixed(2) || 0}
                                suffix="h"
                            />
                        </Card>
                    </Col>
                </Row>
            )}

            {/* 进度显示 */}
            {(job.status === 'RUNNING' || job.status === 'QUEUED') && (
                <Card style={{ marginBottom: '16px' }}>
                    <Progress
                        percent={Math.round(job.progress * 100)}
                        status="active"
                        strokeColor={{ from: '#108ee9', to: '#87d068' }}
                    />
                    <p style={{ marginTop: '8px', color: '#666' }}>
                        任务正在运行中，请稍候...
                    </p>
                </Card>
            )}

            {/* 错误信息 */}
            {job.status === 'FAILED' && job.error_message && (
                <Alert
                    message="任务失败"
                    description={job.error_message}
                    type="error"
                    showIcon
                    style={{ marginBottom: '16px' }}
                />
            )}

            {/* 任务详情 */}
            <Card title="任务信息" style={{ marginBottom: '16px' }}>
                <Descriptions bordered column={2}>
                    <Descriptions.Item label="任务ID">{job.id}</Descriptions.Item>
                    <Descriptions.Item label="状态">{getStatusTag(job.status)}</Descriptions.Item>
                    <Descriptions.Item label="温度">{job.temperature} K</Descriptions.Item>
                    <Descriptions.Item label="电极类型">{job.electrode_type === 'anode' ? '阳极' : '阴极'}</Descriptions.Item>
                    <Descriptions.Item label="电压">{job.voltage} V</Descriptions.Item>
                    <Descriptions.Item label="最大代数">{job.max_generations}</Descriptions.Item>
                    <Descriptions.Item label="最大分子数">{job.max_species}</Descriptions.Item>
                    <Descriptions.Item label="能量截断">{job.energy_cutoff} kcal/mol</Descriptions.Item>
                    <Descriptions.Item label="创建时间" span={2}>
                        {new Date(job.created_at).toLocaleString('zh-CN')}
                    </Descriptions.Item>
                    {job.started_at && (
                        <Descriptions.Item label="开始时间" span={2}>
                            {new Date(job.started_at).toLocaleString('zh-CN')}
                        </Descriptions.Item>
                    )}
                    {job.finished_at && (
                        <Descriptions.Item label="完成时间" span={2}>
                            {new Date(job.finished_at).toLocaleString('zh-CN')}
                        </Descriptions.Item>
                    )}
                    <Descriptions.Item label="初始分子" span={2}>
                        <Space wrap>
                            {job.initial_smiles.map((s, i) => (
                                <Tag key={i} color="blue">{s}</Tag>
                            ))}
                        </Space>
                    </Descriptions.Item>
                    {job.description && (
                        <Descriptions.Item label="描述" span={2}>{job.description}</Descriptions.Item>
                    )}
                </Descriptions>
            </Card>

            {/* 结果标签页 */}
            {job.status === 'COMPLETED' && (
                <Card>
                    <Tabs defaultActiveKey="molecules">
                        <TabPane tab={`分子 (${molecules.length})`} key="molecules">
                            <Table
                                columns={moleculeColumns}
                                dataSource={molecules}
                                rowKey="id"
                                pagination={{ pageSize: 20 }}
                            />
                        </TabPane>
                        <TabPane tab={`反应 (${reactions.length})`} key="reactions">
                            <Table
                                columns={reactionColumns}
                                dataSource={reactions}
                                rowKey="id"
                                pagination={{ pageSize: 20 }}
                            />
                        </TabPane>
                    </Tabs>
                </Card>
            )}
        </div>
    );
};

export default ReactionNetworkDetail;
