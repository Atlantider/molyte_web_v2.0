/**
 * Reaction Network Detail Page
 * 任务详情页 - 增强版，显示更多化学和电化学信息
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
    Tooltip
} from 'antd';
import {
    ArrowLeftOutlined,
    RocketOutlined,
    ReloadOutlined,
    DownloadOutlined,
    EyeOutlined,
    ExperimentOutlined,
    ThunderboltOutlined,
    FireOutlined,
    InfoCircleOutlined,
    BulbOutlined
} from '@ant-design/icons';
import { useParams, useNavigate } from 'react-router-dom';
import {
    getReactionNetworkJob,
    submitReactionNetworkJob,
    getActivatedOperators,
    type ReactionNetworkJob,
    type Molecule,
    type Reaction,
    type ActivatedOperator,
    type OperatorsResponse
} from '../api/reactionNetwork';

const { TabPane } = Tabs;

// 材料信息库
const ANODE_MATERIALS: Record<string, { name: string; carrier: string; description: string; features: string[] }> = {
    'GRAPHITE': {
        name: '石墨',
        carrier: 'Li+',
        description: '层状结构锂离子电池负极材料',
        features: ['成熟稳定', 'SEI膜形成', '可逆嵌锂']
    },
    'LI_METAL': {
        name: '锂金属',
        carrier: 'Li+',
        description: '最高理论容量负极材料',
        features: ['3860 mAh/g', '枝晶风险', '需要稳定SEI']
    },
    'SILICON': {
        name: '硅负极',
        carrier: 'Li+',
        description: '高容量下一代负极',
        features: ['~3600 mAh/g', '体积膨胀大', '需要Si-O-Li SEI']
    },
    'HARD_CARBON': {
        name: '硬碳',
        carrier: 'Na+',
        description: '钠离子电池主流负极',
        features: ['无序结构', '钠离子存储', '~300 mAh/g']
    },
    'NA_METAL': {
        name: '钠金属',
        carrier: 'Na+',
        description: '钠离子电池金属负极',
        features: ['高容量', '低成本', 'SEI形成']
    },
    'K_METAL': {
        name: '钾金属',
        carrier: 'K+',
        description: '钾离子电池金属负极',
        features: ['低电位', '快速动力学']
    },
};

const CATHODE_MATERIALS: Record<string, { name: string; description: string; voltage: string; features: string[] }> = {
    'NMC': {
        name: 'NMC三元',
        description: '镍锰钴三元正极',
        voltage: '3.6-3.8V',
        features: ['平衡性能', '中等电压', '适度氧释放']
    },
    'NMC811': {
        name: 'NMC811高镍',
        description: '高镍三元正极材料',
        voltage: '3.7-4.3V',
        features: ['高能量密度', '高电压氧释放', 'Ni4+活性']
    },
    'LFP': {
        name: '磷酸铁锂',
        description: '最稳定正极材料',
        voltage: '3.2-3.4V',
        features: ['无氧释放', '热稳定', '长寿命']
    },
    'LCO': {
        name: '钴酸锂',
        description: '高电压层状正极',
        voltage: '3.7-4.2V',
        features: ['高电压', '氧释放', '消费电子']
    },
};

// 识别物种类型
const analyzeSpecies = (smiles: string): { type: string; name: string; color: string } => {
    // 阳离子
    if (smiles === '[Li+]') return { type: '载流子', name: 'Li⁺', color: 'blue' };
    if (smiles === '[Na+]') return { type: '载流子', name: 'Na⁺', color: 'blue' };
    if (smiles === '[K+]') return { type: '载流子', name: 'K⁺', color: 'blue' };
    if (smiles === '[Li]') return { type: '金属原子', name: 'Li', color: 'cyan' };

    // 常见阴离子
    if (smiles.includes('PF6') || smiles.includes('P-')) return { type: '阴离子', name: 'PF₆⁻', color: 'red' };
    if (smiles.includes('TFSI') || smiles.includes('N-]S')) return { type: '阴离子', name: 'TFSI⁻', color: 'red' };
    if (smiles.includes('FSI')) return { type: '阴离子', name: 'FSI⁻', color: 'red' };

    // 溶剂
    if (smiles === 'C1COC(=O)O1') return { type: '溶剂', name: 'EC', color: 'green' };
    if (smiles === 'COC(=O)OC') return { type: '溶剂', name: 'DMC', color: 'green' };
    if (smiles === 'CCOC(=O)OCC') return { type: '溶剂', name: 'DEC', color: 'green' };
    if (smiles.includes('OCCOC')) return { type: '溶剂', name: 'DME', color: 'green' };

    // 其他
    return { type: '其他分子', name: smiles.substring(0, 20), color: 'default' };
};

const ReactionNetworkDetail: React.FC = () => {
    const { id } = useParams<{ id: string }>();
    const navigate = useNavigate();
    const [job, setJob] = useState<ReactionNetworkJob | null>(null);
    const [molecules, setMolecules] = useState<Molecule[]>([]);
    const [reactions, setReactions] = useState<Reaction[]>([]);
    const [loading, setLoading] = useState(true);
    const [operatorInfo, setOperatorInfo] = useState<OperatorsResponse | null>(null);
    const [loadingOperators, setLoadingOperators] = useState(false);

    useEffect(() => {
        loadJobDetail();
        loadOperatorInfo();
        const interval = setInterval(() => {
            if (job?.status === 'RUNNING' || job?.status === 'QUEUED') {
                loadJobDetail();
            }
        }, 10000);
        return () => clearInterval(interval);
    }, [id]);

    const loadOperatorInfo = async () => {
        setLoadingOperators(true);
        try {
            const data = await getActivatedOperators(Number(id));
            setOperatorInfo(data);
        } catch (error) {
            console.error('Failed to load operator info:', error);
        } finally {
            setLoadingOperators(false);
        }
    };

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

    // 统计算符使用情况
    const getOperatorStats = () => {
        const stats: Record<string, number> = {};
        reactions.forEach(r => {
            const op = r.operator_name || '未知';
            stats[op] = (stats[op] || 0) + 1;
        });
        return Object.entries(stats)
            .map(([name, count]) => ({ operator: name, count, percentage: ((count / reactions.length) * 100).toFixed(1) }))
            .sort((a, b) => b.count - a.count);
    };

    // 分析初始物种
    const getInitialSpeciesAnalysis = () => {
        if (!job) return [];
        return job.initial_smiles.map(smiles => ({
            smiles,
            ...analyzeSpecies(smiles)
        }));
    };

    const moleculeColumns = [
        { title: 'ID', dataIndex: 'id', key: 'id', width: 80 },
        { title: '名称', dataIndex: 'name', key: 'name' },
        { title: 'SMILES', dataIndex: 'smiles', key: 'smiles', render: (text: string) => <code style={{ fontSize: 12 }}>{text}</code> },
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
            render: (smiles: string[]) => <Space wrap>{smiles.map((s, i) => <Tag key={i} style={{ fontSize: 11 }}>{s.substring(0, 30)}</Tag>)}</Space>
        },
        {
            title: '产物',
            dataIndex: 'product_smiles',
            key: 'products',
            render: (smiles: string[]) => <Space wrap>{smiles.map((s, i) => <Tag key={i} color="green" style={{ fontSize: 11 }}>{s.substring(0, 30)}</Tag>)}</Space>
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

    // 尝试从job.config读取材料信息（如果config存在且为对象）
    const anodeMat = (job as any).config?.anode_material || (job as any).anode_material;
    const cathodeMat = (job as any).config?.cathode_material || (job as any).cathode_material;
    const anodeMaterial = anodeMat ? ANODE_MATERIALS[anodeMat] : null;
    const cathodeMaterial = cathodeMat ? CATHODE_MATERIALS[cathodeMat] : null;
    const operatorStats = reactions.length > 0 ? getOperatorStats() : [];
    const speciesAnalysis = getInitialSpeciesAnalysis();

    return (
        <div style={{ padding: '24px' }}>
            {/* 头部操作栏 */}
            <Card style={{ marginBottom: '16px' }}>
                <Row justify="space-between" align="middle">
                    <Col>
                        <Space size="large">
                            <Button icon={<ArrowLeftOutlined />} onClick={() => navigate('/workspace/liquid-electrolyte/reaction-network')}>
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
                                    onClick={() => navigate(`/workspace/liquid-electrolyte/reaction-network/${id}/visualize`)}
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
                                title="计算耗时"
                                value={(job as any).actual_cpu_hours?.toFixed(2) || '-'}
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

            {/* 电化学环境 */}
            <Card
                title={
                    <Space>
                        <ThunderboltOutlined style={{ color: '#1890ff' }} />
                        <span>电化学环境</span>
                    </Space>
                }
                style={{ marginBottom: '16px' }}
            >
                <Descriptions bordered column={3} size="small">
                    <Descriptions.Item label="温度">
                        {job.temperature} K ({(job.temperature - 273.15).toFixed(1)}°C)
                    </Descriptions.Item>
                    <Descriptions.Item label="电极类型">
                        <Tag color={job.electrode_type === 'anode' ? 'blue' : 'red'}>
                            {job.electrode_type === 'anode' ? '阳极 (负极)' : '阴极 (正极)'}
                        </Tag>
                    </Descriptions.Item>
                    <Descriptions.Item label="电压">
                        <span style={{ fontWeight: 'bold', color: '#1890ff' }}>{job.voltage} V</span>
                    </Descriptions.Item>

                    {anodeMaterial && (
                        <>
                            <Descriptions.Item label="负极材料" span={3}>
                                <Space>
                                    <Tag color="blue">{anodeMaterial.name}</Tag>
                                    <Tag color="cyan">载流子: {anodeMaterial.carrier}</Tag>
                                    <span style={{ color: '#666', fontSize: 12 }}>{anodeMaterial.description}</span>
                                </Space>
                                <div style={{ marginTop: 8 }}>
                                    {anodeMaterial.features.map((f, i) => (
                                        <Tag key={i} color="processing" style={{ marginTop: 4 }}>{f}</Tag>
                                    ))}
                                </div>
                            </Descriptions.Item>
                        </>
                    )}

                    {cathodeMaterial && (
                        <>
                            <Descriptions.Item label="正极材料" span={3}>
                                <Space>
                                    <Tag color="red">{cathodeMaterial.name}</Tag>
                                    <Tag color="orange">电压: {cathodeMaterial.voltage}</Tag>
                                    <span style={{ color: '#666', fontSize: 12 }}>{cathodeMaterial.description}</span>
                                </Space>
                                <div style={{ marginTop: 8 }}>
                                    {cathodeMaterial.features.map((f, i) => (
                                        <Tag key={i} color="warning" style={{ marginTop: 4 }}>{f}</Tag>
                                    ))}
                                </div>
                            </Descriptions.Item>
                        </>
                    )}
                </Descriptions>
            </Card>

            {/* 初始物种分析 */}
            <Card
                title={
                    <Space>
                        <BulbOutlined style={{ color: '#52c41a' }} />
                        <span>初始物种组成</span>
                    </Space>
                }
                style={{ marginBottom: '16px' }}
            >
                <Alert
                    message="系统将从以下初始分子开始生成反应网络"
                    type="info"
                    showIcon
                    style={{ marginBottom: 12 }}
                />
                <Table
                    size="small"
                    pagination={false}
                    columns={[
                        { title: '类型', dataIndex: 'type', key: 'type', render: (t: string) => <Tag color="purple">{t}</Tag> },
                        { title: '名称', dataIndex: 'name', key: 'name' },
                        { title: 'SMILES', dataIndex: 'smiles', key: 'smiles', render: (s: string) => <code style={{ fontSize: 11 }}>{s}</code> }
                    ]}
                    dataSource={speciesAnalysis}
                />
            </Card>


            {/* 激活的算符详情 */}
            {operatorInfo && operatorInfo.operators.length > 0 && (
                <Card
                    title={
                        <Space>
                            <FireOutlined style={{ color: '#faad14' }} />
                            <span>激活的反应算符</span>
                            <Tag color="cyan">{operatorInfo.num_operators} 个算符</Tag>
                            <Tooltip title="根据当前电化学环境自动激活的反应算符">
                                <InfoCircleOutlined style={{ color: '#999' }} />
                            </Tooltip>
                        </Space>
                    }
                    style={{ marginBottom: '16px' }}
                    loading={loadingOperators}
                >
                    <Alert
                        message="算符激活说明"
                        description={
                            <div style={{ fontSize: 12 }}>
                                <div>✓ 根据温度 {operatorInfo.environment.temperature}K、电压 {operatorInfo.environment.voltage}V 自动激活相应算符</div>
                                <div>✓ 不同算符对应不同的化学反应机理（电子转移、质子转移、自由基反应等）</div>
                                <div>✓ 激活驱动力: {operatorInfo.environment.active_drives.join(', ')}</div>
                            </div>
                        }
                        type="info"
                        showIcon
                        style={{ marginBottom: 16 }}
                    />

                    <Table
                        size="small"
                        pagination={{ pageSize: 10 }}
                        expandable={{
                            expandedRowRender: (record: ActivatedOperator) => (
                                <div style={{ padding: '12px', backgroundColor: '#fafafa' }}>
                                    <Row gutter={16}>
                                        <Col span={12}>
                                            <div style={{ marginBottom: 8 }}>
                                                <span style={{ fontWeight: 'bold', color: '#1890ff' }}>描述：</span>
                                                <div style={{ marginTop: 4 }}>{record.description}</div>
                                            </div>
                                            <div style={{ marginBottom: 8 }}>
                                                <span style={{ fontWeight: 'bold', color: '#52c41a' }}>激活原因：</span>
                                                <div style={{ marginTop: 4 }}>
                                                    {record.activation_reasons.map((reason, i) => (
                                                        <Tag key={i} color="green" style={{ marginTop: 4 }}>{reason}</Tag>
                                                    ))}
                                                </div>
                                            </div>
                                        </Col>
                                        <Col span={12}>
                                            {record.conditions.length > 0 && (
                                                <div style={{ marginBottom: 8 }}>
                                                    <span style={{ fontWeight: 'bold', color: '#faad14' }}>适用条件：</span>
                                                    <div style={{ marginTop: 4 }}>
                                                        {record.conditions.map((cond, i) => (
                                                            <div key={i} style={{ fontSize: 12, marginTop: 2 }}>• {cond}</div>
                                                        ))}
                                                    </div>
                                                </div>
                                            )}
                                            {record.molecular_checks.length > 0 && (
                                                <div>
                                                    <span style={{ fontWeight: 'bold', color: '#722ed1' }}>分子特征要求：</span>
                                                    <div style={{ marginTop: 4 }}>
                                                        {record.molecular_checks.map((check, i) => (
                                                            <Tag key={i} color="purple" style={{ marginTop: 4 }}>{check}</Tag>
                                                        ))}
                                                    </div>
                                                </div>
                                            )}
                                        </Col>
                                    </Row>
                                </div>
                            ),
                        }}
                        columns={[
                            {
                                title: '算符名称',
                                dataIndex: 'name',
                                key: 'name',
                                render: (name: string) => <Tag color="orange" style={{ fontSize: 13 }}>{name}</Tag>
                            },
                            {
                                title: '权重',
                                dataIndex: 'weight',
                                key: 'weight',
                                width: 80,
                                render: (w: number) => (w * 10).toFixed(1),
                                sorter: (a, b) => a.weight - b.weight
                            },
                            {
                                title: '必需驱动力',
                                dataIndex: 'required_drives',
                                key: 'required_drives',
                                render: (drives: string[]) => drives.length > 0 ? (
                                    <Space wrap>
                                        {drives.map((d, i) => <Tag key={i} color="red">{d}</Tag>)}
                                    </Space>
                                ) : <span style={{ color: '#999' }}>无</span>
                            },
                            {
                                title: '增强驱动力',
                                dataIndex: 'enhancing_drives',
                                key: 'enhancing_drives',
                                render: (drives: string[]) => drives.length > 0 ? (
                                    <Space wrap>
                                        {drives.map((d, i) => <Tag key={i} color="blue">{d}</Tag>)}
                                    </Space>
                                ) : <span style={{ color: '#999' }}>无</span>
                            },
                        ]}
                        dataSource={operatorInfo.operators}
                        rowKey="name"
                    />
                </Card>
            )}

            {/* 算符使用统计 */}
            {job.status === 'COMPLETED' && operatorStats.length > 0 && (
                <Card
                    title={
                        <Space>
                            <FireOutlined style={{ color: '#faad14' }} />
                            <span>算符使用统计</span>
                            <Tooltip title="统计实际生成的反应中各算符的使用情况">
                                <InfoCircleOutlined style={{ color: '#999' }} />
                            </Tooltip>
                        </Space>
                    }
                    style={{ marginBottom: '16px' }}
                >
                    <Row gutter={16}>
                        <Col span={12}>
                            <Table
                                size="small"
                                pagination={false}
                                columns={[
                                    { title: '算符名称', dataIndex: 'operator', key: 'operator', render: (o: string) => <Tag color="orange">{o}</Tag> },
                                    { title: '使用次数', dataIndex: 'count', key: 'count' },
                                    { title: '占比', dataIndex: 'percentage', key: 'percentage', render: (p: string) => `${p}%` }
                                ]}
                                dataSource={operatorStats}
                            />
                        </Col>
                        <Col span={12}>
                            <Alert
                                message="统计说明"
                                description={
                                    <div style={{ fontSize: 12 }}>
                                        <div>✓ 统计来自实际生成的 {reactions.length} 个反应</div>
                                        <div>✓ 使用频率反映该算符在网络中的重要性</div>
                                        <div>✓ 某些激活的算符可能未被使用（无匹配分子）</div>
                                    </div>
                                }
                                type="info"
                                showIcon
                            />
                        </Col>
                    </Row>
                </Card>
            )}


            {/* 任务详情 */}
            <Card title="任务配置参数" style={{ marginBottom: '16px' }}>
                <Descriptions bordered column={2} size="small">
                    <Descriptions.Item label="任务ID">{job.id}</Descriptions.Item>
                    <Descriptions.Item label="状态">{getStatusTag(job.status)}</Descriptions.Item>
                    <Descriptions.Item label="最大代数">{job.max_generations}</Descriptions.Item>
                    <Descriptions.Item label="最大分子数">{job.max_species}</Descriptions.Item>
                    <Descriptions.Item label="能量截断">{job.energy_cutoff} kcal/mol</Descriptions.Item>
                    <Descriptions.Item label="Slurm分区">{(job as any).slurm_partition || 'cpu'}</Descriptions.Item>
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
                                scroll={{ x: 1000 }}
                            />
                        </TabPane>
                        <TabPane tab={`反应 (${reactions.length})`} key="reactions">
                            <Table
                                columns={reactionColumns}
                                dataSource={reactions}
                                rowKey="id"
                                pagination={{ pageSize: 20 }}
                                scroll={{ x: 1200 }}
                            />
                        </TabPane>
                    </Tabs>
                </Card>
            )}
        </div>
    );
};

export default ReactionNetworkDetail;
