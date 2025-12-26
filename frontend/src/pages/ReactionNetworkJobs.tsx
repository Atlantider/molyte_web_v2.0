/**
 * Reaction Network Jobs List Page
 * 反应网络任务列表 - 现代化设计
 */

import React, { useState, useEffect } from 'react';
import {
    Table,
    Button,
    Tag,
    Space,
    Card,
    Statistic,
    Row,
    Col,
    Input,
    Select,
    message,
    Popconfirm,
    Progress,
    Tooltip,
    Empty
} from 'antd';
import {
    PlusOutlined,
    EyeOutlined,
    DeleteOutlined,
    ReloadOutlined,
    SearchOutlined,
    RocketOutlined,
    ExperimentOutlined
} from '@ant-design/icons';
import { useNavigate } from 'react-router-dom';
import {
    getReactionNetworkJobs,
    submitReactionNetworkJob,
    deleteReactionNetworkJob,
    type ReactionNetworkJob
} from '../api/reactionNetwork';

const { Option } = Select;

const ReactionNetworkJobs: React.FC = () => {
    const navigate = useNavigate();
    const [jobs, setJobs] = useState<ReactionNetworkJob[]>([]);
    const [loading, setLoading] = useState(false);
    const [total, setTotal] = useState(0);
    const [selectedStatus, setSelectedStatus] = useState<string | undefined>();
    const [searchText, setSearchText] = useState('');
    const [pagination, setPagination] = useState({ current: 1, pageSize: 20 });

    // 统计数据
    const [stats, setStats] = useState({
        total: 0,
        running: 0,
        completed: 0,
        failed: 0
    });

    useEffect(() => {
        loadJobs();
        // 自动刷新（每30秒）
        const interval = setInterval(loadJobs, 30000);
        return () => clearInterval(interval);
    }, [selectedStatus, pagination.current, pagination.pageSize]);

    const loadJobs = async () => {
        setLoading(true);
        try {
            const skip = (pagination.current - 1) * pagination.pageSize;
            const data = await getReactionNetworkJobs({
                status: selectedStatus,
                skip,
                limit: pagination.pageSize,
                sort_by: 'created_at',
                sort_desc: true
            });

            setJobs(data.jobs);
            setTotal(data.total);

            // 更新统计
            const newStats = {
                total: data.total,
                running: data.jobs.filter(j => j.status === 'RUNNING' || j.status === 'QUEUED').length,
                completed: data.jobs.filter(j => j.status === 'COMPLETED').length,
                failed: data.jobs.filter(j => j.status === 'FAILED').length
            };
            setStats(newStats);
        } catch (error) {
            message.error('加载任务列表失败');
        } finally {
            setLoading(false);
        }
    };

    const handleSubmit = async (jobId: number) => {
        try {
            await submitReactionNetworkJob(jobId);
            message.success('任务已提交到计算队列');
            loadJobs();
        } catch (error) {
            message.error('提交任务失败');
        }
    };

    const handleDelete = async (jobId: number) => {
        try {
            await deleteReactionNetworkJob(jobId);
            message.success('任务已删除');
            loadJobs();
        } catch (error) {
            message.error('删除任务失败');
        }
    };

    const getStatusTag = (status: string) => {
        const statusConfig: Record<string, { color: string; text: string }> = {
            CREATED: { color: 'default', text: '已创建' },
            QUEUED: { color: 'processing', text: '排队中' },
            RUNNING: { color: 'processing', text: '运行中' },
            POSTPROCESSING: { color: 'processing', text: '后处理' },
            COMPLETED: { color: 'success', text: '已完成' },
            FAILED: { color: 'error', text: '失败' },
            CANCELLED: { color: 'default', text: '已取消' }
        };
        const config = statusConfig[status] || { color: 'default', text: status };
        return <Tag color={config.color}>{config.text}</Tag>;
    };

    const columns = [
        {
            title: 'ID',
            dataIndex: 'id',
            key: 'id',
            width: 80,
            fixed: 'left' as const,
        },
        {
            title: '任务名称',
            dataIndex: 'job_name',
            key: 'job_name',
            width: 200,
            fixed: 'left' as const,
            render: (text: string, record: ReactionNetworkJob) => (
                <a onClick={() => navigate(`/reaction-network/${record.id}`)}>
                    <ExperimentOutlined style={{ marginRight: 8 }} />
                    {text}
                </a>
            ),
        },
        {
            title: '状态',
            dataIndex: 'status',
            key: 'status',
            width: 120,
            render: (status: string) => getStatusTag(status),
        },
        {
            title: '进度',
            dataIndex: 'progress',
            key: 'progress',
            width: 150,
            render: (progress: number, record: ReactionNetworkJob) => {
                if (record.status === 'COMPLETED') {
                    return <Progress percent={100} size="small" status="success" />;
                }
                if (record.status === 'FAILED') {
                    return <Progress percent={100} size="small" status="exception" />;
                }
                if (record.status === 'RUNNING' || record.status === 'POSTPROCESSING') {
                    return <Progress percent={Math.round(progress * 100)} size="small" status="active" />;
                }
                return <Progress percent={0} size="small" />;
            },
        },
        {
            title: '初始分子',
            dataIndex: 'initial_smiles',
            key: 'initial_smiles',
            width: 100,
            render: (smiles: string[]) => (
                <Tooltip title={smiles.join(', ')}>
                    <Tag color="blue">{smiles.length} 个</Tag>
                </Tooltip>
            ),
        },
        {
            title: '温度',
            dataIndex: 'temperature',
            key: 'temperature',
            width: 100,
            render: (temp: number) => `${temp} K`,
        },
        {
            title: '最大代数',
            dataIndex: 'max_generations',
            key: 'max_generations',
            width: 100,
        },
        {
            title: '结果',
            key: 'results',
            width: 180,
            render: (_: any, record: ReactionNetworkJob) => (
                <Space>
                    {record.num_molecules !== null && record.num_molecules !== undefined && (
                        <Tag color="green">{record.num_molecules} 分子</Tag>
                    )}
                    {record.num_reactions !== null && record.num_reactions !== undefined && (
                        <Tag color="orange">{record.num_reactions} 反应</Tag>
                    )}
                </Space>
            ),
        },
        {
            title: '创建时间',
            dataIndex: 'created_at',
            key: 'created_at',
            width: 180,
            render: (date: string) => new Date(date).toLocaleString('zh-CN'),
        },
        {
            title: '操作',
            key: 'actions',
            width: 200,
            fixed: 'right' as const,
            render: (_: any, record: ReactionNetworkJob) => (
                <Space>
                    <Button
                        type="link"
                        size="small"
                        icon={<EyeOutlined />}
                        onClick={() => navigate(`/reaction-network/${record.id}`)}
                    >
                        查看
                    </Button>
                    {record.status === 'CREATED' && (
                        <Button
                            type="link"
                            size="small"
                            icon={<RocketOutlined />}
                            onClick={() => handleSubmit(record.id)}
                        >
                            提交
                        </Button>
                    )}
                    <Popconfirm
                        title="确定要删除这个任务吗？"
                        onConfirm={() => handleDelete(record.id)}
                        okText="确定"
                        cancelText="取消"
                    >
                        <Button
                            type="link"
                            size="small"
                            danger
                            icon={<DeleteOutlined />}
                        >
                            删除
                        </Button>
                    </Popconfirm>
                </Space>
            ),
        },
    ];

    return (
        <div style={{ padding: '24px' }}>
            {/* 统计卡片 */}
            <Row gutter={16} style={{ marginBottom: '24px' }}>
                <Col span={6}>
                    <Card>
                        <Statistic
                            title="总任务数"
                            value={stats.total}
                            prefix={<ExperimentOutlined />}
                            valueStyle={{ color: '#1890ff' }}
                        />
                    </Card>
                </Col>
                <Col span={6}>
                    <Card>
                        <Statistic
                            title="运行中"
                            value={stats.running}
                            prefix={<RocketOutlined />}
                            valueStyle={{ color: '#52c41a' }}
                        />
                    </Card>
                </Col>
                <Col span={6}>
                    <Card>
                        <Statistic
                            title="已完成"
                            value={stats.completed}
                            valueStyle={{ color: '#52c41a' }}
                        />
                    </Card>
                </Col>
                <Col span={6}>
                    <Card>
                        <Statistic
                            title="失败"
                            value={stats.failed}
                            valueStyle={{ color: '#ff4d4f' }}
                        />
                    </Card>
                </Col>
            </Row>

            {/* 操作栏 */}
            <Card style={{ marginBottom: '16px' }}>
                <Row gutter={16} align="middle">
                    <Col flex="auto">
                        <Space>
                            <Button
                                type="primary"
                                icon={<PlusOutlined />}
                                onClick={() => navigate('/reaction-network/create')}
                                size="large"
                            >
                                创建反应网络任务
                            </Button>
                            <Button
                                icon={<ReloadOutlined />}
                                onClick={loadJobs}
                                loading={loading}
                            >
                                刷新
                            </Button>
                        </Space>
                    </Col>
                    <Col>
                        <Space>
                            <Input
                                placeholder="搜索任务名称"
                                prefix={<SearchOutlined />}
                                value={searchText}
                                onChange={(e) => setSearchText(e.target.value)}
                                style={{ width: 200 }}
                            />
                            <Select
                                placeholder="筛选状态"
                                allowClear
                                value={selectedStatus}
                                onChange={setSelectedStatus}
                                style={{ width: 150 }}
                            >
                                <Option value="CREATED">已创建</Option>
                                <Option value="QUEUED">排队中</Option>
                                <Option value="RUNNING">运行中</Option>
                                <Option value="COMPLETED">已完成</Option>
                                <Option value="FAILED">失败</Option>
                            </Select>
                        </Space>
                    </Col>
                </Row>
            </Card>

            {/* 任务表格 */}
            <Card>
                <Table
                    columns={columns}
                    dataSource={jobs.filter(job =>
                        searchText ? job.job_name.toLowerCase().includes(searchText.toLowerCase()) : true
                    )}
                    loading={loading}
                    rowKey="id"
                    pagination={{
                        ...pagination,
                        total,
                        showSizeChanger: true,
                        showQuickJumper: true,
                        showTotal: (total) => `共 ${total} 条`,
                        onChange: (page, pageSize) => {
                            setPagination({ current: page, pageSize: pageSize || 20 });
                        },
                    }}
                    scroll={{ x: 1500 }}
                    locale={{
                        emptyText: (
                            <Empty
                                image={Empty.PRESENTED_IMAGE_SIMPLE}
                                description="暂无任务，点击上方按钮创建第一个反应网络任务"
                            />
                        ),
                    }}
                />
            </Card>
        </div>
    );
};

export default ReactionNetworkJobs;
