/**
 * Reaction Network Comparison Component
 * 反应网络对比分析功能
 */

import React, { useState, useEffect } from 'react';
import {
    Card,
    Select,
    Button,
    Table,
    Row,
    Col,
    Statistic,
    Space,
    Tag,
    Divider,
    Alert,
    Tabs,
    Progress
} from 'antd';
import {
    SwapOutlined,
    BarChartOutlined
} from '@ant-design/icons';
import {
    BarChart,
    Bar,
    LineChart,
    Line,
    PieChart,
    Pie,
    Cell,
    XAxis,
    YAxis,
    CartesianGrid,
    Tooltip,
    Legend,
    ResponsiveContainer
} from 'recharts';
import { getReactionNetworkJob, type ReactionNetworkJob } from '../api/reactionNetwork';

const { Option } = Select;
const { TabPane } = Tabs;

const ReactionNetworkComparison: React.FC = () => {
    const [job1Id, setJob1Id] = useState<number | null>(null);
    const [job2Id, setJob2Id] = useState<number | null>(null);
    const [job1, setJob1] = useState<any>(null);
    const [job2, setJob2] = useState<any>(null);
    const [loading, setLoading] = useState(false);
    const [availableJobs, setAvailableJobs] = useState<ReactionNetworkJob[]>([]);

    useEffect(() => {
        // 加载可用的已完成任务
        loadAvailableJobs();
    }, []);

    const loadAvailableJobs = async () => {
        // TODO: 从API加载已完成的任务列表
    };

    const handleCompare = async () => {
        if (!job1Id || !job2Id) return;

        setLoading(true);
        try {
            const [data1, data2] = await Promise.all([
                getReactionNetworkJob(job1Id),
                getReactionNetworkJob(job2Id)
            ]);
            setJob1(data1);
            setJob2(data2);
        } catch (error) {
            console.error('Failed to load jobs for comparison', error);
        } finally {
            setLoading(false);
        }
    };

    // 对比统计
    const getComparisonStats = () => {
        if (!job1 || !job2) return [];

        return [
            {
                metric: '分子数',
                job1: job1.num_molecules,
                job2: job2.num_molecules,
                diff: ((job1.num_molecules - job2.num_molecules) / job2.num_molecules * 100).toFixed(1) + '%'
            },
            {
                metric: '反应数',
                job1: job1.num_reactions,
                job2: job2.num_reactions,
                diff: ((job1.num_reactions - job2.num_reactions) / job2.num_reactions * 100).toFixed(1) + '%'
            },
            {
                metric: '最大代数',
                job1: job1.max_generation_reached,
                job2: job2.max_generation_reached,
                diff: job1.max_generation_reached - job2.max_generation_reached
            },
            {
                metric: 'CPU核时',
                job1: job1.actual_cpu_hours?.toFixed(2),
                job2: job2.actual_cpu_hours?.toFixed(2),
                diff: ((job1.actual_cpu_hours - job2.actual_cpu_hours) / job2.actual_cpu_hours * 100).toFixed(1) + '%'
            }
        ];
    };

    // 代数分布对比
    const getGenerationDistribution = () => {
        if (!job1 || !job2) return [];

        const maxGen = Math.max(job1.max_generation_reached, job2.max_generation_reached);
        const data = [];

        for (let i = 0; i <= maxGen; i++) {
            data.push({
                generation: `Gen ${i}`,
                job1: job1.molecules?.filter((m: any) => m.generation === i).length || 0,
                job2: job2.molecules?.filter((m: any) => m.generation === i).length || 0
            });
        }

        return data;
    };

    // 能量分布对比
    const getEnergyDistribution = () => {
        if (!job1 || !job2) return [];

        const bins = [
            { range: '< -50', min: -Infinity, max: -50 },
            { range: '-50~0', min: -50, max: 0 },
            { range: '0~50', min: 0, max: 50 },
            { range: '> 50', min: 50, max: Infinity }
        ];

        return bins.map(bin => ({
            range: bin.range,
            job1: job1.reactions?.filter((r: any) =>
                r.reaction_energy >= bin.min && r.reaction_energy < bin.max
            ).length || 0,
            job2: job2.reactions?.filter((r: any) =>
                r.reaction_energy >= bin.min && r.reaction_energy < bin.max
            ).length || 0
        }));
    };

    const columns = [
        { title: '指标', dataIndex: 'metric', key: 'metric', width: 150 },
        {
            title: job1?.job_name || '任务1',
            dataIndex: 'job1',
            key: 'job1',
            render: (val: any) => <strong>{val}</strong>
        },
        {
            title: job2?.job_name || '任务2',
            dataIndex: 'job2',
            key: 'job2',
            render: (val: any) => <strong>{val}</strong>
        },
        {
            title: '差异',
            dataIndex: 'diff',
            key: 'diff',
            render: (val: any) => {
                const isPositive = String(val).includes('+') || parseFloat(val) > 0;
                return (
                    <Tag color={isPositive ? 'green' : 'red'}>
                        {isPositive && '+'}{val}
                    </Tag>
                );
            }
        }
    ];

    const COLORS = ['#0088FE', '#00C49F', '#FFBB28', '#FF8042'];

    return (
        <div style={{ padding: '24px' }}>
            <Card title={<Space><BarChartOutlined />反应网络对比分析</Space>}>
                <Alert
                    message="选择两个已完成的反应网络任务进行对比分析"
                    type="info"
                    showIcon
                    style={{ marginBottom: '24px' }}
                />

                {/* 任务选择 */}
                <Row gutter={16} style={{ marginBottom: '24px' }}>
                    <Col span={10}>
                        <div style={{ marginBottom: 8 }}>任务1:</div>
                        <Select
                            value={job1Id}
                            onChange={setJob1Id}
                            style={{ width: '100%' }}
                            placeholder="选择任务1"
                            showSearch
                            filterOption={(input, option) =>
                                String(option?.children).toLowerCase().includes(input.toLowerCase())
                            }
                        >
                            {availableJobs.map(job => (
                                <Option key={job.id} value={job.id}>
                                    {job.job_name} (ID: {job.id})
                                </Option>
                            ))}
                        </Select>
                    </Col>

                    <Col span={4} style={{ textAlign: 'center', paddingTop: '32px' }}>
                        <SwapOutlined style={{ fontSize: '24px', color: '#1890ff' }} />
                    </Col>

                    <Col span={10}>
                        <div style={{ marginBottom: 8 }}>任务2:</div>
                        <Select
                            value={job2Id}
                            onChange={setJob2Id}
                            style={{ width: '100%' }}
                            placeholder="选择任务2"
                            showSearch
                            filterOption={(input, option) =>
                                String(option?.children).toLowerCase().includes(input.toLowerCase())
                            }
                        >
                            {availableJobs.map(job => (
                                <Option key={job.id} value={job.id}>
                                    {job.job_name} (ID: {job.id})
                                </Option>
                            ))}
                        </Select>
                    </Col>
                </Row>

                <Button
                    type="primary"
                    icon={<BarChartOutlined />}
                    onClick={handleCompare}
                    disabled={!job1Id || !job2Id}
                    loading={loading}
                    block
                    size="large"
                >
                    开始对比
                </Button>

                {job1 && job2 && (
                    <>
                        <Divider />

                        {/* 对比统计表 */}
                        <Card
                            title="📊 统计对比"
                            style={{ marginBottom: '24px' }}
                            type="inner"
                        >
                            <Table
                                columns={columns}
                                dataSource={getComparisonStats()}
                                pagination={false}
                                size="middle"
                            />
                        </Card>

                        {/* 图表对比 */}
                        <Tabs defaultActiveKey="generation">
                            <TabPane tab="代数分布" key="generation">
                                <ResponsiveContainer width="100%" height={400}>
                                    <BarChart data={getGenerationDistribution()}>
                                        <CartesianGrid strokeDasharray="3 3" />
                                        <XAxis dataKey="generation" />
                                        <YAxis />
                                        <Tooltip />
                                        <Legend />
                                        <Bar dataKey="job1" fill="#8884d8" name={job1.job_name} />
                                        <Bar dataKey="job2" fill="#82ca9d" name={job2.job_name} />
                                    </BarChart>
                                </ResponsiveContainer>
                            </TabPane>

                            <TabPane tab="能量分布" key="energy">
                                <ResponsiveContainer width="100%" height={400}>
                                    <BarChart data={getEnergyDistribution()}>
                                        <CartesianGrid strokeDasharray="3 3" />
                                        <XAxis dataKey="range" />
                                        <YAxis />
                                        <Tooltip />
                                        <Legend />
                                        <Bar dataKey="job1" fill="#8884d8" name={job1.job_name} />
                                        <Bar dataKey="job2" fill="#82ca9d" name={job2.job_name} />
                                    </BarChart>
                                </ResponsiveContainer>
                            </TabPane>

                            <TabPane tab="参数对比" key="params">
                                <Row gutter={[16, 16]}>
                                    <Col span={12}>
                                        <Card title={job1.job_name} type="inner">
                                            <Statistic title="温度" value={job1.temperature} suffix="K" />
                                            <Statistic title="电压" value={job1.voltage} suffix="V" />
                                            <Statistic title="最大代数" value={job1.max_generations} />
                                            <Statistic title="能量截断" value={job1.energy_cutoff} suffix="kcal/mol" />
                                        </Card>
                                    </Col>
                                    <Col span={12}>
                                        <Card title={job2.job_name} type="inner">
                                            <Statistic title="温度" value={job2.temperature} suffix="K" />
                                            <Statistic title="电压" value={job2.voltage} suffix="V" />
                                            <Statistic title="最大代数" value={job2.max_generations} />
                                            <Statistic title="能量截断" value={job2.energy_cutoff} suffix="kcal/mol" />
                                        </Card>
                                    </Col>
                                </Row>
                            </TabPane>
                        </Tabs>
                    </>
                )}
            </Card>
        </div>
    );
};

export default ReactionNetworkComparison;
