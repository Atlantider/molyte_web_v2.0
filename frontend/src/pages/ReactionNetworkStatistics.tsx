/**
 * Reaction Network Statistics Dashboard
 * 反应网络统计仪表板
 */

import React, { useState, useEffect } from 'react';
import {
    Card,
    Row,
    Col,
    Statistic,
    Progress,
    Table,
    Tag,
    Space,
    Spin
} from 'antd';
import {
    BarChart,
    Bar,
    PieChart,
    Pie,
    Cell,
    LineChart,
    Line,
    XAxis,
    YAxis,
    CartesianGrid,
    Tooltip,
    Legend,
    ResponsiveContainer
} from 'recharts';
import {
    RiseOutlined,
    FallOutlined,
    ExperimentOutlined,
    FireOutlined
} from '@ant-design/icons';
import { useParams } from 'react-router-dom';

const ReactionNetworkStatistics: React.FC = () => {
    const { id } = useParams<{ id: string }>();
    const [loading, setLoading] = useState(true);
    const [statistics, setStatistics] = useState<any>(null);
    const [energyAnalysis, setEnergyAnalysis] = useState<any>(null);

    useEffect(() => {
        loadStatistics();
        loadEnergyAnalysis();
    }, [id]);

    const loadStatistics = async () => {
        try {
            // TODO: 从API加载统计数据
            // const data = await getJobStatistics(Number(id));
            // setStatistics(data);
        } catch (error) {
            console.error('Failed to load statistics', error);
        } finally {
            setLoading(false);
        }
    };

    const loadEnergyAnalysis = async () => {
        try {
            // TODO: 从API加载能量分析
            // const data = await getEnergyAnalysis(Number(id));
            // setEnergyAnalysis(data);
        } catch (error) {
            console.error('Failed to load energy analysis', error);
        }
    };

    const COLORS = ['#0088FE', '#00C49F', '#FFBB28', '#FF8042', '#8884d8', '#82ca9d'];

    // 模拟数据
    const mockGenerationData = [
        { generation: 'Gen 0', count: 3 },
        { generation: 'Gen 1', count: 12 },
        { generation: 'Gen 2', count: 25 },
        { generation: 'Gen 3', count: 18 }
    ];

    const mockEnergyData = [
        { name: '< -50 kcal/mol', value: 5 },
        { name: '-50~0 kcal/mol', value: 35 },
        { name: '0~50 kcal/mol', value: 48 },
        { name: '> 50 kcal/mol', value: 12 }
    ];

    const mockOperatorData = [
        { operator: '环开裂', count: 15 },
        { operator: 'H转移', count: 22 },
        { operator: '电子注入', count: 18 },
        { operator: '自由基加成', count: 12 },
        { operator: '取代', count: 8 }
    ];

    const topReactionsColumns = [
        { title: '反应ID', dataIndex: 'id', key: 'id', width: 80 },
        {
            title: '反应物',
            dataIndex: 'reactants',
            key: 'reactants',
            render: (smiles: string[]) => smiles.map((s, i) => <Tag key={i}>{s.substring(0, 20)}...</Tag>)
        },
        {
            title: '产物',
            dataIndex: 'products',
            key: 'products',
            render: (smiles: string[]) => smiles.map((s, i) => <Tag key={i} color="green">{s.substring(0, 20)}...</Tag>)
        },
        {
            title: '算符',
            dataIndex: 'operator',
            key: 'operator',
            render: (op: string) => <Tag color="orange">{op}</Tag>
        },
        {
            title: '反应能',
            dataIndex: 'energy',
            key: 'energy',
            render: (energy: number) => (
                <span style={{ fontWeight: 'bold', color: energy < 0 ? '#52c41a' : '#f5222d' }}>
                    {energy.toFixed(2)} kcal/mol
                </span>
            ),
            sorter: (a: any, b: any) => a.energy - b.energy
        }
    ];

    if (loading) {
        return <div style={{ textAlign: 'center', padding: '100px' }}><Spin size="large" /></div>;
    }

    return (
        <div style={{ padding: '24px' }}>
            {/* 概览统计 */}
            <Row gutter={16} style={{ marginBottom: '24px' }}>
                <Col span={6}>
                    <Card>
                        <Statistic
                            title="总分子数"
                            value={58}
                            prefix={<ExperimentOutlined />}
                            valueStyle={{ color: '#1890ff' }}
                        />
                    </Card>
                </Col>
                <Col span={6}>
                    <Card>
                        <Statistic
                            title="总反应数"
                            value={100}
                            prefix={<FireOutlined />}
                            valueStyle={{ color: '#faad14' }}
                        />
                    </Card>
                </Col>
                <Col span={6}>
                    <Card>
                        <Statistic
                            title="放能反应"
                            value={35}
                            suffix="/ 100"
                            prefix={<FallOutlined />}
                            valueStyle={{ color: '#52c41a' }}
                        />
                        <Progress percent={35} strokeColor="#52c41a" showInfo={false} />
                    </Card>
                </Col>
                <Col span={6}>
                    <Card>
                        <Statistic
                            title="吸能反应"
                            value={65}
                            suffix="/ 100"
                            prefix={<RiseOutlined />}
                            valueStyle={{ color: '#f5222d' }}
                        />
                        <Progress percent={65} strokeColor="#f5222d" showInfo={false} />
                    </Card>
                </Col>
            </Row>

            {/* 图表区域 */}
            <Row gutter={16} style={{ marginBottom: '24px' }}>
                <Col span={12}>
                    <Card title="📊 代数分布">
                        <ResponsiveContainer width="100%" height={300}>
                            <BarChart data={mockGenerationData}>
                                <CartesianGrid strokeDasharray="3 3" />
                                <XAxis dataKey="generation" />
                                <YAxis />
                                <Tooltip />
                                <Bar dataKey="count" fill="#1890ff" />
                            </BarChart>
                        </ResponsiveContainer>
                    </Card>
                </Col>
                <Col span={12}>
                    <Card title="🎯 能量分布">
                        <ResponsiveContainer width="100%" height={300}>
                            <PieChart>
                                <Pie
                                    data={mockEnergyData}
                                    cx="50%"
                                    cy="50%"
                                    labelLine={false}
                                    label={(entry) => `${entry.name}: ${entry.value}`}
                                    outerRadius={100}
                                    fill="#8884d8"
                                    dataKey="value"
                                >
                                    {mockEnergyData.map((entry, index) => (
                                        <Cell key={`cell-${index}`} fill={COLORS[index % COLORS.length]} />
                                    ))}
                                </Pie>
                                <Tooltip />
                            </PieChart>
                        </ResponsiveContainer>
                    </Card>
                </Col>
            </Row>

            <Row gutter={16} style={{ marginBottom: '24px' }}>
                <Col span={24}>
                    <Card title="🔧 算符使用统计">
                        <ResponsiveContainer width="100%" height={300}>
                            <BarChart data={mockOperatorData} layout="horizontal">
                                <CartesianGrid strokeDasharray="3 3" />
                                <XAxis dataKey="operator" />
                                <YAxis />
                                <Tooltip />
                                <Bar dataKey="count" fill="#52c41a" />
                            </BarChart>
                        </ResponsiveContainer>
                    </Card>
                </Col>
            </Row>

            {/* 关键反应 */}
            <Row gutter={16}>
                <Col span={24}>
                    <Card title="⚡ 能量最低的反应（Top 10）">
                        <Table
                            columns={topReactionsColumns}
                            dataSource={[]} // TODO: 从API加载
                            pagination={false}
                            size="small"
                        />
                    </Card>
                </Col>
            </Row>
        </div>
    );
};

export default ReactionNetworkStatistics;
