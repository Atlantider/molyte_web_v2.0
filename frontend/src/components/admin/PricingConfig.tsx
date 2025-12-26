/**
 * 全局定价配置页面
 * 管理核时定价和任务类型定价
 */
import React, { useEffect, useState } from 'react';
import { Card, Row, Col, Typography, InputNumber, Button, Table, message, Space, Tag, Spin, theme } from 'antd';
import { DollarOutlined, SaveOutlined, ReloadOutlined, TeamOutlined, AppstoreOutlined } from '@ant-design/icons';
import {
    adminGetTaskTypePrices,
    adminUpdateTaskTypePrice,
    adminGetUserTypePrices,
    adminUpdateUserTypePrice,
    TaskTypePrice
} from '../../api/pricing';
import './PricingConfig.css';

const { Title, Text } = Typography;

// 用户类型定价接口
interface UserTypePrice {
    user_type: string;
    core_hour_price: number;
}

const PricingConfig: React.FC = () => {
    const { token } = theme.useToken();
    const [loading, setLoading] = useState(true);
    const [saving, setSaving] = useState(false);
    const [taskTypePrices, setTaskTypePrices] = useState<TaskTypePrice[]>([]);
    const [userTypePrices, setUserTypePrices] = useState<UserTypePrice[]>([]);
    const [editedTaskPrices, setEditedTaskPrices] = useState<Record<string, number>>({});
    const [editedUserPrices, setEditedUserPrices] = useState<Record<string, number>>({});

    // 任务类型名称映射
    const taskTypeNames: Record<string, string> = {
        'MD': '分子动力学模拟 (MD)',
        'QC': '量子化学计算 (QC)',
        'RESP': 'RESP电荷拟合',
        'POSTPROCESS': 'MD后处理分析',
        'CLUSTER_ANALYSIS': '团簇分析',
        'REACTION_NETWORK': '反应网络分析',
    };

    // 用户类型名称映射
    const userTypeNames: Record<string, string> = {
        'STUDENT': '学生用户',
        'RESEARCHER': '研究人员',
        'COMPANY': '企业用户',
    };

    // 加载数据
    const loadData = async () => {
        setLoading(true);
        try {
            const [taskPrices, userPrices] = await Promise.all([
                adminGetTaskTypePrices().catch(() => []),
                adminGetUserTypePrices().catch(() => [])
            ]);
            setTaskTypePrices(taskPrices);
            setUserTypePrices(userPrices);

            // 初始化编辑状态
            const taskPriceMap: Record<string, number> = {};
            taskPrices.forEach(p => { taskPriceMap[p.task_type] = p.price_per_hour; });
            setEditedTaskPrices(taskPriceMap);

            const userPriceMap: Record<string, number> = {};
            userPrices.forEach(p => { userPriceMap[p.user_type] = p.core_hour_price; });
            setEditedUserPrices(userPriceMap);
        } catch (error) {
            message.error('加载定价配置失败');
        } finally {
            setLoading(false);
        }
    };

    useEffect(() => {
        loadData();
    }, []);

    // 保存任务类型价格
    const handleSaveTaskPrice = async (taskType: string) => {
        const price = editedTaskPrices[taskType];
        if (price === undefined) return;

        setSaving(true);
        try {
            await adminUpdateTaskTypePrice(taskType, price);
            message.success(`${taskTypeNames[taskType] || taskType} 价格已更新`);
            loadData();
        } catch (error) {
            message.error('保存失败');
        } finally {
            setSaving(false);
        }
    };

    // 保存用户类型价格
    const handleSaveUserPrice = async (userType: string) => {
        const price = editedUserPrices[userType];
        if (price === undefined) return;

        setSaving(true);
        try {
            await adminUpdateUserTypePrice(userType, price);
            message.success(`${userTypeNames[userType] || userType} 价格已更新`);
            loadData();
        } catch (error) {
            message.error('保存失败');
        } finally {
            setSaving(false);
        }
    };

    // 任务类型定价表格列
    const taskColumns = [
        {
            title: '任务类型',
            dataIndex: 'task_type',
            key: 'task_type',
            render: (type: string) => (
                <Space>
                    <AppstoreOutlined style={{ color: token.colorPrimary }} />
                    <Text strong>{taskTypeNames[type] || type}</Text>
                </Space>
            ),
        },
        {
            title: '当前价格',
            dataIndex: 'price_per_hour',
            key: 'current_price',
            render: (price: number) => (
                <Tag color="blue">¥{price.toFixed(2)}/任务</Tag>
            ),
        },
        {
            title: '新价格',
            key: 'new_price',
            render: (_: any, record: TaskTypePrice) => (
                <InputNumber
                    min={0.01}
                    step={0.1}
                    precision={2}
                    value={editedTaskPrices[record.task_type]}
                    onChange={(val) => setEditedTaskPrices({ ...editedTaskPrices, [record.task_type]: val || 0 })}
                    addonBefore="¥"
                    addonAfter="/任务"
                    style={{ width: 160 }}
                />
            ),
        },
        {
            title: '操作',
            key: 'action',
            render: (_: any, record: TaskTypePrice) => (
                <Button
                    type="primary"
                    icon={<SaveOutlined />}
                    onClick={() => handleSaveTaskPrice(record.task_type)}
                    loading={saving}
                    disabled={editedTaskPrices[record.task_type] === record.price_per_hour}
                >
                    保存
                </Button>
            ),
        },
    ];

    // 用户类型定价表格列
    const userColumns = [
        {
            title: '用户类型',
            dataIndex: 'user_type',
            key: 'user_type',
            render: (type: string) => (
                <Space>
                    <TeamOutlined style={{ color: token.colorSuccess }} />
                    <Text strong>{userTypeNames[type] || type}</Text>
                </Space>
            ),
        },
        {
            title: '当前核时单价',
            dataIndex: 'core_hour_price',
            key: 'current_price',
            render: (price: number) => (
                <Tag color="green">¥{price?.toFixed(4) || '0.0000'}/核时</Tag>
            ),
        },
        {
            title: '新价格',
            key: 'new_price',
            render: (_: any, record: UserTypePrice) => (
                <InputNumber
                    min={0.0001}
                    step={0.01}
                    precision={4}
                    value={editedUserPrices[record.user_type]}
                    onChange={(val) => setEditedUserPrices({ ...editedUserPrices, [record.user_type]: val || 0 })}
                    addonBefore="¥"
                    addonAfter="/核时"
                    style={{ width: 180 }}
                />
            ),
        },
        {
            title: '操作',
            key: 'action',
            render: (_: any, record: UserTypePrice) => (
                <Button
                    type="primary"
                    icon={<SaveOutlined />}
                    onClick={() => handleSaveUserPrice(record.user_type)}
                    loading={saving}
                    disabled={editedUserPrices[record.user_type] === record.core_hour_price}
                >
                    保存
                </Button>
            ),
        },
    ];

    if (loading) {
        return (
            <div style={{ display: 'flex', justifyContent: 'center', alignItems: 'center', minHeight: 400 }}>
                <Spin size="large" tip="加载中..." />
            </div>
        );
    }

    return (
        <div style={{ padding: '24px', background: token.colorBgLayout, minHeight: 'calc(100vh - 64px)' }}>
            {/* 页面标题 */}
            <div style={{ marginBottom: 24 }}>
                <Title level={3} style={{ margin: 0, marginBottom: 4 }}>
                    <DollarOutlined style={{ marginRight: 10, color: token.colorPrimary }} />
                    全局定价配置
                </Title>
                <Text type="secondary">设置核时计费和任务计费的全局价格</Text>
            </div>

            <Row gutter={[24, 24]}>
                {/* 核时定价（按用户类型） */}
                <Col span={24}>
                    <Card
                        title={
                            <Space>
                                <TeamOutlined style={{ color: token.colorSuccess }} />
                                <span>核时定价（按用户类型）</span>
                            </Space>
                        }
                        extra={
                            <Button icon={<ReloadOutlined />} onClick={loadData}>
                                刷新
                            </Button>
                        }
                    >
                        <Text type="secondary" style={{ display: 'block', marginBottom: 16 }}>
                            设置不同用户类型的核时单价，用于"按核时计费"模式
                        </Text>
                        <Table
                            columns={userColumns}
                            dataSource={userTypePrices}
                            rowKey="user_type"
                            pagination={false}
                            size="middle"
                        />
                    </Card>
                </Col>

                {/* 任务定价 */}
                <Col span={24}>
                    <Card
                        title={
                            <Space>
                                <AppstoreOutlined style={{ color: token.colorPrimary }} />
                                <span>任务定价（按任务类型）</span>
                            </Space>
                        }
                        extra={
                            <Button icon={<ReloadOutlined />} onClick={loadData}>
                                刷新
                            </Button>
                        }
                    >
                        <Text type="secondary" style={{ display: 'block', marginBottom: 16 }}>
                            设置不同任务类型的固定价格，用于"按任务计费"模式
                        </Text>
                        <Table
                            columns={taskColumns}
                            dataSource={taskTypePrices}
                            rowKey="task_type"
                            pagination={false}
                            size="middle"
                        />
                    </Card>
                </Col>
            </Row>
        </div>
    );
};

export default PricingConfig;
