/**
 * User & Billing Management Page
 * 统一管理用户、定价、充值和消费
 */
import React, { useState, useEffect } from 'react';
import {
  Card,
  Table,
  Button,
  Tag,
  Space,
  Modal,
  Form,
  Input,
  InputNumber,
  message,
  Popconfirm,
  Row,
  Col,
  Statistic,
  Typography,
  theme,
  Tabs,
  Tooltip,
  Divider,
  Spin,
  Segmented,
} from 'antd';

const { Title, Text } = Typography;
import {
  EditOutlined,
  DeleteOutlined,
  UserOutlined,
  DollarOutlined,
  TeamOutlined,
  CrownOutlined,
  HistoryOutlined,
  BarChartOutlined,
  AppstoreOutlined,
  ExperimentOutlined,
} from '@ant-design/icons';
import AdminNav from '../../components/AdminNav';
import { useThemeStore } from '../../stores/themeStore';
import apiClient from '../../api/client';
import {
  adminGetBillingConfig,
  adminUpdateBillingConfig,
  adminGetTaskTypePrices,
  adminUpdateTaskTypePrice,
  adminGetUserTypePrices,
  adminUpdateUserTypePrice,
  type BillingConfig,
  type TaskTypePrice,
  type UserTypePrice,
} from '../../api/pricing';
import type { User } from '../../types';

interface UserBillingInfo extends User {
  custom_cpu_hour_price?: number;
  price_updated_at?: string;
  total_recharge?: number;
  total_consumption?: number;
  balance?: number;
}

interface PricingConfig {
  global_price: number;
  role_prices: Record<string, number>;
  custom_prices: Record<number, number>;
}

interface RechargeRecord {
  id: number;
  user_id: number;
  amount: number;
  cpu_hours: number;
  created_at: string;
}

interface ConsumptionRecord {
  user_id: number;
  total_consumption: number;
  total_cpu_hours: number;
  last_consumption_at: string;
}

interface DetailedConsumptionRecord {
  id: number;
  user_id: number;
  username: string;
  amount: number;
  cpu_hours: number;
  created_at: string;
}

const UserBillingManagement: React.FC = () => {
  const { token } = theme.useToken();
  useThemeStore();  // Subscribe to theme changes

  // 状态管理
  const [loading, setLoading] = useState(false);
  const [users, setUsers] = useState<UserBillingInfo[]>([]);
  const [pricingConfig, setPricingConfig] = useState<PricingConfig>({
    global_price: 0.1,
    role_prices: { ADMIN: 0.1, PREMIUM: 0.08, USER: 0.1, GUEST: 0.15 },
    custom_prices: {},
  });
  const [rechargeRecords, setRechargeRecords] = useState<RechargeRecord[]>([]);
  const [consumptionRecords, setConsumptionRecords] = useState<ConsumptionRecord[]>([]);
  const [detailedConsumption, setDetailedConsumption] = useState<DetailedConsumptionRecord[]>([]);
  const [consumptionStats, setConsumptionStats] = useState({
    total: 0,
    totalMoney: 0,
    count: 0,
    average: 0,
    averageMoney: 0,
  });

  // 新计费模型状态
  const [billingConfig, setBillingConfig] = useState<BillingConfig>({
    pricing_mode: 'CORE_HOUR',
  });
  const [taskPrices, setTaskPrices] = useState<TaskTypePrice[]>([]);
  const [userPrices, setUserPrices] = useState<UserTypePrice[]>([]);

  // 编辑状态
  const [editingUser, setEditingUser] = useState<UserBillingInfo | null>(null);
  const [editingPrice, setEditingPrice] = useState<number | null>(null);
  const [editModalVisible, setEditModalVisible] = useState(false);
  const [editingRolePrices, setEditingRolePrices] = useState<Record<string, number>>({});
  const [isEditingRolePrices, setIsEditingRolePrices] = useState(false);
  const [editingGlobalPrice, setEditingGlobalPrice] = useState<number>(0.1);
  const [isEditingGlobalPrice, setIsEditingGlobalPrice] = useState(false);
  const [form] = Form.useForm();

  // 标签页
  const [activeTab, setActiveTab] = useState('users');

  useEffect(() => {
    loadData();
  }, []);

  const loadData = async () => {
    setLoading(true);
    try {
      // 加载用户列表
      const usersRes = await apiClient.get('/admin/users');
      setUsers(usersRes.data);

      // 加载定价配置
      const pricingRes = await apiClient.get('/billing/admin/pricing-config');
      setPricingConfig(pricingRes.data);

      // 加载充值记录
      const rechargeRes = await apiClient.get('/billing/admin/recharge-records');
      setRechargeRecords(rechargeRes.data);

      // 加载消费记录
      const consumptionRes = await apiClient.get('/billing/admin/consumption-records');
      setConsumptionRecords(consumptionRes.data);

      // 加载详细消费记录
      const detailedRes = await apiClient.get('/billing/admin/consumption-details');
      setDetailedConsumption(detailedRes.data);

      // 计算消费统计
      const totalConsumption = detailedRes.data.reduce((sum: number, item: any) => sum + item.cpu_hours, 0);
      const totalMoney = detailedRes.data.reduce((sum: number, item: any) => sum + (item.money_amount || 0), 0);
      const consumptionCount = detailedRes.data.length;
      const averageConsumption = consumptionCount > 0 ? totalConsumption / consumptionCount : 0;
      const averageMoney = consumptionCount > 0 ? totalMoney / consumptionCount : 0;

      setConsumptionStats({
        total: totalConsumption,
        totalMoney: totalMoney,
        count: consumptionCount,
        average: averageConsumption,
        averageMoney: averageMoney,
      });

      // 加载新计费模型数据
      try {
        const [configRes, tasksRes, usersRes] = await Promise.all([
          adminGetBillingConfig(),
          adminGetTaskTypePrices(),
          adminGetUserTypePrices(),
        ]);
        setBillingConfig(configRes);
        setTaskPrices(tasksRes);
        setUserPrices(usersRes);
      } catch (pricingError) {
        console.error('❌ Failed to load pricing config:', pricingError);
      }
    } catch (error: any) {
      message.error(error.response?.data?.detail || '加载数据失败');
    } finally {
      setLoading(false);
    }
  };

  const handleEditPrice = (user: UserBillingInfo) => {
    setEditingUser(user);
    setEditingPrice(user.custom_cpu_hour_price || pricingConfig.global_price);
    setEditModalVisible(true);
  };

  const handleSavePrice = async () => {
    if (!editingUser || editingPrice === null) return;

    try {
      await apiClient.put(`/billing/admin/user-pricing/${editingUser.id}`, {
        user_id: editingUser.id,
        custom_cpu_hour_price: editingPrice,
      });
      message.success(`已更新 ${editingUser.username} 的定价`);
      setEditModalVisible(false);
      loadData();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '更新定价失败');
    }
  };

  const handleDeletePrice = async (userId: number, username: string) => {
    try {
      await apiClient.delete(`/billing/admin/user-pricing/${userId}`);
      message.success(`已删除 ${username} 的自定义定价`);
      loadData();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '删除定价失败');
    }
  };

  const handleEditRolePrices = () => {
    setEditingRolePrices({ ...pricingConfig.role_prices });
    setIsEditingRolePrices(true);
  };

  const handleSaveRolePrices = async () => {
    try {
      await apiClient.put('/billing/admin/role-pricing', {
        role_prices: editingRolePrices,
      });
      message.success('角色定价已更新');
      setIsEditingRolePrices(false);
      loadData();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '更新定价失败');
    }
  };

  const handleCancelEditRolePrices = () => {
    setIsEditingRolePrices(false);
    setEditingRolePrices({ ...pricingConfig.role_prices });
  };

  const handleEditGlobalPrice = () => {
    setEditingGlobalPrice(pricingConfig.global_price);
    setIsEditingGlobalPrice(true);
  };

  const handleSaveGlobalPrice = async () => {
    try {
      await apiClient.put('/billing/admin/pricing', {
        cpu_hour_price: editingGlobalPrice,
        min_recharge_amount: 10, // 保持默认值
        max_debt_cpu_hours: 100, // 保持默认值
      });
      setPricingConfig(prev => ({ ...prev, global_price: editingGlobalPrice }));
      setIsEditingGlobalPrice(false);
      message.success('全局定价更新成功');
      loadData(); // 重新加载数据
    } catch (error: any) {
      message.error(error.response?.data?.detail || '更新全局定价失败');
    }
  };

  const handleCancelEditGlobalPrice = () => {
    setIsEditingGlobalPrice(false);
    setEditingGlobalPrice(pricingConfig.global_price);
  };

  // 新计费模型处理函数
  const handleBillingModeChange = async (mode: string) => {
    try {
      await adminUpdateBillingConfig({ pricing_mode: mode as 'CORE_HOUR' | 'TASK_TYPE' });
      setBillingConfig({ ...billingConfig, pricing_mode: mode as 'CORE_HOUR' | 'TASK_TYPE' });
      message.success('计费模式已更新');
    } catch (error: any) {
      message.error('更新失败: ' + (error.message || '未知错误'));
    }
  };

  const handleTaskPriceUpdate = async (taskType: string, price: number) => {
    try {
      await adminUpdateTaskTypePrice(taskType, price);
      message.success('任务类型价格已更新');
    } catch (error: any) {
      message.error('更新失败: ' + (error.message || '未知错误'));
    }
  };

  const handleUserPriceUpdate = async (userType: string, price: number) => {
    try {
      await adminUpdateUserTypePrice(userType, price);
      message.success('用户类型核时单价已更新');
    } catch (error: any) {
      message.error('更新失败: ' + (error.message || '未知错误'));
    }
  };

  // 用户列表列定义
  const userColumns = [
    {
      title: '用户名',
      dataIndex: 'username',
      key: 'username',
      width: 120,
      render: (text: string, record: UserBillingInfo) => (
        <Space>
          <UserOutlined />
          <span>{text}</span>
        </Space>
      ),
    },
    {
      title: '邮箱',
      dataIndex: 'email',
      key: 'email',
      width: 180,
      ellipsis: true,
    },
    {
      title: '角色',
      dataIndex: 'role',
      key: 'role',
      width: 100,
      render: (role: string) => {
        const roleConfig: any = {
          ADMIN: { text: '管理员', color: 'red' },
          PREMIUM: { text: '高级用户', color: 'gold' },
          USER: { text: '普通用户', color: 'blue' },
          GUEST: { text: '访客', color: 'default' },
        };
        const config = roleConfig[role] || { text: role, color: 'default' };
        return <Tag color={config.color}>{config.text}</Tag>;
      },
    },
    {
      title: '计费模式',
      key: 'billing_mode',
      width: 120,
      render: (_: any, record: UserBillingInfo) => {
        const mode = (record as any).billing_mode || 'CORE_HOUR';
        if (mode === 'TASK_TYPE') {
          return <Tag color="purple">按任务</Tag>;
        }
        return <Tag color="blue">按核时</Tag>;
      },
    },
    {
      title: '生效价格',
      key: 'effective_price',
      width: 140,
      render: (_: any, record: UserBillingInfo) => {
        const mode = (record as any).billing_mode || 'CORE_HOUR';
        if (mode === 'TASK_TYPE') {
          const customPrices = (record as any).custom_task_prices;
          if (customPrices && Object.keys(customPrices).length > 0) {
            return <Tag color="orange">任务(自定义)</Tag>;
          }
          return <Tag color="purple">任务(标准)</Tag>;
        }
        // CORE_HOUR mode
        if (record.custom_cpu_hour_price) {
          return <Tag color="orange">¥{record.custom_cpu_hour_price}/h (自定义)</Tag>;
        }
        const rolePrice = pricingConfig.role_prices[record.role] || pricingConfig.global_price;
        return <span>¥{rolePrice}/h (标准)</span>;
      },
    },
    {
      title: '余额',
      key: 'balance',
      width: 100,
      render: (_: any, record: UserBillingInfo) => (
        <span style={{ color: token.colorSuccess }}>¥{record.balance || 0}</span>
      ),
    },
    {
      title: '操作',
      key: 'actions',
      width: 150,
      fixed: 'right' as const,
      render: (_: any, record: UserBillingInfo) => (
        <Space size={4}>
          <Tooltip title="编辑定价">
            <Button
              type="link"
              size="small"
              icon={<EditOutlined />}
              onClick={() => handleEditPrice(record)}
            />
          </Tooltip>
          {record.custom_cpu_hour_price && (
            <Popconfirm
              title="确定删除自定义定价？"
              description={`删除后，${record.username} 将恢复使用标准定价`}
              onConfirm={() => handleDeletePrice(record.id, record.username)}
              okText="确定删除"
              cancelText="取消"
            >
              <Tooltip title="删除自定义定价，恢复为标准定价">
                <Button type="link" size="small" danger icon={<DeleteOutlined />} />
              </Tooltip>
            </Popconfirm>
          )}
        </Space>
      ),
    },
  ];

  return (
    <div style={{ padding: '20px 24px', background: token.colorBgLayout, minHeight: 'calc(100vh - 64px)' }}>
      {/* 页面标题 */}
      <div style={{ marginBottom: 16 }}>
        <Title level={3} style={{ margin: 0, marginBottom: 4 }}>
          <DollarOutlined style={{ marginRight: 10, color: token.colorPrimary }} />
          计费管理
        </Title>
        <Text type="secondary">管理用户定价、充值和消费统计</Text>
      </div>

      <AdminNav />

      {/* 统计卡片 - 与用户管理页面统一风格 */}
      <Row gutter={16} style={{ marginBottom: 20 }}>
        {[
          { label: '总用户数', value: users.length, color: '#85a5ff', icon: <TeamOutlined /> },
          { label: '全局单价', value: `¥${pricingConfig.global_price}/核时`, color: '#95de64', icon: <DollarOutlined />, isText: true },
          { label: '自定义定价', value: Object.keys(pricingConfig.custom_prices).length, color: '#faad14', icon: <CrownOutlined /> },
          { label: '总消费核时', value: `${consumptionStats.total.toFixed(1)}h`, color: '#ff4d4f', icon: <HistoryOutlined />, isText: true },
          { label: '总充值核时', value: `${rechargeRecords.reduce((sum, r) => sum + r.cpu_hours, 0).toFixed(1)}h`, color: '#52c41a', icon: <HistoryOutlined />, isText: true },
          { label: '消费用户', value: consumptionRecords.length, color: '#b37feb', icon: <UserOutlined /> },
        ].map((item, idx) => (
          <Col xs={12} sm={8} md={4} key={idx}>
            <div style={{
              padding: '16px 20px',
              background: token.colorBgContainer,
              borderRadius: 10,
              border: `1px solid ${token.colorBorder}`,
              borderLeft: `4px solid ${item.color}`,
              display: 'flex',
              alignItems: 'center',
              gap: 14,
            }}>
              <div style={{
                width: 42,
                height: 42,
                borderRadius: 10,
                background: `${item.color}15`,
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                fontSize: 20,
                color: item.color,
              }}>
                {item.icon}
              </div>
              <div>
                <Text type="secondary" style={{ fontSize: 13, display: 'block' }}>{item.label}</Text>
                <Text strong style={{ fontSize: item.isText ? 18 : 22, color: item.color }}>{item.value}</Text>
              </div>
            </div>
          </Col>
        ))}
      </Row>

      {/* 标签页 */}
      <Card bordered={false} style={{ borderRadius: 12 }}>
        <Spin spinning={loading}>
          <Tabs
            activeKey={activeTab}
            onChange={setActiveTab}
            items={[
              {
                key: 'users',
                label: <span><UserOutlined /> 用户列表</span>,
                children: (
                  <Table
                    columns={userColumns}
                    dataSource={users}
                    rowKey="id"
                    pagination={{ pageSize: 20 }}
                    scroll={{ x: 1200 }}
                  />
                ),
              },
              {
                key: 'pricing',
                label: <span><DollarOutlined /> 定价管理</span>,
                children: (() => {
                  const taskNames: Record<string, { name: string; icon: string; color: string }> = {
                    FORCEFIELD: { name: '力场生成', icon: '🧲', color: '#722ed1' },
                    MD: { name: 'MD计算', icon: '⚛️', color: '#13c2c2' },
                    POSTPROCESS: { name: '后处理', icon: '📊', color: '#52c41a' },
                    QC: { name: 'QC计算', icon: '🔬', color: '#1890ff' },
                    REACTION_NETWORK: { name: '反应网络', icon: '🔗', color: '#fa8c16' },
                  };

                  const userNames: Record<string, { name: string; icon: string; color: string }> = {
                    ADMIN: { name: '管理员', icon: '👑', color: '#722ed1' },
                    GUEST: { name: '访客', icon: '👤', color: '#8c8c8c' },
                    USER: { name: '普通用户', icon: '👥', color: '#1890ff' },
                    PREMIUM: { name: '高级用户', icon: '⭐', color: '#faad14' },
                  };

                  return (
                    <div style={{ padding: '20px 0' }}>
                      {/* 计费模式选择 */}
                      <Card bordered={false} style={{ marginBottom: 24 }}>
                        <Space direction="vertical" size="large" style={{ width: '100%' }}>
                          <div>
                            <Text strong style={{ display: 'block', marginBottom: 12, fontSize: 16 }}>
                              计费模式
                            </Text>
                            <Segmented
                              value={billingConfig.pricing_mode}
                              onChange={handleBillingModeChange}
                              block
                              size="large"
                              options={[
                                {
                                  label: (
                                    <div style={{ padding: '8px 24px', textAlign: 'center' }}>
                                      <DollarOutlined style={{ fontSize: 20, display: 'block', marginBottom: 8 }} />
                                      <div>按核时计费</div>
                                      <Text type="secondary" style={{ fontSize: 12 }}>
                                        不同用户类型不同单价
                                      </Text>
                                    </div>
                                  ),
                                  value: 'CORE_HOUR',
                                },
                                {
                                  label: (
                                    <div style={{ padding: '8px 24px', textAlign: 'center' }}>
                                      <AppstoreOutlined style={{ fontSize: 20, display: 'block', marginBottom: 8 }} />
                                      <div>按任务计费</div>
                                      <Text type="secondary" style={{ fontSize: 12 }}>
                                        每种任务固定价格
                                      </Text>
                                    </div>
                                  ),
                                  value: 'TASK_TYPE',
                                },
                              ]}
                            />
                          </div>
                        </Space>
                      </Card>

                      {/* 价格配置 */}
                      {billingConfig.pricing_mode === 'CORE_HOUR' ? (
                        <Card
                          title={
                            <Space>
                              <TeamOutlined />
                              <span>用户类型核时单价</span>
                            </Space>
                          }
                          bordered={false}
                        >
                          <Row gutter={[16, 16]}>
                            {userPrices.map((item) => {
                              const config = userNames[item.user_type];
                              return (
                                <Col xs={24} sm={12} lg={6} key={item.user_type}>
                                  <Card
                                    size="small"
                                    style={{
                                      borderRadius: 8,
                                      borderLeft: `4px solid ${config?.color || '#1890ff'}`,
                                    }}
                                  >
                                    <Statistic
                                      title={
                                        <Space>
                                          <span style={{ fontSize: 20 }}>{config?.icon}</span>
                                          <span>{config?.name || item.user_type}</span>
                                        </Space>
                                      }
                                      value={item.core_hour_price}
                                      precision={2}
                                      prefix="¥"
                                      suffix="/核时"
                                      valueStyle={{ fontSize: 24, fontWeight: 600 }}
                                    />
                                    <div style={{ marginTop: 12 }}>
                                      <InputNumber
                                        value={item.core_hour_price}
                                        onChange={(v) => {
                                          if (v !== null && v >= 0) {
                                            setUserPrices((prev) =>
                                              prev.map((p) =>
                                                p.user_type === item.user_type
                                                  ? { ...p, core_hour_price: v }
                                                  : p
                                              )
                                            );
                                          }
                                        }}
                                        onBlur={() => handleUserPriceUpdate(item.user_type, item.core_hour_price)}
                                        prefix="¥"
                                        suffix="/核时"
                                        min={0}
                                        step={0.1}
                                        precision={2}
                                        style={{ width: '100%' }}
                                      />
                                    </div>
                                  </Card>
                                </Col>
                              );
                            })}
                          </Row>
                        </Card>
                      ) : (
                        <Card
                          title={
                            <Space>
                              <ExperimentOutlined />
                              <span>任务类型价格</span>
                            </Space>
                          }
                          bordered={false}
                        >
                          <Row gutter={[16, 16]}>
                            {taskPrices.map((item) => {
                              const config = taskNames[item.task_type];
                              return (
                                <Col xs={24} sm={12} md={8} lg={6} key={item.task_type}>
                                  <Card
                                    size="small"
                                    style={{
                                      borderRadius: 8,
                                      borderLeft: `4px solid ${config?.color || '#1890ff'}`,
                                    }}
                                  >
                                    <Statistic
                                      title={
                                        <Space>
                                          <span style={{ fontSize: 20 }}>{config?.icon}</span>
                                          <span>{config?.name || item.task_type}</span>
                                        </Space>
                                      }
                                      value={item.price_per_hour}
                                      precision={2}
                                      prefix="¥"
                                      suffix="/任务"
                                      valueStyle={{ fontSize: 24, fontWeight: 600 }}
                                    />
                                    <div style={{ marginTop: 12 }}>
                                      <InputNumber
                                        value={item.price_per_hour}
                                        onChange={(v) => {
                                          if (v !== null && v > 0) {
                                            setTaskPrices((prev) =>
                                              prev.map((p) =>
                                                p.task_type === item.task_type
                                                  ? { ...p, price_per_hour: v }
                                                  : p
                                              )
                                            );
                                          }
                                        }}
                                        onBlur={() => handleTaskPriceUpdate(item.task_type, item.price_per_hour)}
                                        prefix="¥"
                                        suffix="/任务"
                                        min={0.01}
                                        step={1}
                                        precision={2}
                                        style={{ width: '100%' }}
                                      />
                                    </div>
                                  </Card>
                                </Col>
                              );
                            })}
                          </Row>
                        </Card>
                      )}
                    </div>
                  );
                })(),
              },
              {
                key: 'recharge',
                label: <span><HistoryOutlined /> 充值记录</span>,
                children: (
                  <Table
                    columns={[
                      {
                        title: '用户ID',
                        dataIndex: 'user_id',
                        key: 'user_id',
                        width: 80,
                        sorter: (a, b) => a.user_id - b.user_id,
                      },
                      {
                        title: '用户名',
                        key: 'username',
                        width: 120,
                        render: (_, record) => {
                          const user = users.find(u => u.id === record.user_id);
                          return user?.username || '-';
                        }
                      },
                      {
                        title: '充值金额',
                        dataIndex: 'amount',
                        key: 'amount',
                        width: 100,
                        render: (val) => `¥${val.toFixed(2)}`,
                        sorter: (a, b) => a.amount - b.amount,
                      },
                      {
                        title: '获得核时',
                        dataIndex: 'cpu_hours',
                        key: 'cpu_hours',
                        width: 100,
                        render: (val) => `${val.toFixed(2)}h`,
                        sorter: (a, b) => a.cpu_hours - b.cpu_hours,
                      },
                      {
                        title: '充值时间',
                        dataIndex: 'created_at',
                        key: 'created_at',
                        width: 180,
                        render: (val) => new Date(val).toLocaleString('zh-CN'),
                        sorter: (a, b) => new Date(a.created_at).getTime() - new Date(b.created_at).getTime(),
                      },
                    ]}
                    dataSource={rechargeRecords}
                    rowKey="id"
                    pagination={{ pageSize: 20 }}
                    scroll={{ x: 800 }}
                  />
                ),
              },
              {
                key: 'consumption',
                label: <span><BarChartOutlined /> 消费统计</span>,
                children: (
                  <Tabs
                    items={[
                      {
                        key: 'summary',
                        label: '消费汇总',
                        children: (
                          <Table
                            columns={[
                              {
                                title: '用户ID',
                                dataIndex: 'user_id',
                                key: 'user_id',
                                width: 80,
                                sorter: (a, b) => a.user_id - b.user_id,
                              },
                              {
                                title: '用户名',
                                key: 'username',
                                width: 120,
                                render: (_, record) => {
                                  const user = users.find(u => u.id === record.user_id);
                                  return user?.username || '-';
                                }
                              },
                              {
                                title: '总消费金额',
                                dataIndex: 'total_consumption',
                                key: 'total_consumption',
                                width: 120,
                                render: (val) => `¥${val.toFixed(2)}`,
                                sorter: (a, b) => a.total_consumption - b.total_consumption,
                              },
                              {
                                title: '消耗核时',
                                dataIndex: 'total_cpu_hours',
                                key: 'total_cpu_hours',
                                width: 100,
                                render: (val) => `${val.toFixed(2)}h`,
                                sorter: (a, b) => a.total_cpu_hours - b.total_cpu_hours,
                              },
                              {
                                title: '最后消费时间',
                                dataIndex: 'last_consumption_at',
                                key: 'last_consumption_at',
                                width: 180,
                                render: (val) => val ? new Date(val).toLocaleString('zh-CN') : '-',
                                sorter: (a, b) => new Date(a.last_consumption_at).getTime() - new Date(b.last_consumption_at).getTime(),
                              },
                            ]}
                            dataSource={consumptionRecords}
                            rowKey="user_id"
                            pagination={{ pageSize: 20 }}
                            scroll={{ x: 800 }}
                          />
                        ),
                      },
                      {
                        key: 'details',
                        label: '消费详情',
                        children: (
                          <Table
                            columns={[
                              {
                                title: '用户ID',
                                dataIndex: 'user_id',
                                key: 'user_id',
                                width: 80,
                                sorter: (a, b) => a.user_id - b.user_id,
                              },
                              {
                                title: '用户名',
                                dataIndex: 'username',
                                key: 'username',
                                width: 120,
                              },
                              {
                                title: '消费核时',
                                dataIndex: 'cpu_hours',
                                key: 'cpu_hours',
                                width: 100,
                                render: (val) => `${val.toFixed(2)}h`,
                                sorter: (a, b) => a.cpu_hours - b.cpu_hours,
                              },
                              {
                                title: '消费金额',
                                dataIndex: 'amount',
                                key: 'amount',
                                width: 100,
                                render: (val) => `¥${val.toFixed(2)}`,
                                sorter: (a, b) => a.amount - b.amount,
                              },
                              {
                                title: '消费时间',
                                dataIndex: 'created_at',
                                key: 'created_at',
                                width: 180,
                                render: (val) => new Date(val).toLocaleString('zh-CN'),
                                sorter: (a, b) => new Date(a.created_at).getTime() - new Date(b.created_at).getTime(),
                              },
                            ]}
                            dataSource={detailedConsumption}
                            rowKey="id"
                            pagination={{ pageSize: 20 }}
                            scroll={{ x: 800 }}
                            loading={loading}
                          />
                        ),
                      },
                    ]}
                  />
                ),
              },
            ]}
          />
        </Spin>
      </Card>

      {/* 编辑定价模态框 */}
      <Modal
        title={`编辑 ${editingUser?.username} 的定价`}
        open={editModalVisible}
        onOk={handleSavePrice}
        onCancel={() => setEditModalVisible(false)}
        okText="保存"
        cancelText="取消"
        width={600}
      >
        <Form layout="vertical">
          <Form.Item label="用户名">
            <Input value={editingUser?.username} disabled />
          </Form.Item>
          <Form.Item label="邮箱">
            <Input value={editingUser?.email} disabled />
          </Form.Item>
          <Form.Item label="角色">
            <Input value={editingUser?.role} disabled />
          </Form.Item>
          <Form.Item label="全局价格">
            <Input value={`¥${pricingConfig.role_prices[editingUser?.role || 'USER'] || pricingConfig.global_price}`} disabled />
          </Form.Item>
          <Divider />
          <Form.Item
            label="自定义核时单价 (元/核时)"
            tooltip="设置此用户的自定义价格，将覆盖全局定价。新价格将在下一次充值时生效。"
          >
            <InputNumber
              value={editingPrice}
              onChange={(val) => setEditingPrice(val)}
              min={0.001}
              step={0.01}
              precision={4}
              style={{ width: '100%' }}
              placeholder="输入自定义价格，留空则使用全局定价"
            />
          </Form.Item>
          <Text type="secondary" style={{ fontSize: 12 }}>
            💡 提示：修改价格后，新价格将在用户下一次充值时生效。已消费的任务不会重新计费。
          </Text>
        </Form>
      </Modal>
    </div>
  );
};

export default UserBillingManagement;

