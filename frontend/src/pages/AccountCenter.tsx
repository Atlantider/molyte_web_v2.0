/**
 * 统一的账户中心页面
 * 
 * 使用 Tab 模式整合所有账户相关功能：
 * - 账户总览：显示账户信息、配额状态、快速操作
 * - 充值中心：根据账户类型显示不同的充值选项
 * - 账户设置：个人信息、密码、偏好设置
 * - 使用统计：配额使用历史、任务统计
 */

import React, { useState, useEffect } from 'react';
import {
  Card,
  Tabs,
  Row,
  Col,
  Statistic,
  Button,
  Space,
  Alert,
  Progress,
  Tag,
  Avatar,
  Typography,
  Divider,
  Form,
  Input,
  InputNumber,
  Select,
  message,
  Spin,
  Table,
  DatePicker,
  Radio,
  Modal,
  theme
} from 'antd';
import {
  UserOutlined,
  WalletOutlined,
  SettingOutlined,
  BarChartOutlined,
  CrownOutlined,
  TeamOutlined,
  BankOutlined,
  ReloadOutlined,
  EditOutlined,
  LockOutlined,
  HistoryOutlined,
  DollarOutlined,
  GiftOutlined,
  RocketOutlined,
  CheckCircleOutlined,
  ThunderboltOutlined,
  ClusterOutlined,
  LineChartOutlined,
  DatabaseOutlined,
} from '@ant-design/icons';
import { useNavigate, useSearchParams } from 'react-router-dom';
import { useThemeStore } from '../stores/themeStore';
import { useAuthStore } from '../stores/authStore';
import { useAccountType } from '../hooks/useAccountType';
import { useQuota } from '../hooks/useQuota';
import SubAccountManagementPanel from '../components/SubAccountManagementPanel';
import { getBalance, createOrder, BalanceInfo, RechargeOrder } from '../api/billing';
import { formatCpuHours, formatBalance, QUOTA_PRECISION, safeToFixed } from '../utils/formatQuotaDisplay';
import { getMySubAccountInfo } from '../api/accounts';
// 主账号管理和用户定价已移至管理面板
// TODO: 创建新的账户管理API文件
// import { getMyOrganizationInfo } from '../api/accounts';

const { Title, Text, Paragraph } = Typography;
const { TabPane } = Tabs;
const { Option } = Select;

// 预设充值金额
const PRESET_AMOUNTS = [10, 50, 100, 200, 500, 1000];

const AccountCenter: React.FC = () => {
  const navigate = useNavigate();
  const [searchParams] = useSearchParams();
  const { mode } = useThemeStore();
  const { user } = useAuthStore();
  const { token } = theme.useToken();
  const isDark = mode === 'dark';

  const { accountInfo, loading: accountLoading, refetch: refetchAccount } = useAccountType();
  const { quota, loading: quotaLoading, refetch: refetchQuota } = useQuota();

  // 从 URL 参数中获取初始 tab，如果没有则默认为 'overview'
  const initialTab = searchParams.get('tab') || 'overview';
  const [activeTab, setActiveTab] = useState(initialTab);
  const [form] = Form.useForm();
  const [passwordForm] = Form.useForm();
  const [editOrgForm] = Form.useForm();

  // 使用统计相关状态
  const [statsLoading, setStatsLoading] = useState(false);
  const [dailyStats, setDailyStats] = useState<any[]>([]);
  const [statsPeriod, setStatsPeriod] = useState(7);

  // 充值相关状态
  const [balance, setBalance] = useState<BalanceInfo | null>(null);
  const [rechargeAmount, setRechargeAmount] = useState<number>(100);
  const [paymentMethod, setPaymentMethod] = useState<string>('simulated');
  const [rechargeLoading, setRechargeLoading] = useState(false);
  const [rechargeModalVisible, setRechargeModalVisible] = useState(false);

  // 子账号信息
  const [subAccountInfo, setSubAccountInfo] = useState<any>(null);
  const [subAccountLoading, setSubAccountLoading] = useState(false);

  const loading = accountLoading || quotaLoading;

  // 监听 URL 参数变化，自动切换 tab
  useEffect(() => {
    const tabParam = searchParams.get('tab');
    if (tabParam && tabParam !== activeTab) {
      setActiveTab(tabParam);
    }
  }, [searchParams]);

  // 获取账户类型配置
  const getAccountConfig = () => {
    if (!accountInfo) return { label: '未知', color: 'default', icon: <UserOutlined /> };

    const { account_type } = accountInfo;
    const configs = {
      'PERSONAL': { label: '个人用户', color: 'blue', icon: <UserOutlined /> },
      'personal': { label: '个人用户', color: 'blue', icon: <UserOutlined /> },
      'MASTER_ACCOUNT': { label: '主账号', color: 'gold', icon: <CrownOutlined /> },
      'master_account': { label: '主账号', color: 'gold', icon: <CrownOutlined /> },
      'SUB_ACCOUNT': { label: '子账号', color: 'purple', icon: <BankOutlined /> },
      'sub_account': { label: '子账号', color: 'purple', icon: <BankOutlined /> }
    };
    return configs[account_type as keyof typeof configs] || configs['personal'];
  };

  // 获取配额状态
  const getQuotaStatus = () => {
    if (!quota) return { status: 'normal', color: '#52c41a' };
    
    const available = quota.available_quota;
    if (available < 1) return { status: 'critical', color: '#ff4d4f' };
    if (available < 10) return { status: 'warning', color: '#faad14' };
    return { status: 'normal', color: '#52c41a' };
  };

  // 刷新数据
  const handleRefresh = () => {
    refetchAccount();
    refetchQuota();
  };

  // 加载余额信息
  const loadBalance = async () => {
    try {
      const balanceData = await getBalance();
      setBalance(balanceData);
    } catch (error) {
      console.error('Failed to load balance:', error);
    }
  };

  // 加载子账号信息（包括主账号信息）
  const loadSubAccountInfo = async () => {
    if (!accountInfo) {
      console.warn('accountInfo is not available yet');
      return;
    }

    const isSubAccount = accountInfo.account_type === 'SUB_ACCOUNT' || accountInfo.account_type === 'sub_account';
    console.log('loadSubAccountInfo - accountInfo:', accountInfo, 'isSubAccount:', isSubAccount);

    if (!isSubAccount) {
      return;
    }

    setSubAccountLoading(true);
    try {
      console.log('Fetching sub-account info...');
      const info = await getMySubAccountInfo();
      console.log('Sub-account info loaded successfully:', info);
      setSubAccountInfo(info);
    } catch (error: any) {
      console.error('Failed to load sub-account info:', error);
      message.error(error.response?.data?.detail || '加载子账号信息失败');
    } finally {
      setSubAccountLoading(false);
    }
  };

  // 处理充值
  const handleRecharge = async () => {
    if (!rechargeAmount || rechargeAmount < 10) {
      message.error('充值金额必须大于等于 10 元');
      return;
    }

    setRechargeLoading(true);
    try {
      const order = await createOrder({
        amount: rechargeAmount,
        payment_method: paymentMethod
      });

      message.success(`订单创建成功，订单号：${order.order_no}`);
      setRechargeModalVisible(false);
      setRechargeAmount(100);

      // 刷新余额和配额
      await loadBalance();
      await refetchQuota();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '充值失败，请重试');
    } finally {
      setRechargeLoading(false);
    }
  };

  // 处理组织配额申请
  const handleOrganizationQuotaRequest = () => {
    message.info('组织配额申请功能已提交，请联系管理员审核。');
    // TODO: 实现实际的组织配额申请逻辑
    // 可以打开一个模态框让用户填写申请信息
  };

  // 处理个人信息更新
  const handleUpdateProfile = async (values: any) => {
    try {
      // TODO: 调用更新个人信息的API
      message.success('个人信息更新成功');
    } catch (error) {
      message.error('更新失败，请重试');
    }
  };

  // 处理密码修改
  const handleChangePassword = async (values: any) => {
    if (values.newPassword !== values.confirmPassword) {
      message.error('两次输入的密码不一致');
      return;
    }

    try {
      // TODO: 调用修改密码的API
      message.success('密码修改成功');
      passwordForm.resetFields();
    } catch (error) {
      message.error('密码修改失败，请重试');
    }
  };

  // 加载每日统计数据
  const loadDailyStats = async () => {
    setStatsLoading(true);
    try {
      const token = localStorage.getItem('access_token');
      console.log('Loading daily stats with token:', token ? 'present' : 'missing', 'days:', statsPeriod);

      const response = await fetch(`/api/v1/users/me/daily-stats?days=${statsPeriod}`, {
        headers: {
          'Authorization': `Bearer ${token}`
        }
      });

      console.log('Stats API response status:', response.status);

      if (response.ok) {
        const data = await response.json();
        console.log('Daily stats loaded successfully:', data);
        setDailyStats(Array.isArray(data) ? data : []);
      } else {
        const errorData = await response.text();
        console.warn('Failed to load stats:', response.status, errorData);
        setDailyStats([]);
        if (response.status === 401) {
          message.error('身份验证失败，请重新登录');
        }
      }
    } catch (error) {
      console.error('Error loading daily stats:', error);
      message.error('加载统计数据失败');
      setDailyStats([]);
    } finally {
      setStatsLoading(false);
    }
  };

  // 当统计周期改变时重新加载数据
  React.useEffect(() => {
    if (activeTab === 'statistics') {
      loadDailyStats();
    }
  }, [statsPeriod, activeTab]);

  // 当账户信息改变时加载子账号信息
  React.useEffect(() => {
    if (accountInfo && (accountInfo.account_type === 'SUB_ACCOUNT' || accountInfo.account_type === 'sub_account')) {
      loadSubAccountInfo();
    }
  }, [accountInfo]);

  // 当切换到账号管理 tab 时也加载子账号信息
  React.useEffect(() => {
    if (activeTab === 'account-management') {
      if (accountInfo && (accountInfo.account_type === 'SUB_ACCOUNT' || accountInfo.account_type === 'sub_account')) {
        loadSubAccountInfo();
      }
    }
  }, [activeTab]);

  // 当切换到充值 Tab 时加载余额信息
  React.useEffect(() => {
    if (activeTab === 'recharge') {
      loadBalance();
    }
  }, [activeTab]);

  // 账户总览 Tab
  const renderOverviewTab = () => {
    const accountConfig = getAccountConfig();
    const quotaStatus = getQuotaStatus();

    return (
      <div>
        {/* 配额警告 */}
        {quota && quota.available_quota < 10 && (
          <Alert
            message={quota.available_quota < 1 ? "配额即将用尽" : "配额不足警告"}
            description={`您的可用配额仅剩 ${formatCpuHours(quota.available_quota)} 核时，建议及时充值。`}
            type={quota.available_quota < 1 ? "error" : "warning"}
            showIcon
            style={{ marginBottom: 24 }}
            action={
              <Button 
                size="small" 
                type="primary"
                onClick={() => setActiveTab('recharge')}
              >
                立即充值
              </Button>
            }
          />
        )}

        <Row gutter={[24, 24]}>
          {/* 账户信息卡片 */}
          <Col xs={24} lg={8}>
            <Card title="账户信息" extra={<Button icon={<EditOutlined />} size="small">编辑</Button>}>
              <div style={{ textAlign: 'center', marginBottom: 16 }}>
                <Avatar size={64} icon={<UserOutlined />} style={{ marginBottom: 8 }} />
                <div>
                  <Title level={4} style={{ margin: 0 }}>{user?.username}</Title>
                  <Text type="secondary">{user?.email}</Text>
                </div>
              </div>
              
              <Divider />
              
              <div style={{ marginBottom: 16 }}>
                <Tag 
                  icon={accountConfig.icon} 
                  color={accountConfig.color}
                  style={{ width: '100%', textAlign: 'center', padding: '8px 0' }}
                >
                  {accountConfig.label}
                </Tag>
              </div>

              <Space direction="vertical" style={{ width: '100%' }}>
                <Button 
                  block 
                  icon={<SettingOutlined />}
                  onClick={() => setActiveTab('settings')}
                >
                  账户设置
                </Button>
                <Button 
                  block 
                  icon={<HistoryOutlined />}
                  onClick={() => setActiveTab('statistics')}
                >
                  使用统计
                </Button>
              </Space>
            </Card>
          </Col>

          {/* 配额状态卡片 */}
          <Col xs={24} lg={16}>
            <Card 
              title="配额状态" 
              extra={
                <Space>
                  <Button icon={<ReloadOutlined />} onClick={handleRefresh}>刷新</Button>
                  <Button 
                    type="primary" 
                    icon={<WalletOutlined />}
                    onClick={() => setActiveTab('recharge')}
                  >
                    充值配额
                  </Button>
                </Space>
              }
            >
              <Row gutter={16}>
                <Col xs={24} sm={8}>
                  <Statistic
                    title="可用余额"
                    value={quota?.balance_cpu_hours || 0}
                    suffix="核时"
                    valueStyle={{ color: quota?.balance_cpu_hours && quota.balance_cpu_hours < 0 ? '#ff4d4f' : '#52c41a' }}
                    precision={QUOTA_PRECISION}
                  />
                </Col>
                <Col xs={24} sm={8}>
                  <Statistic
                    title="冻结核时"
                    value={quota?.frozen_cpu_hours || 0}
                    suffix="核时"
                    valueStyle={{ color: '#faad14' }}
                    precision={QUOTA_PRECISION}
                  />
                </Col>
                <Col xs={24} sm={8}>
                  <Statistic
                    title="已消费"
                    value={quota?.used_cpu_hours || 0}
                    suffix="核时"
                    valueStyle={{ color: '#1890ff' }}
                    precision={QUOTA_PRECISION}
                  />
                </Col>
              </Row>

              {/* 核时来源 */}
              {quota?.quota_sources && (
                <div style={{ marginTop: 24 }}>
                  <Text strong>核时来源</Text>
                  <div style={{ marginTop: 8, display: 'flex', gap: 8, flexWrap: 'wrap' }}>
                    <Tag color="blue">
                      初始赠送: {((quota.quota_sources.free_granted ?? 0) || 0).toFixed(1)} h
                    </Tag>
                    <Tag color="green">
                      充值获得: {((quota.quota_sources.recharge ?? 0) || 0).toFixed(1)} h
                    </Tag>
                    <Tag color="gold">
                      管理员赠送: {((quota.quota_sources.admin_granted ?? 0) || 0).toFixed(1)} h
                    </Tag>
                  </div>
                  <div style={{ marginTop: 8 }}>
                    <Text type="secondary" style={{ fontSize: 12 }}>
                      总核时: {formatCpuHours(quota.total_cpu_hours)} h =
                      {formatCpuHours(quota.quota_sources.free_granted)} +
                      {formatCpuHours(quota.quota_sources.recharge)} +
                      {formatCpuHours(quota.quota_sources.admin_granted)} -
                      {formatCpuHours(quota.used_cpu_hours)} (已消费)
                    </Text>
                  </div>
                </div>
              )}

              {/* 配额使用进度 */}
              {quota && (
                <div style={{ marginTop: 24 }}>
                  <Text strong>配额使用情况</Text>
                  <div style={{ marginTop: 8 }}>
                    <Progress
                      percent={Math.round(((quota.used_cpu_hours || 0) / (quota.total_cpu_hours || 1)) * 100)}
                      status={quotaStatus.status === 'critical' ? 'exception' :
                             quotaStatus.status === 'warning' ? 'active' : 'success'}
                    />
                    <div style={{ display: 'flex', justifyContent: 'space-between', marginTop: 8, fontSize: 12 }}>
                      <Text type="secondary">
                        已使用: {formatCpuHours(quota.used_cpu_hours)} h / {formatCpuHours(quota.total_cpu_hours)} h
                      </Text>
                      <Text type="secondary">
                        {Math.round(((quota.used_cpu_hours || 0) / (quota.total_cpu_hours || 1)) * 100)}%
                      </Text>
                    </div>
                  </div>
                </div>
              )}

              {/* 核时系统说明 */}
              <Alert
                message="核时系统说明"
                description={
                  (quota?.balance_cpu_hours ?? 0) < 0
                    ? `您当前欠费 ${safeToFixed(quota?.balance_cpu_hours, 2)} 核时。欠费期间无法提交新任务。请立即充值以继续使用服务。`
                    : `您当前可用余额为 ${safeToFixed(quota?.balance_cpu_hours, 2)} 核时。充值会增加可用余额，消费会减少可用余额。`
                }
                type={(quota?.balance_cpu_hours ?? 0) < 0 ? 'error' : 'info'}
                showIcon
                style={{ marginTop: 24 }}
              />
            </Card>
          </Col>
        </Row>
      </div>
    );
  };

  // 充值中心 Tab
  const renderRechargeTab = () => {
    const accountConfig = getAccountConfig();

    return (
      <div>
        {/* 余额概览 */}
        <Row gutter={16} style={{ marginBottom: 24 }}>
          <Col xs={24} sm={12} lg={6}>
            <Card style={{ borderRadius: 12, boxShadow: '0 2px 8px rgba(0,0,0,0.06)', border: 'none' }}>
              <Statistic
                title="可用余额"
                value={balance?.available || quota?.balance_cpu_hours || 0}
                suffix="核时"
                precision={2}
                valueStyle={{ color: (balance?.available || quota?.balance_cpu_hours || 0) > 10 ? '#52c41a' : '#ff4d4f' }}
              />
            </Card>
          </Col>
          <Col xs={24} sm={12} lg={6}>
            <Card style={{ borderRadius: 12, boxShadow: '0 2px 8px rgba(0,0,0,0.06)', border: 'none' }}>
              <Statistic
                title="冻结中"
                value={balance?.frozen || quota?.frozen_cpu_hours || 0}
                suffix="核时"
                precision={2}
              />
            </Card>
          </Col>
          <Col xs={24} sm={12} lg={6}>
            <Card style={{ borderRadius: 12, boxShadow: '0 2px 8px rgba(0,0,0,0.06)', border: 'none' }}>
              <Statistic
                title="欠费"
                value={balance?.debt || 0}
                suffix="核时"
                precision={2}
                valueStyle={{ color: (balance?.debt || 0) > 0 ? '#ff4d4f' : '#52c41a' }}
              />
            </Card>
          </Col>
          <Col xs={24} sm={12} lg={6}>
            <Card style={{ borderRadius: 12, boxShadow: '0 2px 8px rgba(0,0,0,0.06)', border: 'none' }}>
              <Statistic
                title="当前单价"
                value={balance?.price_per_hour || 0.1}
                prefix="¥"
                suffix="/核时"
                precision={QUOTA_PRECISION}
              />
            </Card>
          </Col>
        </Row>

        {/* 欠费警告 */}
        {(balance?.debt || 0) > 0 && (
          <Alert
            type="error"
            showIcon
            message="您有欠费未还清"
            description={`当前欠费 ${(balance?.debt || 0).toFixed(2)} 核时，请尽快充值。欠费期间无法提交新任务。`}
            style={{ marginBottom: 24, borderRadius: 8 }}
          />
        )}

        {/* 充值说明 */}
        <Alert
          message="充值说明"
          description={`您的当前单价为 ¥${(balance?.price_per_hour || 0.1).toFixed(4)}/核时。充值金额将按此单价转换为核时数。如需修改单价，请联系管理员。`}
          type="info"
          showIcon
          style={{ marginBottom: 24, borderRadius: 8 }}
        />

        {/* 充值表单 */}
        <Card
          title="立即充值"
          extra={<Text type="secondary">1元 = {(1 / (balance?.price_per_hour || 0.1)).toFixed(0)} 核时</Text>}
          style={{ borderRadius: 12, boxShadow: '0 2px 8px rgba(0,0,0,0.06)', border: 'none' }}
        >
          <div style={{ marginBottom: 16 }}>
            <Text strong>选择金额：</Text>
            <div style={{ marginTop: 8 }}>
              <Radio.Group value={rechargeAmount} onChange={e => setRechargeAmount(e.target.value)}>
                <Space wrap>
                  {PRESET_AMOUNTS.map(a => (
                    <Radio.Button key={a} value={a}>¥{a}</Radio.Button>
                  ))}
                </Space>
              </Radio.Group>
            </div>
          </div>

          <div style={{ marginBottom: 16 }}>
            <Text strong>自定义金额：</Text>
            <InputNumber
              style={{ width: 200, marginLeft: 8 }}
              min={10}
              max={10000}
              value={rechargeAmount}
              onChange={v => setRechargeAmount(v || 10)}
              prefix="¥"
              precision={0}
            />
          </div>

          <div style={{ marginBottom: 24 }}>
            <Text strong>支付方式：</Text>
            <div style={{ marginTop: 8 }}>
              <Radio.Group value={paymentMethod} onChange={e => setPaymentMethod(e.target.value)}>
                <Space>
                  <Radio.Button value="simulated">模拟支付</Radio.Button>
                  <Radio.Button value="wechat" disabled>微信支付</Radio.Button>
                  <Radio.Button value="alipay" disabled>支付宝</Radio.Button>
                </Space>
              </Radio.Group>
            </div>
            <Text type="secondary" style={{ fontSize: 12, marginTop: 4, display: 'block' }}>
              * 微信/支付宝支付即将上线
            </Text>
          </div>

          <Divider />

          <div style={{ textAlign: 'center' }}>
            <div style={{ marginBottom: 16 }}>
              <Text>充值金额：</Text>
              <Text strong style={{ fontSize: 24, color: '#1677ff' }}>¥{rechargeAmount}</Text>
              <Text style={{ marginLeft: 16 }}>可获得：</Text>
              <Text strong style={{ fontSize: 24, color: '#52c41a' }}>
                {(rechargeAmount / (balance?.price_per_hour || 0.1)).toFixed(2)} 核时
              </Text>
            </div>
            <Button
              type="primary"
              size="large"
              loading={rechargeLoading}
              onClick={handleRecharge}
              style={{ borderRadius: 8, boxShadow: '0 2px 8px rgba(22, 119, 255, 0.3)' }}
            >
              立即充值
            </Button>
          </div>
        </Card>
      </div>
    );
  };

  // 账户设置 Tab
  const renderSettingsTab = () => {
    return (
      <div>
        <Row gutter={24}>
          <Col xs={24} lg={12}>
            <Card title="个人信息" style={{ marginBottom: 24 }}>
              <Form
                form={form}
                layout="vertical"
                onFinish={handleUpdateProfile}
                initialValues={{
                  username: user?.username,
                  email: user?.email
                }}
              >
                <Form.Item
                  label="用户名"
                  name="username"
                  rules={[{ required: true, message: '请输入用户名' }]}
                >
                  <Input />
                </Form.Item>
                <Form.Item
                  label="邮箱"
                  name="email"
                  rules={[
                    { required: true, message: '请输入邮箱' },
                    { type: 'email', message: '请输入有效的邮箱地址' }
                  ]}
                >
                  <Input />
                </Form.Item>
                <Form.Item label="真实姓名" name="full_name">
                  <Input placeholder="请输入真实姓名" />
                </Form.Item>
                <Form.Item
                  label="手机号"
                  name="phone"
                  rules={[
                    { pattern: /^1[3-9]\d{9}$/, message: '请输入有效的手机号' }
                  ]}
                >
                  <Input placeholder="请输入手机号" />
                </Form.Item>
                <Form.Item>
                  <Button type="primary" htmlType="submit">保存修改</Button>
                </Form.Item>
              </Form>
            </Card>
          </Col>

          <Col xs={24} lg={12}>
            <Card title="安全设置" style={{ marginBottom: 24 }}>
              <Form
                form={passwordForm}
                layout="vertical"
                onFinish={handleChangePassword}
              >
                <Form.Item
                  label="当前密码"
                  name="oldPassword"
                  rules={[{ required: true, message: '请输入当前密码' }]}
                >
                  <Input.Password />
                </Form.Item>
                <Form.Item
                  label="新密码"
                  name="newPassword"
                  rules={[
                    { required: true, message: '请输入新密码' },
                    { min: 6, message: '密码长度至少6位' }
                  ]}
                >
                  <Input.Password />
                </Form.Item>
                <Form.Item
                  label="确认新密码"
                  name="confirmPassword"
                  rules={[
                    { required: true, message: '请确认新密码' },
                    ({ getFieldValue }) => ({
                      validator(_, value) {
                        if (!value || getFieldValue('newPassword') === value) {
                          return Promise.resolve();
                        }
                        return Promise.reject(new Error('两次输入的密码不一致'));
                      },
                    }),
                  ]}
                >
                  <Input.Password />
                </Form.Item>
                <Form.Item>
                  <Button type="primary" htmlType="submit" icon={<LockOutlined />}>
                    修改密码
                  </Button>
                </Form.Item>
              </Form>
            </Card>
          </Col>
        </Row>
      </div>
    );
  };

  // 账号管理 Tab
  const renderOrganizationTab = () => {
    if (!accountInfo) return null;

    const { account_type } = accountInfo;

    // 个人用户不显示账号管理
    if (account_type === 'PERSONAL' || account_type === 'personal') {
      return (
        <Alert
          message="个人用户"
          description="个人用户无需账号管理功能。"
          type="info"
          showIcon
        />
      );
    }

    // 主账号显示子账号管理面板
    if (account_type === 'MASTER_ACCOUNT' || account_type === 'master_account') {
      return <SubAccountManagementPanel compact={true} />;
    }

    // 子账号显示主账号信息
    if (account_type === 'SUB_ACCOUNT' || account_type === 'sub_account') {
      // 计算总可用配额
      const personalBalance = subAccountInfo?.balance_cpu_hours || 0;
      const allocatedQuota = subAccountInfo?.allocated_quota || 0;
      const totalAvailableQuota = personalBalance + allocatedQuota;

      return (
        <Spin spinning={subAccountLoading}>
          <div>
            {!subAccountInfo ? (
              <Alert
                message="加载中..."
                description="正在加载您的子账号信息，请稍候"
                type="info"
                showIcon
                style={{ marginBottom: 24 }}
              />
            ) : null}

            {/* 主账号信息卡片 */}
            {subAccountInfo && (
              <Card title="主账号信息" style={{ marginBottom: 24 }}>
                <Row gutter={[16, 16]}>
                  <Col xs={24} sm={12} md={6}>
                    <Card style={{ borderRadius: 12, boxShadow: '0 2px 8px rgba(0,0,0,0.06)', border: 'none' }}>
                      <Statistic
                        title="主账号"
                        value={subAccountInfo.master_username || '未知'}
                        prefix={<CrownOutlined />}
                      />
                    </Card>
                  </Col>
                  <Col xs={24} sm={12} md={6}>
                    <Card style={{ borderRadius: 12, boxShadow: '0 2px 8px rgba(0,0,0,0.06)', border: 'none' }}>
                      <Statistic
                        title="主账号邮箱"
                        value={subAccountInfo.master_email || '未知'}
                        prefix={<UserOutlined />}
                      />
                    </Card>
                  </Col>
                  <Col xs={24} sm={12} md={6}>
                    <Card style={{ borderRadius: 12, boxShadow: '0 2px 8px rgba(0,0,0,0.06)', border: 'none' }}>
                      <Statistic
                        title="分配配额"
                        value={formatCpuHours(allocatedQuota)}
                        suffix="核时"
                        prefix={<WalletOutlined />}
                        valueStyle={{ color: '#1890ff' }}
                      />
                    </Card>
                  </Col>
                  <Col xs={24} sm={12} md={6}>
                    <Card style={{ borderRadius: 12, boxShadow: '0 2px 8px rgba(0,0,0,0.06)', border: 'none' }}>
                      <Statistic
                        title="个人余额"
                        value={formatCpuHours(personalBalance)}
                        suffix="核时"
                        prefix={<WalletOutlined />}
                        valueStyle={{ color: personalBalance > 0 ? '#52c41a' : '#ff4d4f' }}
                      />
                    </Card>
                  </Col>
                </Row>

                {/* 总可用配额 */}
                <Divider style={{ margin: '20px 0' }} />
                <Row gutter={[16, 16]}>
                  <Col xs={24}>
                    <Card
                      style={{
                        borderRadius: 12,
                        boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
                        border: 'none',
                        background: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)'
                      }}
                    >
                      <div style={{ color: '#fff' }}>
                        <div style={{ fontSize: '12px', opacity: 0.9, marginBottom: '8px' }}>
                          📊 总可用配额
                        </div>
                        <div style={{ fontSize: '28px', fontWeight: 'bold', marginBottom: '4px' }}>
                          {formatCpuHours(totalAvailableQuota)}
                        </div>
                        <div style={{ fontSize: '12px', opacity: 0.85 }}>
                          核时 = {formatCpuHours(personalBalance)} (个人) + {formatCpuHours(allocatedQuota)} (分配)
                        </div>
                      </div>
                    </Card>
                  </Col>
                </Row>
              </Card>
            )}

            {/* 配额说明 */}
            <Card title="配额说明" style={{ marginBottom: 24 }}>
              <Alert
                message="子账号配额来源"
                description={
                  <div>
                    <p>您的可用配额来自两个部分：</p>
                    <ul style={{ marginBottom: 0 }}>
                      <li><strong>个人余额</strong>：您自己充值获得的核时</li>
                      <li><strong>分配配额</strong>：主账号管理员分配给您的核时</li>
                    </ul>
                    <p style={{ marginTop: 12, marginBottom: 0 }}>
                      <strong>总可用配额 = 个人余额 + 分配配额 = {formatCpuHours(totalAvailableQuota)} 核时</strong>
                    </p>
                  </div>
                }
                type="info"
                showIcon
              />
            </Card>
          </div>
        </Spin>
      );
    }

    return null;
  };

  // 使用统计 Tab
  const renderStatisticsTab = () => {
    // 计算汇总统计
    const totalStats = dailyStats.reduce((acc, day) => ({
      jobs_submitted: acc.jobs_submitted + (day.jobs_submitted || 0),
      jobs_completed: acc.jobs_completed + (day.jobs_completed || 0),
      jobs_failed: acc.jobs_failed + (day.jobs_failed || 0),
      jobs_cancelled: acc.jobs_cancelled + (day.jobs_cancelled || 0),
      cpu_hours_used: acc.cpu_hours_used + (day.cpu_hours_used || 0),
      cluster_analysis_cpu_hours: acc.cluster_analysis_cpu_hours + (day.cluster_analysis_cpu_hours || 0),
      cluster_analysis_task_count: acc.cluster_analysis_task_count + (day.cluster_analysis_task_count || 0),
      storage_used_gb: acc.storage_used_gb + (day.storage_used_gb || 0),
      max_concurrent_jobs: Math.max(acc.max_concurrent_jobs, day.max_concurrent_jobs || 0),
    }), {
      jobs_submitted: 0,
      jobs_completed: 0,
      jobs_failed: 0,
      jobs_cancelled: 0,
      cpu_hours_used: 0,
      cluster_analysis_cpu_hours: 0,
      cluster_analysis_task_count: 0,
      storage_used_gb: 0,
      max_concurrent_jobs: 0,
    });

    const successRate = totalStats.jobs_submitted > 0
      ? ((totalStats.jobs_completed / totalStats.jobs_submitted) * 100).toFixed(1)
      : '0';

    // 日期表格列定义
    const dailyTableColumns = [
      {
        title: '日期',
        dataIndex: 'date',
        key: 'date',
        render: (date: string) => new Date(date).toLocaleDateString('zh-CN'),
        width: 100,
      },
      {
        title: '提交',
        dataIndex: 'jobs_submitted',
        key: 'jobs_submitted',
        align: 'right' as const,
        width: 80,
      },
      {
        title: '完成',
        dataIndex: 'jobs_completed',
        key: 'jobs_completed',
        align: 'right' as const,
        width: 80,
      },
      {
        title: '失败',
        dataIndex: 'jobs_failed',
        key: 'jobs_failed',
        align: 'right' as const,
        width: 80,
      },
      {
        title: '取消',
        dataIndex: 'jobs_cancelled',
        key: 'jobs_cancelled',
        align: 'right' as const,
        width: 80,
      },
      {
        title: 'MD核时',
        dataIndex: 'cpu_hours_used',
        key: 'cpu_hours_used',
        align: 'right' as const,
        width: 100,
        render: (value: number) => value ? value.toFixed(2) : '0.00',
      },
      {
        title: '分析核时',
        dataIndex: 'cluster_analysis_cpu_hours',
        key: 'cluster_analysis_cpu_hours',
        align: 'right' as const,
        width: 100,
        render: (value: number) => value ? value.toFixed(2) : '0.00',
      },
      {
        title: '分析任务',
        dataIndex: 'cluster_analysis_task_count',
        key: 'cluster_analysis_task_count',
        align: 'right' as const,
        width: 100,
      },
      {
        title: '存储(GB)',
        dataIndex: 'storage_used_gb',
        key: 'storage_used_gb',
        align: 'right' as const,
        width: 100,
        render: (value: number) => value ? value.toFixed(2) : '0.00',
      },
    ];

    return (
      <div>
        {/* 统计概览卡片 */}
        <Card title={`使用统计概览 (${statsPeriod}天)`} style={{ marginBottom: 24 }}>
          <div style={{ marginBottom: 20 }}>
            <Space>
              <Text strong>统计周期：</Text>
              <Radio.Group value={statsPeriod} onChange={e => setStatsPeriod(e.target.value)}>
                <Radio.Button value={7}>最近7天</Radio.Button>
                <Radio.Button value={30}>最近30天</Radio.Button>
                <Radio.Button value={90}>最近90天</Radio.Button>
              </Radio.Group>
              <Button icon={<ReloadOutlined />} onClick={loadDailyStats} loading={statsLoading}>
                刷新
              </Button>
            </Space>
          </div>

          {dailyStats.length === 0 ? (
            <Alert
              message="暂无统计数据"
              description="该时间段内没有任何使用记录"
              type="info"
              showIcon
              style={{ marginBottom: 20 }}
            />
          ) : null}

          <Row gutter={[16, 16]}>
            <Col xs={24} sm={12} md={6}>
              <Card bordered={false} style={{ textAlign: 'center', borderRadius: 8 }}>
                <Statistic
                  title="提交任务"
                  value={totalStats.jobs_submitted}
                  suffix="个"
                  valueStyle={{ color: '#1890ff' }}
                  prefix={<RocketOutlined />}
                />
              </Card>
            </Col>
            <Col xs={24} sm={12} md={6}>
              <Card bordered={false} style={{ textAlign: 'center', borderRadius: 8 }}>
                <Statistic
                  title="完成任务"
                  value={totalStats.jobs_completed}
                  suffix="个"
                  valueStyle={{ color: '#52c41a' }}
                  prefix={<CheckCircleOutlined />}
                />
              </Card>
            </Col>
            <Col xs={24} sm={12} md={6}>
              <Card bordered={false} style={{ textAlign: 'center', borderRadius: 8 }}>
                <Statistic
                  title="成功率"
                  value={successRate}
                  suffix="%"
                  valueStyle={{ color: totalStats.jobs_completed > 0 ? '#52c41a' : '#faad14' }}
                  prefix={<BarChartOutlined />}
                />
              </Card>
            </Col>
            <Col xs={24} sm={12} md={6}>
              <Card bordered={false} style={{ textAlign: 'center', borderRadius: 8 }}>
                <Statistic
                  title="消耗核时"
                  value={totalStats.cpu_hours_used.toFixed(2)}
                  suffix="h"
                  valueStyle={{ color: '#722ed1' }}
                  prefix={<ThunderboltOutlined />}
                />
              </Card>
            </Col>
          </Row>

          <Divider />

          <Row gutter={[16, 16]} style={{ marginTop: 16 }}>
            <Col xs={24} sm={12} md={6}>
              <Card bordered={false} style={{ textAlign: 'center', borderRadius: 8 }}>
                <Statistic
                  title="分析核时"
                  value={totalStats.cluster_analysis_cpu_hours.toFixed(2)}
                  suffix="h"
                  valueStyle={{ color: '#fa8c16' }}
                  prefix={<LineChartOutlined />}
                />
              </Card>
            </Col>
            <Col xs={24} sm={12} md={6}>
              <Card bordered={false} style={{ textAlign: 'center', borderRadius: 8 }}>
                <Statistic
                  title="分析任务"
                  value={totalStats.cluster_analysis_task_count}
                  suffix="个"
                  valueStyle={{ color: '#13c2c2' }}
                  prefix={<ClusterOutlined />}
                />
              </Card>
            </Col>
            <Col xs={24} sm={12} md={6}>
              <Card bordered={false} style={{ textAlign: 'center', borderRadius: 8 }}>
                <Statistic
                  title="总存储"
                  value={totalStats.storage_used_gb.toFixed(2)}
                  suffix="GB"
                  valueStyle={{ color: '#eb2f96' }}
                  prefix={<DatabaseOutlined />}
                />
              </Card>
            </Col>
            <Col xs={24} sm={12} md={6}>
              <Card bordered={false} style={{ textAlign: 'center', borderRadius: 8 }}>
                <Statistic
                  title="最大并发"
                  value={totalStats.max_concurrent_jobs}
                  suffix="个"
                  valueStyle={{ color: '#2f54eb' }}
                  prefix={<ClusterOutlined />}
                />
              </Card>
            </Col>
          </Row>
        </Card>

        {/* 每日详情表格 */}
        <Card title="每日统计详情" style={{ marginBottom: 24 }}>
          <Table
            columns={dailyTableColumns}
            dataSource={dailyStats.map((item, index) => ({ ...item, key: index }))}
            rowKey="key"
            pagination={{ pageSize: 15, showSizeChanger: true }}
            scroll={{ x: 1200 }}
            loading={statsLoading}
          />
        </Card>
      </div>
    );
  };

  if (loading) {
    return (
      <div style={{ padding: 24, textAlign: 'center' }}>
        <Spin size="large" />
        <div style={{ marginTop: 16 }}>
          <Text type="secondary">加载账户信息中...</Text>
        </div>
      </div>
    );
  }

  return (
    <div style={{
      padding: 24,
      background: isDark ? '#141414' : '#f5f7fa',
      minHeight: 'calc(100vh - 64px - 48px)', // 减去 header 和 footer
      paddingBottom: 48, // 额外底部间距
    }}>
      <div style={{ marginBottom: 24 }}>
        <Title level={2} style={{ margin: 0 }}>账户中心</Title>
        <Text type="secondary">管理您的账户信息、配额和设置</Text>
      </div>

      <Card>
        <Tabs
          activeKey={activeTab}
          onChange={setActiveTab}
          items={[
            {
              key: 'overview',
              label: (
                <span>
                  <UserOutlined />
                  账户总览
                </span>
              ),
              children: renderOverviewTab(),
            },
            {
              key: 'recharge',
              label: (
                <span>
                  <WalletOutlined />
                  充值中心
                </span>
              ),
              children: renderRechargeTab(),
            },
            // 账号管理 Tab - 根据账户类型显示
            ...(accountInfo && accountInfo.account_type !== 'PERSONAL' && accountInfo.account_type !== 'personal' ? [{
              key: 'organization',
              label: (
                <span>
                  {(accountInfo.account_type === 'MASTER_ACCOUNT' || accountInfo.account_type === 'master_account') ? <BankOutlined /> : <TeamOutlined />}
                  {(accountInfo.account_type === 'MASTER_ACCOUNT' || accountInfo.account_type === 'master_account') ? '子账号管理' : '账号信息'}
                </span>
              ),
              children: renderOrganizationTab(),
            }] : []),
            // 管理员专用 Tabs 已移至管理面板
            {
              key: 'settings',
              label: (
                <span>
                  <SettingOutlined />
                  账户设置
                </span>
              ),
              children: renderSettingsTab(),
            },
            {
              key: 'statistics',
              label: (
                <span>
                  <BarChartOutlined />
                  使用统计
                </span>
              ),
              children: renderStatisticsTab(),
            },
          ]}
        />
      </Card>
    </div>
  );
};

export default AccountCenter;
