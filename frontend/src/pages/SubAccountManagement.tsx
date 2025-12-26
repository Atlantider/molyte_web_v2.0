/**
 * 子账号管理页面
 * 
 * 仅主账号用户可以访问，用于管理子账号
 */

import React, { useState, useEffect } from 'react';
import {
  Card,
  Table,
  Button,
  Space,
  Modal,
  Form,
  Input,
  InputNumber,
  message,
  Spin,
  Tag,
  Statistic,
  Row,
  Col,
  Popconfirm,
  Alert,
  Typography,
  Divider
} from 'antd';
import {
  PlusOutlined,
  EditOutlined,
  DeleteOutlined,
  UserOutlined,
  WalletOutlined,
  TeamOutlined,
  ReloadOutlined
} from '@ant-design/icons';
import { useNavigate } from 'react-router-dom';
import { useThemeStore } from '../stores/themeStore';
import {
  getMySubAccounts,
  createMySubAccount,
  deleteMySubAccount,
  SubAccount,
  CreateSubAccountRequest,
} from '../api/accounts';

const { Title, Text } = Typography;

const SubAccountManagement: React.FC = () => {
  const navigate = useNavigate();
  const { mode } = useThemeStore();
  const isDark = mode === 'dark';
  
  const [subAccounts, setSubAccounts] = useState<SubAccount[]>([]);
  const [loading, setLoading] = useState(true);
  const [createModalVisible, setCreateModalVisible] = useState(false);
  const [form] = Form.useForm();

  // 加载子账号列表
  const loadSubAccounts = async () => {
    try {
      setLoading(true);
      const data = await getMySubAccounts();
      setSubAccounts(data);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '加载子账号列表失败');
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    loadSubAccounts();
  }, []);

  // 创建子账号
  const handleCreateSubAccount = async (values: CreateSubAccountRequest) => {
    try {
      await createMySubAccount(values);
      message.success('子账号创建成功');
      setCreateModalVisible(false);
      form.resetFields();
      loadSubAccounts();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '创建子账号失败');
    }
  };

  // 删除子账号
  const handleDeleteSubAccount = async (id: number) => {
    try {
      await deleteMySubAccount(id);
      message.success('子账号删除成功');
      loadSubAccounts();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '删除子账号失败');
    }
  };

  // 表格列定义
  const columns = [
    {
      title: '用户ID',
      dataIndex: 'user_id',
      key: 'user_id',
      width: 100,
    },
    {
      title: '用户名',
      dataIndex: 'username',
      key: 'username',
      render: (text: string) => (
        <Space>
          <UserOutlined />
          {text || '-'}
        </Space>
      ),
    },
    {
      title: '邮箱',
      dataIndex: 'email',
      key: 'email',
    },
    {
      title: '个人配额',
      dataIndex: 'personal_quota',
      key: 'personal_quota',
      render: (value: number | null) => (
        <Tag color="blue">{value === null ? '无限制' : value.toFixed(2) + ' 核时'}</Tag>
      ),
    },
    {
      title: '已用配额',
      dataIndex: 'personal_used',
      key: 'personal_used',
      render: (value: number) => (
        <Tag color="orange">{value.toFixed(2)} 核时</Tag>
      ),
    },
    {
      title: '可用配额',
      key: 'balance',
      render: (record: SubAccount) => {
        const balance = record.balance_cpu_hours;
        return (
          <Tag color={balance > 10 ? 'green' : balance > 1 ? 'orange' : 'red'}>
            {balance.toFixed(2)} 核时
          </Tag>
        );
      },
    },
    {
      title: '状态',
      dataIndex: 'is_active',
      key: 'is_active',
      render: (active: boolean) => (
        <Tag color={active ? 'green' : 'red'}>
          {active ? '活跃' : '停用'}
        </Tag>
      ),
    },
    {
      title: '创建时间',
      dataIndex: 'created_at',
      key: 'created_at',
      render: (date: string) => new Date(date).toLocaleDateString(),
    },
    {
      title: '操作',
      key: 'actions',
      render: (record: SubAccount) => (
        <Space>
          <Button
            type="text"
            icon={<EditOutlined />}
            size="small"
            onClick={() => {
              // TODO: 实现编辑功能
              message.info('编辑功能开发中');
            }}
          >
            编辑
          </Button>
          <Popconfirm
            title="确定要删除这个子账号吗？"
            description="删除后无法恢复，请谨慎操作。"
            onConfirm={() => handleDeleteSubAccount(record.id)}
            okText="确定"
            cancelText="取消"
          >
            <Button
              type="text"
              danger
              icon={<DeleteOutlined />}
              size="small"
            >
              删除
            </Button>
          </Popconfirm>
        </Space>
      ),
    },
  ];

  // 计算统计数据
  const totalPersonalQuota = subAccounts.reduce((sum, account) => sum + (account.personal_quota || 0), 0);
  const totalPersonalUsed = subAccounts.reduce((sum, account) => sum + (account.personal_used || 0), 0);
  const totalBalance = subAccounts.reduce((sum, account) => sum + (account.balance_cpu_hours || 0), 0);
  const activeCount = subAccounts.filter(account => account.is_active).length;

  return (
    <div style={{ padding: 24, background: isDark ? '#141414' : '#f5f7fa', minHeight: '100vh' }}>
      {/* 页面标题 */}
      <div style={{ marginBottom: 24, display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
        <div>
          <Title level={2} style={{ margin: 0, marginBottom: 8 }}>子账号管理</Title>
          <Text type="secondary">管理您的子账号和配额分配</Text>
        </div>
        <Space>
          <Button icon={<ReloadOutlined />} onClick={loadSubAccounts}>
            刷新
          </Button>
          <Button 
            type="primary" 
            icon={<PlusOutlined />} 
            onClick={() => setCreateModalVisible(true)}
          >
            创建子账号
          </Button>
        </Space>
      </div>

      {/* 统计卡片 */}
      <Row gutter={[16, 16]} style={{ marginBottom: 24 }}>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic
              title="子账号总数"
              value={subAccounts.length}
              prefix={<UserOutlined />}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic
              title="活跃账号"
              value={activeCount}
              prefix={<TeamOutlined />}
              valueStyle={{ color: '#52c41a' }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic
              title="总个人配额"
              value={totalPersonalQuota}
              suffix="核时"
              prefix={<WalletOutlined />}
              valueStyle={{ color: '#1890ff' }}
              precision={2}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic
              title="总可用配额"
              value={totalBalance}
              suffix="核时"
              prefix={<WalletOutlined />}
              valueStyle={{
                color: totalBalance > 50 ? '#52c41a' : totalBalance > 10 ? '#faad14' : '#ff4d4f'
              }}
              precision={2}
            />
          </Card>
        </Col>
      </Row>

      {/* 配额警告 */}
      {totalBalance < 10 && (
        <Alert
          message="配额不足警告"
          description={`您的子账号总可用配额不足 10 核时，建议及时充值。`}
          type="warning"
          showIcon
          style={{ marginBottom: 24 }}
          action={
            <Button size="small" onClick={() => navigate('/workspace/account-center')}>
              去充值
            </Button>
          }
        />
      )}

      {/* 子账号列表 */}
      <Card title="子账号列表">
        <Table
          columns={columns}
          dataSource={subAccounts}
          rowKey="id"
          loading={loading}
          pagination={{
            pageSize: 10,
            showSizeChanger: true,
            showQuickJumper: true,
            showTotal: (total) => `共 ${total} 个子账号`,
          }}
        />
      </Card>

      {/* 创建子账号对话框 */}
      <Modal
        title="创建子账号"
        open={createModalVisible}
        onOk={() => form.submit()}
        onCancel={() => {
          setCreateModalVisible(false);
          form.resetFields();
        }}
        width={600}
      >
        <Form
          form={form}
          layout="vertical"
          onFinish={handleCreateSubAccount}
        >
          <Form.Item
            label="用户名"
            name="username"
            rules={[{ required: true, message: '请输入用户名' }]}
            extra="子账号的登录用户名"
          >
            <Input
              placeholder="用户名"
            />
          </Form.Item>

          <Form.Item
            label="邮箱"
            name="email"
            rules={[
              { required: true, message: '请输入邮箱' },
              { type: 'email', message: '请输入有效的邮箱地址' }
            ]}
            extra="子账号的邮箱地址"
          >
            <Input
              placeholder="邮箱地址"
              type="email"
            />
          </Form.Item>

          <Form.Item
            label="密码"
            name="password"
            rules={[
              { required: true, message: '请输入密码' },
              { min: 6, message: '密码至少6个字符' }
            ]}
            extra="子账号的登录密码"
          >
            <Input.Password
              placeholder="密码"
            />
          </Form.Item>

          <Form.Item
            label="个人配额（可选）"
            name="personal_quota"
            extra="为该子账号设置的个人配额限制，留空表示无限制"
          >
            <InputNumber
              style={{ width: '100%' }}
              placeholder="配额（核时）"
              min={0.1}
              step={0.1}
              precision={2}
            />
          </Form.Item>

          <Divider />

          <Alert
            message="创建说明"
            description={
              <ul style={{ margin: 0, paddingLeft: 20 }}>
                <li>将创建一个新的子账号用户</li>
                <li>子账号可以使用自己的用户名和密码登录</li>
                <li>可以为子账号设置个人配额限制</li>
                <li>子账号的配额来自主账号的配额</li>
                <li>您可以随时调整或删除子账号</li>
              </ul>
            }
            type="info"
            showIcon
          />
        </Form>
      </Modal>
    </div>
  );
};

export default SubAccountManagement;
