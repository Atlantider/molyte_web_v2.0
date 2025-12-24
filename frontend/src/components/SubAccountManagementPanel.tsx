/**
 * 子账号管理面板组件
 * 
 * 可以嵌入到 AccountCenter 或作为独立页面使用
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
  Divider,
  Switch
} from 'antd';
import {
  PlusOutlined,
  EditOutlined,
  DeleteOutlined,
  UserOutlined,
  WalletOutlined,
  TeamOutlined,
  ReloadOutlined,
  SearchOutlined,
  EyeOutlined,
  UnorderedListOutlined,
  UserAddOutlined
} from '@ant-design/icons';
import { useThemeStore } from '../stores/themeStore';
import { useNavigate } from 'react-router-dom';
import {
  getMySubAccounts,
  createMySubAccount,
  deleteMySubAccount,
  updateMySubAccount,
  addExistingUserAsSubAccount,
  getSubAccountJobs,
  getAllSubAccountsJobs,
  SubAccount,
  CreateSubAccountRequest,
  UpdateSubAccountRequest,
  SubAccountJob,
} from '../api/accounts';

const { Title, Text } = Typography;

interface SubAccountManagementPanelProps {
  compact?: boolean; // 是否为紧凑模式（嵌入到 tab 中）
}

const SubAccountManagementPanel: React.FC<SubAccountManagementPanelProps> = ({ compact = false }) => {
  const { mode } = useThemeStore();
  const isDark = mode === 'dark';
  const navigate = useNavigate();

  const [subAccounts, setSubAccounts] = useState<SubAccount[]>([]);
  const [loading, setLoading] = useState(true);
  const [createModalVisible, setCreateModalVisible] = useState(false);
  const [addExistingModalVisible, setAddExistingModalVisible] = useState(false);
  const [editModalVisible, setEditModalVisible] = useState(false);
  const [jobsModalVisible, setJobsModalVisible] = useState(false);
  const [editingSubAccount, setEditingSubAccount] = useState<SubAccount | null>(null);
  const [viewingSubAccount, setViewingSubAccount] = useState<SubAccount | null>(null);
  const [subAccountJobs, setSubAccountJobs] = useState<SubAccountJob[]>([]);
  const [jobsLoading, setJobsLoading] = useState(false);
  const [jobsTotal, setJobsTotal] = useState(0);
  const [form] = Form.useForm();
  const [addExistingForm] = Form.useForm();
  const [editForm] = Form.useForm();

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

  // 添加现有用户为子账号
  const handleAddExistingUser = async (values: { username_or_email: string; allocated_quota?: number }) => {
    try {
      await addExistingUserAsSubAccount(values);
      message.success('成功添加现有用户为子账号');
      setAddExistingModalVisible(false);
      addExistingForm.resetFields();
      loadSubAccounts();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '添加失败');
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

  // 打开编辑对话框
  const handleEditClick = (subAccount: SubAccount) => {
    setEditingSubAccount(subAccount);
    editForm.setFieldsValue({
      personal_quota: subAccount.personal_quota,
      is_active: subAccount.is_active,
    });
    setEditModalVisible(true);
  };

  // 查看子账号任务
  const handleViewJobs = async (subAccount: SubAccount) => {
    setViewingSubAccount(subAccount);
    setJobsModalVisible(true);
    setJobsLoading(true);
    try {
      const result = await getSubAccountJobs(subAccount.id, 0, 50);
      setSubAccountJobs(result.jobs);
      setJobsTotal(result.total);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '加载任务列表失败');
    } finally {
      setJobsLoading(false);
    }
  };

  // 更新子账号
  const handleUpdateSubAccount = async (values: UpdateSubAccountRequest) => {
    if (!editingSubAccount) return;
    try {
      await updateMySubAccount(editingSubAccount.id, values);
      message.success('子账号更新成功');
      setEditModalVisible(false);
      editForm.resetFields();
      setEditingSubAccount(null);
      loadSubAccounts();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '更新子账号失败');
    }
  };

  // 表格列定义
  const columns = [
    {
      title: '用户名',
      dataIndex: 'username',
      key: 'username',
      render: (text: string, record: SubAccount) => (
        <Button
          type="link"
          onClick={() => navigate(`/workspace/sub-accounts/${record.id}`)}
          style={{ padding: 0 }}
        >
          <Text strong>{text}</Text>
        </Button>
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
      render: (value: number | null | undefined) => (
        <Tag color="blue">{value === null || value === undefined ? '无限制' : value.toFixed(2) + ' 核时'}</Tag>
      ),
    },
    {
      title: '已用配额',
      dataIndex: 'personal_used',
      key: 'personal_used',
      render: (value: number | null | undefined) => (
        <Tag color="orange">{(value ?? 0).toFixed(2)} 核时</Tag>
      ),
    },
    {
      title: '可用配额',
      key: 'balance',
      render: (record: SubAccount) => {
        const balance = record.balance_cpu_hours ?? 0;
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
          {active ? '活跃' : '禁用'}
        </Tag>
      ),
    },
    {
      title: '操作',
      key: 'action',
      render: (record: SubAccount) => (
        <Space size="small">
          <Button
            type="text"
            size="small"
            icon={<EditOutlined />}
            onClick={() => handleEditClick(record)}
          />
          <Popconfirm
            title="删除子账号"
            description="确定要删除这个子账号吗？"
            onConfirm={() => handleDeleteSubAccount(record.id)}
            okText="确定"
            cancelText="取消"
          >
            <Button
              type="text"
              size="small"
              danger
              icon={<DeleteOutlined />}
            />
          </Popconfirm>
          <Button
            type="text"
            size="small"
            icon={<EyeOutlined />}
            onClick={() => handleViewJobs(record)}
            title="查看任务"
          />
        </Space>
      ),
    },
  ];

  // 计算统计数据
  const totalPersonalQuota = subAccounts.reduce((sum, account) => sum + (account.personal_quota || 0), 0);
  const totalPersonalUsed = subAccounts.reduce((sum, account) => sum + (account.personal_used || 0), 0);
  const totalBalance = subAccounts.reduce((sum, account) => sum + (account.balance_cpu_hours || 0), 0);
  const activeCount = subAccounts.filter(account => account.is_active).length;

  if (loading) {
    return (
      <div style={{ textAlign: 'center', padding: '24px' }}>
        <Spin size="large" />
      </div>
    );
  }

  return (
    <div>
      {/* 页面标题 - 仅在非紧凑模式显示 */}
      {!compact && (
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
              icon={<UserAddOutlined />}
              onClick={() => setAddExistingModalVisible(true)}
            >
              添加现有用户
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
      )}

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
        />
      )}

      {/* 子账号表格 */}
      <Card
        title="子账号列表"
        extra={
          <Space>
            <Button icon={<ReloadOutlined />} onClick={loadSubAccounts}>
              刷新
            </Button>
            <Button
              icon={<UserAddOutlined />}
              onClick={() => setAddExistingModalVisible(true)}
            >
              添加现有用户
            </Button>
            <Button
              type="primary"
              icon={<PlusOutlined />}
              onClick={() => setCreateModalVisible(true)}
            >
              创建子账号
            </Button>
          </Space>
        }
      >
        <Table
          columns={columns}
          dataSource={subAccounts}
          rowKey="id"
          pagination={{ pageSize: 10 }}
          loading={loading}
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

      {/* 编辑子账号对话框 */}
      <Modal
        title={`编辑子账号 - ${editingSubAccount?.username}`}
        open={editModalVisible}
        onOk={() => editForm.submit()}
        onCancel={() => {
          setEditModalVisible(false);
          editForm.resetFields();
          setEditingSubAccount(null);
        }}
        width={600}
      >
        <Form
          form={editForm}
          layout="vertical"
          onFinish={handleUpdateSubAccount}
        >
          <Form.Item label="用户名">
            <Text strong>{editingSubAccount?.username}</Text>
          </Form.Item>

          <Form.Item label="邮箱">
            <Text>{editingSubAccount?.email}</Text>
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

          <Form.Item
            label="账号状态"
            name="is_active"
            valuePropName="checked"
            extra="禁用后子账号将无法登录和使用"
          >
            <Switch checkedChildren="启用" unCheckedChildren="禁用" />
          </Form.Item>

          <Divider />

          <Alert
            message="编辑说明"
            description={
              <ul style={{ margin: 0, paddingLeft: 20 }}>
                <li>用户名和邮箱不可修改</li>
                <li>可以修改个人配额限制</li>
                <li>可以启用或禁用子账号</li>
                <li>修改将立即生效</li>
              </ul>
            }
            type="info"
            showIcon
          />
        </Form>
      </Modal>

      {/* 添加现有用户对话框 */}
      <Modal
        title="添加现有用户为子账号"
        open={addExistingModalVisible}
        onOk={() => addExistingForm.submit()}
        onCancel={() => {
          setAddExistingModalVisible(false);
          addExistingForm.resetFields();
        }}
        width={500}
      >
        <Form
          form={addExistingForm}
          layout="vertical"
          onFinish={handleAddExistingUser}
        >
          <Form.Item
            label="用户名或邮箱"
            name="username_or_email"
            rules={[{ required: true, message: '请输入要添加的用户名或邮箱' }]}
            extra="输入已注册用户的用户名或邮箱地址"
          >
            <Input
              prefix={<SearchOutlined />}
              placeholder="请输入用户名或邮箱"
            />
          </Form.Item>

          <Form.Item
            label="分配配额（可选）"
            name="allocated_quota"
            extra="为该子账号分配的初始配额"
          >
            <InputNumber
              style={{ width: '100%' }}
              placeholder="配额（核时）"
              min={0}
              step={1}
              precision={2}
            />
          </Form.Item>

          <Alert
            message="添加说明"
            description={
              <ul style={{ margin: 0, paddingLeft: 20 }}>
                <li>只能添加账户类型为"个人用户"的现有账号</li>
                <li>添加后该用户将成为您的子账号</li>
                <li>子账号可以使用主账号分配的配额</li>
                <li>子账号也可以自己充值获得额外配额</li>
              </ul>
            }
            type="info"
            showIcon
          />
        </Form>
      </Modal>

      {/* 查看子账号任务对话框 */}
      <Modal
        title={`${viewingSubAccount?.username || ''} 的任务列表`}
        open={jobsModalVisible}
        onCancel={() => {
          setJobsModalVisible(false);
          setViewingSubAccount(null);
          setSubAccountJobs([]);
        }}
        footer={null}
        width={900}
      >
        <Spin spinning={jobsLoading}>
          <div style={{ marginBottom: 16 }}>
            <Text type="secondary">共 {jobsTotal} 个任务</Text>
          </div>
          <Table
            columns={[
              {
                title: '任务ID',
                dataIndex: 'id',
                key: 'id',
                width: 80,
              },
              {
                title: '状态',
                dataIndex: 'status',
                key: 'status',
                width: 100,
                render: (status: string) => {
                  const colorMap: Record<string, string> = {
                    COMPLETED: 'green',
                    RUNNING: 'blue',
                    QUEUED: 'orange',
                    FAILED: 'red',
                    CANCELLED: 'default',
                  };
                  return <Tag color={colorMap[status] || 'default'}>{status}</Tag>;
                },
              },
              {
                title: '创建时间',
                dataIndex: 'created_at',
                key: 'created_at',
                width: 180,
                render: (val: string) => val ? new Date(val).toLocaleString('zh-CN') : '-',
              },
              {
                title: '完成时间',
                dataIndex: 'finished_at',
                key: 'finished_at',
                width: 180,
                render: (val: string) => val ? new Date(val).toLocaleString('zh-CN') : '-',
              },
              {
                title: '消耗核时',
                dataIndex: 'actual_cpu_hours',
                key: 'actual_cpu_hours',
                width: 100,
                render: (val: number) => val ? val.toFixed(2) : '-',
              },
              {
                title: '操作',
                key: 'action',
                width: 80,
                render: (record: SubAccountJob) => (
                  <Button
                    type="link"
                    size="small"
                    onClick={() => {
                      setJobsModalVisible(false);
                      navigate(`/workspace/md-jobs/${record.id}`);
                    }}
                  >
                    详情
                  </Button>
                ),
              },
            ]}
            dataSource={subAccountJobs}
            rowKey="id"
            pagination={{ pageSize: 10 }}
            size="small"
          />
        </Spin>
      </Modal>
    </div>
  );
};

export default SubAccountManagementPanel;
