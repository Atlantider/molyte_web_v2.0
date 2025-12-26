/**
 * User Management Page
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
  Select,
  InputNumber,
  message,
  Popconfirm,
  Switch,
  Row,
  Col,
  Statistic,
  Typography,
  theme,
  Tabs,
} from 'antd';

const { Title, Text } = Typography;
import {
  PlusOutlined,
  DeleteOutlined,
  UserOutlined,
  SearchOutlined,
  FilterOutlined,
  ReloadOutlined,
  TeamOutlined,
  CrownOutlined,
  UserSwitchOutlined,
  CheckCircleOutlined,
  DollarOutlined,
  SettingOutlined,
  AppstoreOutlined,
  DatabaseOutlined,
  ThunderboltOutlined,
} from '@ant-design/icons';
import { useNavigate } from 'react-router-dom';
import AdminNav from '../../components/AdminNav';
import { useThemeStore } from '../../stores/themeStore';
import {
  getAllUsers,
  createUser,
  updateUser,
  deleteUser,
  updateUserStatus,
  UserListItem,
  UserCreate,
  UserUpdate,
  getAllPartitions,
  PartitionInfo,
} from '../../api/admin';
import { adminGrantCpuHours } from '../../api/billing';
import apiClient from '../../api/client';

const UserManagement: React.FC = () => {
  const navigate = useNavigate();
  const { token } = theme.useToken();
  const { mode } = useThemeStore();
  const isDark = mode === 'dark';
  const [loading, setLoading] = useState(false);
  const [users, setUsers] = useState<UserListItem[]>([]);
  const [filteredUsers, setFilteredUsers] = useState<UserListItem[]>([]);
  const [modalVisible, setModalVisible] = useState(false);
  const [editingUser, setEditingUser] = useState<UserListItem | null>(null);
  const [partitions, setPartitions] = useState<PartitionInfo[]>([]);
  const [form] = Form.useForm();

  // 定价配置状态（用于显示角色定价规则）
  const [pricingConfig, setPricingConfig] = useState<any>({
    global_price: 0.1,
    role_prices: { ADMIN: 0.1, PREMIUM: 0.08, USER: 0.1, GUEST: 0.15 },
  });

  // 批量操作状态
  const [selectedRowKeys, setSelectedRowKeys] = useState<number[]>([]);
  const [batchModalVisible, setBatchModalVisible] = useState(false);
  const [batchForm] = Form.useForm();

  // 赠送核时状态
  const [grantModalVisible, setGrantModalVisible] = useState(false);
  const [grantingUser, setGrantingUser] = useState<UserListItem | null>(null);
  const [grantForm] = Form.useForm();

  // 筛选和搜索状态
  const [searchText, setSearchText] = useState('');
  const [roleFilter, setRoleFilter] = useState<string | undefined>(undefined);
  const [userTypeFilter, setUserTypeFilter] = useState<string | undefined>(undefined);
  const [statusFilter, setStatusFilter] = useState<boolean | undefined>(undefined);
  const [organizationFilter, setOrganizationFilter] = useState('');

  useEffect(() => {
    loadUsers();
    loadPartitions();
    loadPricingConfig();
  }, []);

  // 应用筛选
  useEffect(() => {
    applyFilters();
  }, [users, searchText, roleFilter, userTypeFilter, statusFilter, organizationFilter]);

  const loadUsers = async () => {
    setLoading(true);
    try {
      const data = await getAllUsers();
      setUsers(data);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '加载用户列表失败');
    } finally {
      setLoading(false);
    }
  };

  const loadPricingConfig = async () => {
    try {
      const response = await apiClient.get('/billing/admin/pricing-config');
      setPricingConfig(response.data);
    } catch (error: any) {
      console.warn('加载定价配置失败:', error);
    }
  };



  // 批量操作处理
  const handleBatchOperation = () => {
    if (selectedRowKeys.length === 0) {
      message.warning('请先选择要操作的用户');
      return;
    }
    setBatchModalVisible(true);
  };

  const handleBatchSubmit = async () => {
    try {
      const values = batchForm.getFieldsValue();
      const updates: any = {};

      // 处理模块权限
      if (values.modules !== undefined && values.modules && values.modules.length > 0) {
        updates.allowed_modules = values.modules;
      }

      // 处理队列权限
      if (values.partitions !== undefined && values.partitions && values.partitions.length > 0) {
        updates.allowed_partitions = values.partitions;
      }

      // 处理计费方式
      if (values.billing_mode) {
        if (values.billing_mode === 'custom') {
          if (values.custom_price === undefined || values.custom_price === null) {
            message.error('自定义计费模式下，请输入价格');
            return;
          }
          updates.custom_cpu_hour_price = values.custom_price;
        } else if (values.billing_mode === 'role') {
          updates.custom_cpu_hour_price = null;
        }
      }

      // 检查是否有任何更新
      if (Object.keys(updates).length === 0) {
        message.warning('请至少选择一项要修改的内容');
        return;
      }

      // 批量更新选中的用户
      for (const userId of selectedRowKeys) {
        await updateUser(userId, updates);
      }

      message.success(`已成功更新 ${selectedRowKeys.length} 个用户`);
      setBatchModalVisible(false);
      batchForm.resetFields();
      setSelectedRowKeys([]);
      loadUsers();
    } catch (error: any) {
      message.error('批量更新失败: ' + (error.response?.data?.detail || '请稍后重试'));
    }
  };



  const applyFilters = () => {
    let filtered = [...users];

    // 搜索过滤
    if (searchText) {
      const search = searchText.toLowerCase();
      filtered = filtered.filter(
        (user) =>
          user.username.toLowerCase().includes(search) ||
          user.email.toLowerCase().includes(search) ||
          (user.organization && user.organization.toLowerCase().includes(search))
      );
    }

    // 角色过滤
    if (roleFilter) {
      filtered = filtered.filter((user) => user.role === roleFilter);
    }

    // 用户类型过滤
    if (userTypeFilter) {
      filtered = filtered.filter((user) => user.user_type === userTypeFilter);
    }

    // 状态过滤
    if (statusFilter !== undefined) {
      filtered = filtered.filter((user) => user.is_active === statusFilter);
    }

    // 组织过滤
    if (organizationFilter) {
      const org = organizationFilter.toLowerCase();
      filtered = filtered.filter(
        (user) => user.organization && user.organization.toLowerCase().includes(org)
      );
    }

    setFilteredUsers(filtered);
  };

  const handleResetFilters = () => {
    setSearchText('');
    setRoleFilter(undefined);
    setUserTypeFilter(undefined);
    setStatusFilter(undefined);
    setOrganizationFilter('');
  };

  const loadPartitions = async () => {
    try {
      const data = await getAllPartitions();
      setPartitions(data);
    } catch (error: any) {
      console.error('Failed to load partitions:', error);
      // 如果获取失败，使用默认值
      setPartitions([]);
    }
  };

  const handleCreate = () => {
    setEditingUser(null);
    form.resetFields();
    setModalVisible(true);
  };

  const handleEdit = (user: UserListItem) => {
    setEditingUser(user);
    const billingModeType = user.billing_mode || 'CORE_HOUR';
    const isCustomCoreHour = user.custom_cpu_hour_price && user.custom_cpu_hour_price > 0;
    form.setFieldsValue({
      ...user,
      billing_mode_type: billingModeType,
      billing_mode: isCustomCoreHour ? 'custom' : 'role',
      custom_price: isCustomCoreHour ? user.custom_cpu_hour_price : undefined,
    });
    setModalVisible(true);
  };

  const handleSubmit = async () => {
    try {
      const values = await form.validateFields();

      // 处理计费方式
      const updateData: any = { ...values };

      // 设置 billing_mode (CORE_HOUR or TASK_TYPE)
      updateData.billing_mode = values.billing_mode_type;

      // 根据计费模式处理价格
      if (values.billing_mode_type === 'CORE_HOUR') {
        if (values.billing_mode === 'custom') {
          if (values.custom_price === undefined || values.custom_price === null) {
            message.error('请输入自定义价格');
            return;
          }
          updateData.custom_cpu_hour_price = values.custom_price;
        } else {
          // 按用户类型计费时，清除自定义价格
          updateData.custom_cpu_hour_price = null;
        }
        // 清除任务价格自定义
        updateData.custom_task_prices = null;
      } else {
        // 任务类型计费
        updateData.custom_cpu_hour_price = null;
        // custom_task_prices 暂时保持原样，后续可单独配置
      }

      // 删除临时表单字段
      delete updateData.billing_mode_type;
      delete updateData.custom_price;

      if (editingUser) {
        // Update user
        await updateUser(editingUser.id, updateData as UserUpdate);
        message.success('用户更新成功');
      } else {
        // Create user
        await createUser(updateData as UserCreate);
        message.success('用户创建成功');
      }

      setModalVisible(false);
      loadUsers();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '操作失败');
    }
  };

  const handleDelete = async (userId: number) => {
    try {
      await deleteUser(userId);
      message.success('用户删除成功');
      loadUsers();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '删除失败');
    }
  };

  const handleStatusChange = async (userId: number, isActive: boolean) => {
    try {
      await updateUserStatus(userId, isActive);
      message.success(`用户已${isActive ? '启用' : '禁用'}`);
      loadUsers();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '操作失败');
    }
  };

  // 打开赠送核时对话框
  const handleOpenGrantModal = (user: UserListItem) => {
    setGrantingUser(user);
    grantForm.resetFields();
    setGrantModalVisible(true);
  };

  // 提交赠送核时
  const handleGrantCpuHours = async (values: any) => {
    if (!grantingUser) return;

    try {
      await adminGrantCpuHours(grantingUser.id, values.amount, values.reason);
      message.success(`成功赠送 ${values.amount} 核时给用户 ${grantingUser.username}`);
      setGrantModalVisible(false);
      grantForm.resetFields();
      setGrantingUser(null);
      loadUsers();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '赠送失败');
    }
  };

  const columns = [
    {
      title: 'ID',
      dataIndex: 'id',
      key: 'id',
      width: 60,
      fixed: 'left' as const,
      sorter: (a: UserListItem, b: UserListItem) => a.id - b.id,
      render: (id: number) => <span style={{ whiteSpace: 'nowrap' }}>{id}</span>,
    },
    {
      title: '用户名',
      dataIndex: 'username',
      key: 'username',
      width: 120,
      fixed: 'left' as const,
      sorter: (a: UserListItem, b: UserListItem) => a.username.localeCompare(b.username),
      render: (text: string, record: UserListItem) => (
        <a onClick={() => navigate(`/workspace/admin/users/${record.id}`)}>
          <UserOutlined /> {text}
        </a>
      ),
    },
    {
      title: '邮箱',
      dataIndex: 'email',
      key: 'email',
      width: 180,
      ellipsis: true,
      sorter: (a: UserListItem, b: UserListItem) => a.email.localeCompare(b.email),
    },
    {
      title: '角色',
      dataIndex: 'role',
      key: 'role',
      width: 100,
      sorter: (a: UserListItem, b: UserListItem) => a.role.localeCompare(b.role),
      render: (role: string) => {
        // 统一的颜色方案
        const roleConfig: any = {
          ADMIN: { color: '#85a5ff', label: '管理员', icon: <CrownOutlined /> },
          PREMIUM: { color: '#95de64', label: '高级用户', icon: <UserSwitchOutlined /> },
          USER: { color: '#b37feb', label: '普通用户', icon: <UserOutlined /> },
          GUEST: { color: '#ffc069', label: '访客', icon: <TeamOutlined /> },
        };
        const config = roleConfig[role] || { color: '#d9d9d9', label: role, icon: <UserOutlined /> };
        return (
          <Tag
            color={config.color}
            icon={config.icon}
            style={{
              fontSize: '12px',
              fontWeight: 500,
              padding: '2px 8px',
              borderRadius: '6px',
            }}
          >
            {config.label}
          </Tag>
        );
      },
    },
    {
      title: '状态',
      dataIndex: 'is_active',
      key: 'is_active',
      width: 60,
      sorter: (a: UserListItem, b: UserListItem) => Number(a.is_active) - Number(b.is_active),
      render: (isActive: boolean, record: UserListItem) => (
        <Switch
          checked={isActive}
          onChange={(checked) => handleStatusChange(record.id, checked)}
          size="small"
          style={{ minWidth: 'auto' }}
        />
      ),
    },
    {
      title: '可用余额',
      dataIndex: 'balance_cpu_hours',
      key: 'balance_cpu_hours',
      width: 100,
      sorter: (a: UserListItem, b: UserListItem) => a.balance_cpu_hours - b.balance_cpu_hours,
      render: (hours: number) => {
        const color = hours < 0 ? '#f5222d' : hours === 0 ? '#faad14' : '#52c41a';
        return <span style={{ color, fontWeight: 'bold' }}>{hours.toFixed(1)}h</span>;
      },
    },
    {
      title: '日/并发',
      key: 'limits',
      width: 80,
      render: (_: any, record: UserListItem) => `${record.daily_job_limit}/${record.concurrent_job_limit}`,
    },
    {
      title: '类型',
      dataIndex: 'user_type',
      key: 'user_type',
      width: 90,
      sorter: (a: UserListItem, b: UserListItem) => (a.user_type || '').localeCompare(b.user_type || ''),
      render: (type: string) => {
        const typeConfig: any = {
          STUDENT: { text: '学生', color: '#722ed1', icon: <UserOutlined /> },
          RESEARCHER: { text: '研究者', color: '#13c2c2', icon: <DatabaseOutlined /> },
          COMPANY: { text: '企业', color: '#eb2f96', icon: <TeamOutlined /> },
        };
        const config = typeConfig[type] || { text: type || '未设置', color: '#d9d9d9', icon: <UserOutlined /> };
        return (
          <Tag
            color={config.color}
            icon={config.icon}
            style={{
              fontSize: '12px',
              fontWeight: 500,
              padding: '2px 8px',
              borderRadius: '6px',
            }}
          >
            {config.text}
          </Tag>
        );
      },
    },
    {
      title: '组织',
      dataIndex: 'organization',
      key: 'organization',
      width: 120,
      ellipsis: true,
    },
    {
      title: '队列',
      dataIndex: 'allowed_partitions',
      key: 'allowed_partitions',
      width: 80,
      render: (partitions: string[] | null) => {
        if (!partitions || partitions.length === 0) {
          return <Tag color="default">全部</Tag>;
        }
        if (partitions.length === 1) {
          return <Tag color="blue">{partitions[0]}</Tag>;
        }
        return <Tag color="blue">{partitions.length}个</Tag>;
      },
    },
    {
      title: '模块权限',
      dataIndex: 'allowed_modules',
      key: 'allowed_modules',
      width: 100,
      render: (modules: string[] | null) => {
        if (!modules || modules.length === 0) {
          return <Tag color="default">全部</Tag>;
        }
        if (modules.length === 6) {
          return <Tag color="default">全部</Tag>;
        }
        return <Tag color="blue">{modules.length}个</Tag>;
      },
    },

    {
      title: '计费方式',
      key: 'pricing',
      width: 140,
      render: (_: any, record: UserListItem) => {
        const billingMode = record.billing_mode || 'CORE_HOUR';
        if (billingMode === 'TASK_TYPE') {
          return (
            <Tag color="purple" style={{ margin: 0, fontSize: '12px' }}>
              按任务
              {record.custom_task_prices && Object.keys(record.custom_task_prices).length > 0 && (
                <span style={{ marginLeft: 4, opacity: 0.7 }}>(自定义)</span>
              )}
            </Tag>
          );
        }
        // CORE_HOUR mode
        return (
          <>
            {record.custom_cpu_hour_price ? (
              <Tag color="orange" style={{ margin: 0, fontSize: '12px' }}>
                核时(¥{record.custom_cpu_hour_price})
              </Tag>
            ) : (
              <Tag color="blue" style={{ margin: 0, fontSize: '12px' }}>核时(标准)</Tag>
            )}
          </>
        );
      },
    },
    {
      title: '创建时间',
      dataIndex: 'created_at',
      key: 'created_at',
      width: 160,
      render: (time: string) => time ? new Date(time).toLocaleString('zh-CN') : '-',
    },
    {
      title: '操作',
      key: 'actions',
      width: 180,
      fixed: 'right' as const,
      render: (_: any, record: UserListItem) => (
        <Space size={4}>
          <Button
            type="link"
            size="small"
            style={{ padding: '0 4px' }}
            onClick={() => handleEdit(record)}
          >
            编辑
          </Button>
          <Button
            type="link"
            size="small"
            style={{ padding: '0 4px', color: token.colorSuccess }}
            onClick={() => handleOpenGrantModal(record)}
            icon={<DollarOutlined />}
          >
            赠送
          </Button>
          <Popconfirm
            title="确定要删除这个用户吗？"
            onConfirm={() => handleDelete(record.id)}
            okText="确定"
            cancelText="取消"
          >
            <Button type="link" size="small" danger style={{ padding: '0 4px' }}>
              删除
            </Button>
          </Popconfirm>
        </Space>
      ),
    },
  ];

  // 统计数据
  const totalUsers = users.length;
  const activeUsers = users.filter((u) => u.is_active).length;
  const adminUsers = users.filter((u) => u.role === 'ADMIN').length;
  const premiumUsers = users.filter((u) => u.role === 'PREMIUM').length;

  return (
    <div style={{ padding: '20px 24px', background: token.colorBgLayout, minHeight: 'calc(100vh - 64px)' }}>
      {/* 页面标题 */}
      <div style={{ marginBottom: 16 }}>
        <Title level={3} style={{ margin: 0, marginBottom: 4 }}>
          <TeamOutlined style={{ marginRight: 10, color: token.colorPrimary }} />
          用户管理
        </Title>
        <Text type="secondary">管理系统用户、角色配置和权限分配</Text>
      </div>

      <AdminNav />

      {/* 统计卡片 - 适中尺寸 */}
      <Row gutter={16} style={{ marginBottom: 20 }}>
        {[
          { label: '总用户数', value: totalUsers, color: '#85a5ff', icon: <TeamOutlined /> },
          { label: '活跃用户', value: activeUsers, color: '#95de64', icon: <CheckCircleOutlined /> },
          { label: '管理员', value: adminUsers, color: '#85a5ff', icon: <CrownOutlined /> },
          { label: '高级用户', value: premiumUsers, color: '#95de64', icon: <CrownOutlined /> },
        ].map((item, idx) => (
          <Col xs={12} sm={6} key={idx}>
            <div style={{
              padding: '16px 20px',
              background: isDark ? 'rgba(255,255,255,0.03)' : 'white',
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
                <Text strong style={{ fontSize: 22, color: item.color }}>{item.value}</Text>
              </div>
            </div>
          </Col>
        ))}
      </Row>

      <Card
        title={
          <Space size={10}>
            <UserOutlined style={{ color: token.colorPrimary }} />
            <span>用户列表</span>
            {selectedRowKeys.length > 0 && (
              <Tag color="blue" style={{ fontSize: 13 }}>{selectedRowKeys.length} 已选中</Tag>
            )}
          </Space>
        }
        bordered={false}
        extra={
          <Space size={10}>
            <Button
              icon={<SettingOutlined />}
              onClick={handleBatchOperation}
              disabled={selectedRowKeys.length === 0}
              title={selectedRowKeys.length === 0 ? '请先选择用户' : `已选择 ${selectedRowKeys.length} 个用户`}
            >
              批量设置 ({selectedRowKeys.length})
            </Button>
            <Button type="primary" icon={<PlusOutlined />} onClick={handleCreate}>
              创建用户
            </Button>
          </Space>
        }
        style={{
          borderRadius: 12,
          boxShadow: isDark ? '0 2px 12px rgba(0,0,0,0.3)' : '0 2px 12px rgba(0,0,0,0.08)',
          background: token.colorBgContainer,
        }}
        styles={{ body: { padding: '16px 20px' } }}
      >
        {/* 筛选栏 */}
        <div style={{
          marginBottom: 16,
          padding: '16px',
          background: isDark ? 'rgba(255,255,255,0.02)' : '#fafafa',
          borderRadius: 10,
          border: `1px solid ${token.colorBorder}`,
        }}>
          {/* 第一行：搜索框 + 操作按钮 */}
          <div style={{ display: 'flex', gap: 12, marginBottom: 14, alignItems: 'center' }}>
            <Input
              placeholder="搜索用户名、邮箱、组织..."
              prefix={<SearchOutlined style={{ color: token.colorTextSecondary }} />}
              value={searchText}
              onChange={(e) => setSearchText(e.target.value)}
              allowClear
              style={{ flex: 1, maxWidth: 300 }}
            />
            <Button icon={<FilterOutlined />} onClick={handleResetFilters}>重置筛选</Button>
            <Button icon={<ReloadOutlined />} onClick={loadUsers}>刷新</Button>
            <Text type="secondary" style={{ marginLeft: 'auto', fontSize: 13 }}>
              显示 {filteredUsers.length} / {users.length} 用户
            </Text>
          </div>

          {/* 第二行：筛选标签 */}
          <div style={{ display: 'flex', gap: 24, flexWrap: 'wrap', alignItems: 'center' }}>
            {/* 角色筛选 */}
            <div style={{ display: 'flex', alignItems: 'center', gap: 10 }}>
              <Text type="secondary" style={{ fontSize: 13 }}>角色:</Text>
              <Space size={6}>
                {[
                  { value: 'ADMIN', label: '管理员', color: '#85a5ff' },
                  { value: 'PREMIUM', label: '高级', color: '#95de64' },
                  { value: 'USER', label: '普通', color: '#b37feb' },
                  { value: 'GUEST', label: '访客', color: '#ffc069' },
                ].map((role) => {
                  const isSelected = roleFilter === role.value;
                  return (
                    <div
                      key={role.value}
                      onClick={() => setRoleFilter(isSelected ? undefined : role.value)}
                      style={{
                        cursor: 'pointer',
                        fontSize: 13,
                        padding: '4px 12px',
                        borderRadius: '4px',
                        border: `1px solid ${role.color}`,
                        color: isSelected ? 'white' : role.color,
                        background: isSelected ? role.color : 'transparent',
                        transition: 'all 0.2s',
                        whiteSpace: 'nowrap',
                      }}
                    >
                      {role.label} ({users.filter(u => u.role === role.value).length})
                    </div>
                  );
                })}
              </Space>
            </div>

            {/* 类型筛选 */}
            <div style={{ display: 'flex', alignItems: 'center', gap: 10 }}>
              <Text type="secondary" style={{ fontSize: 13 }}>类型:</Text>
              <Space size={6}>
                {[
                  { value: 'STUDENT', label: '学生', color: '#722ed1' },
                  { value: 'RESEARCHER', label: '研究者', color: '#13c2c2' },
                  { value: 'COMPANY', label: '企业', color: '#eb2f96' },
                ].map((type) => {
                  const isSelected = userTypeFilter === type.value;
                  return (
                    <div
                      key={type.value}
                      onClick={() => setUserTypeFilter(isSelected ? undefined : type.value)}
                      style={{
                        cursor: 'pointer',
                        fontSize: 13,
                        padding: '4px 12px',
                        borderRadius: '4px',
                        border: `1px solid ${type.color}`,
                        color: isSelected ? 'white' : type.color,
                        background: isSelected ? type.color : 'transparent',
                        transition: 'all 0.2s',
                        whiteSpace: 'nowrap',
                      }}
                    >
                      {type.label}
                    </div>
                  );
                })}
              </Space>
            </div>

            {/* 状态筛选 */}
            <div style={{ display: 'flex', alignItems: 'center', gap: 10 }}>
              <Text type="secondary" style={{ fontSize: 13 }}>状态:</Text>
              <Space size={6}>
                <div
                  onClick={() => setStatusFilter(statusFilter === true ? undefined : true)}
                  style={{
                    cursor: 'pointer',
                    fontSize: 13,
                    padding: '4px 12px',
                    borderRadius: '4px',
                    border: `1px solid #52c41a`,
                    color: statusFilter === true ? 'white' : '#52c41a',
                    background: statusFilter === true ? '#52c41a' : 'transparent',
                    transition: 'all 0.2s',
                    whiteSpace: 'nowrap',
                  }}
                >
                  激活 ({users.filter(u => u.is_active).length})
                </div>
                <div
                  onClick={() => setStatusFilter(statusFilter === false ? undefined : false)}
                  style={{
                    cursor: 'pointer',
                    fontSize: 13,
                    padding: '4px 12px',
                    borderRadius: '4px',
                    border: `1px solid #f5222d`,
                    color: statusFilter === false ? 'white' : '#f5222d',
                    background: statusFilter === false ? '#f5222d' : 'transparent',
                    transition: 'all 0.2s',
                    whiteSpace: 'nowrap',
                  }}
                >
                  禁用 ({users.filter(u => !u.is_active).length})
                </div>
              </Space>
            </div>
          </div>
        </div>
        <Table
          dataSource={filteredUsers}
          columns={columns}
          rowKey="id"
          loading={loading}
          scroll={{ x: 1200 }}
          rowSelection={{
            selectedRowKeys,
            onChange: (keys) => {
              // 只保留在 filteredUsers 中存在的 key
              const validKeys = keys.filter(key =>
                filteredUsers.some(user => user.id === key)
              );
              setSelectedRowKeys(validKeys as number[]);
            },
            selections: [
              Table.SELECTION_ALL,
              Table.SELECTION_INVERT,
              Table.SELECTION_NONE,
            ],
          }}
          pagination={{
            pageSize: 20,
            showSizeChanger: true,
            showTotal: (total) => `共 ${total} 个用户`,
            pageSizeOptions: ['10', '20', '50', '100'],
            style: { marginTop: 16 },
          }}
          style={{
            borderRadius: 10,
          }}
          className="user-management-table"
        />
      </Card>

      {/* Create/Edit User Modal */}
      <Modal
        title={
          <Space>
            <UserOutlined style={{ color: token.colorPrimary }} />
            <span style={{ fontWeight: 600, color: token.colorText }}>
              {editingUser ? '编辑用户' : '创建用户'}
            </span>
          </Space>
        }
        open={modalVisible}
        onOk={handleSubmit}
        onCancel={() => setModalVisible(false)}
        width={600}
        centered
        okText="确定"
        cancelText="取消"
        styles={{
          body: {
            maxHeight: '70vh',
            overflowY: 'auto',
            padding: '24px',
            background: token.colorBgContainer,
          },
          header: {
            borderBottom: `1px solid ${token.colorBorder}`,
            paddingBottom: 16,
            background: token.colorBgContainer,
          },
          content: {
            borderRadius: 12,
            boxShadow: isDark ? '0 8px 32px rgba(0,0,0,0.4)' : '0 8px 32px rgba(0,0,0,0.12)',
          },
        }}
      >
        <Form form={form} layout="vertical">
          <Form.Item
            name="username"
            label="用户名"
            rules={[
              { required: true, message: '请输入用户名' },
              { min: 3, message: '用户名至少 3 个字符' },
            ]}
          >
            <Input placeholder="请输入用户名" disabled={!!editingUser} />
          </Form.Item>

          <Form.Item
            name="email"
            label="邮箱"
            rules={[
              { required: true, message: '请输入邮箱' },
              { type: 'email', message: '请输入有效的邮箱地址' },
            ]}
          >
            <Input placeholder="请输入邮箱" />
          </Form.Item>

          {!editingUser && (
            <Form.Item
              name="password"
              label="密码"
              rules={[
                { required: true, message: '请输入密码' },
                { min: 6, message: '密码至少 6 个字符' },
              ]}
            >
              <Input.Password placeholder="请输入密码" />
            </Form.Item>
          )}

          <Form.Item
            name="role"
            label="角色"
            rules={[{ required: true, message: '请选择角色' }]}
            initialValue="USER"
          >
            <Select>
              <Select.Option value="ADMIN">管理员</Select.Option>
              <Select.Option value="PREMIUM">高级用户</Select.Option>
              <Select.Option value="USER">普通用户</Select.Option>
              <Select.Option value="GUEST">访客</Select.Option>
            </Select>
          </Form.Item>

          <Form.Item
            name="balance_cpu_hours"
            label="初始可用余额 (核时)"
            rules={[{ required: true, message: '请输入初始可用余额' }]}
            initialValue={100}
          >
            <InputNumber min={0} style={{ width: '100%' }} />
          </Form.Item>

          <Form.Item
            name="free_cpu_hours_granted"
            label="初始赠送核时"
            rules={[{ required: true, message: '请输入初始赠送核时' }]}
            initialValue={100}
          >
            <InputNumber min={0} style={{ width: '100%' }} />
          </Form.Item>

          <Form.Item
            name="daily_job_limit"
            label="每日任务限制"
            rules={[{ required: true, message: '请输入每日任务限制' }]}
            initialValue={10}
          >
            <InputNumber min={0} style={{ width: '100%' }} />
          </Form.Item>

          <Form.Item
            name="concurrent_job_limit"
            label="并发任务限制"
            rules={[{ required: true, message: '请输入并发任务限制' }]}
            initialValue={3}
          >
            <InputNumber min={0} style={{ width: '100%' }} />
          </Form.Item>

          <Form.Item
            name="storage_quota_gb"
            label="存储配额 (GB)"
            rules={[{ required: true, message: '请输入存储配额' }]}
            initialValue={10}
          >
            <InputNumber min={0} style={{ width: '100%' }} />
          </Form.Item>

          <Form.Item
            name="allowed_partitions"
            label="可用队列"
            tooltip="选择用户可以使用的队列。管理员留空表示可以使用所有队列"
            initialValue={['cpu']}
          >
            <Select
              mode="multiple"
              placeholder="选择可用队列（留空表示全部队列）"
              allowClear
            >
              {partitions.length > 0 ? (
                partitions.map((p) => (
                  <Select.Option key={p.name} value={p.name}>
                    {p.name} ({p.state === 'up' ? '可用' : '不可用'})
                  </Select.Option>
                ))
              ) : (
                <>
                  <Select.Option value="cpu">cpu</Select.Option>
                  <Select.Option value="gpu">gpu</Select.Option>
                  <Select.Option value="debug">debug</Select.Option>
                </>
              )}
            </Select>
          </Form.Item>

          <Form.Item
            name="allowed_modules"
            label="可用模块"
            tooltip="选择用户可以访问的功能模块。留空表示可以访问所有模块"
            initialValue={['electrolytes', 'md', 'analysis', 'qc', 'ai-discovery', 'anion-generation']}
          >
            <Select
              mode="multiple"
              placeholder="选择可用模块（留空表示全部模块）"
              allowClear
            >
              <Select.Option value="electrolytes">溶液配方管理</Select.Option>
              <Select.Option value="md">溶液MD分析</Select.Option>
              <Select.Option value="analysis">溶鞘QC分析（后处理）</Select.Option>
              <Select.Option value="qc">溶元QC分析</Select.Option>
              <Select.Option value="ai-discovery">溶元AI推荐</Select.Option>
              <Select.Option value="anion-generation">溶盐FF开发</Select.Option>
            </Select>
          </Form.Item>

          {/* 计费模式设置 */}
          <Form.Item
            name="billing_mode_type"
            label={
              <Space size={4}>
                <DollarOutlined style={{ color: token.colorPrimary }} />
                <span>计费模式</span>
              </Space>
            }
            initialValue="CORE_HOUR"
          >
            <Select onChange={() => form.validateFields(['custom_price', 'billing_mode'])}>
              <Select.Option value="CORE_HOUR">按核时计费</Select.Option>
              <Select.Option value="TASK_TYPE">按任务类型计费</Select.Option>
            </Select>
          </Form.Item>

          {/* 核时计费选项 */}
          <Form.Item
            noStyle
            shouldUpdate
          >
            {({ getFieldValue }) => {
              const billingModeType = getFieldValue('billing_mode_type');
              if (billingModeType !== 'CORE_HOUR') return null;
              return (
                <>
                  <Form.Item
                    name="billing_mode"
                    label="核时定价方式"
                    initialValue="role"
                  >
                    <Select onChange={() => form.validateFields(['custom_price'])}>
                      <Select.Option value="role">按用户类型（标准）</Select.Option>
                      <Select.Option value="custom">自定义核时单价</Select.Option>
                    </Select>
                  </Form.Item>

                  {/* 显示用户类型定价规则 */}
                  <Form.Item noStyle shouldUpdate>
                    {({ getFieldValue: gfv }) => {
                      const billingMode = gfv('billing_mode');
                      const userRole = gfv('role');
                      if (billingMode === 'role' && userRole) {
                        const rolePrice = pricingConfig.role_prices?.[userRole] || pricingConfig.global_price;
                        return (
                          <div style={{
                            padding: '12px',
                            backgroundColor: isDark ? 'rgba(24, 144, 255, 0.1)' : '#f0f5ff',
                            borderRadius: '4px',
                            marginBottom: '16px',
                            fontSize: '12px',
                            color: isDark ? '#69b1ff' : '#0050b3'
                          }}>
                            <div style={{ marginBottom: '4px' }}>
                              当前用户类型 <strong>{userRole === 'ADMIN' ? '管理员' : userRole === 'PREMIUM' ? '高级' : userRole === 'USER' ? '普通' : '访客'}</strong> 定价：
                            </div>
                            <div style={{ fontSize: '14px', fontWeight: 'bold' }}>
                              ¥{rolePrice.toFixed(4)}/核时
                            </div>
                          </div>
                        );
                      }
                      return null;
                    }}
                  </Form.Item>

                  {/* 自定义核时价格输入 */}
                  <Form.Item noStyle shouldUpdate>
                    {({ getFieldValue: gfv }) => {
                      const billingMode = gfv('billing_mode');
                      return billingMode === 'custom' ? (
                        <Form.Item
                          name="custom_price"
                          label="自定义核时单价"
                          rules={[{ required: true, message: '请输入自定义价格' }]}
                        >
                          <InputNumber
                            min={0.01}
                            step={0.01}
                            precision={4}
                            style={{ width: '100%' }}
                            addonAfter="元/核时"
                            placeholder="请输入价格"
                          />
                        </Form.Item>
                      ) : null;
                    }}
                  </Form.Item>
                </>
              );
            }}
          </Form.Item>

          {/* 任务类型计费选项 */}
          <Form.Item
            noStyle
            shouldUpdate
          >
            {({ getFieldValue }) => {
              const billingModeType = getFieldValue('billing_mode_type');
              if (billingModeType !== 'TASK_TYPE') return null;
              return (
                <div style={{
                  padding: '16px',
                  backgroundColor: isDark ? 'rgba(114, 46, 209, 0.1)' : '#f9f0ff',
                  borderRadius: '8px',
                  marginBottom: '16px',
                }}>
                  <div style={{ marginBottom: '12px', fontWeight: 500, color: isDark ? '#b37feb' : '#722ed1' }}>
                    任务类型计费说明
                  </div>
                  <div style={{ fontSize: '12px', color: isDark ? '#d3adf7' : '#531dab' }}>
                    用户将按照任务类型的固定价格计费，而非核时。价格由系统全局设置，如需为此用户自定义，请在用户详情页单独配置。
                  </div>
                </div>
              );
            }}
          </Form.Item>
        </Form>
      </Modal>

      {/* 批量操作 Modal */}
      <Modal
        title={
          <Space>
            <SettingOutlined style={{ color: token.colorPrimary }} />
            <span style={{ fontWeight: 600, color: token.colorText }}>
              批量管理 ({selectedRowKeys.length} 个用户)
            </span>
          </Space>
        }
        open={batchModalVisible}
        onOk={handleBatchSubmit}
        onCancel={() => {
          setBatchModalVisible(false);
          batchForm.resetFields();
        }}
        width={650}
        centered
        okText="确定"
        cancelText="取消"
        styles={{
          body: {
            padding: '24px',
            background: token.colorBgContainer,
            maxHeight: '70vh',
            overflowY: 'auto',
          },
          header: {
            borderBottom: `1px solid ${token.colorBorder}`,
            paddingBottom: 16,
            background: token.colorBgContainer,
          },
          content: {
            borderRadius: 12,
            boxShadow: isDark ? '0 8px 32px rgba(0,0,0,0.4)' : '0 8px 32px rgba(0,0,0,0.12)',
          },
        }}
      >
        <Form form={batchForm} layout="vertical" style={{ marginTop: 24 }}>
          {/* 价格设置 */}
          <Form.Item
            name="billing_mode"
            label={
              <Space size={4}>
                <DollarOutlined style={{ color: token.colorPrimary }} />
                <span>价格设置</span>
              </Space>
            }
            tooltip="选择计费方式（不选择则不修改）"
          >
            <Select placeholder="选择计费方式（留空表示不修改）" allowClear onChange={() => batchForm.validateFields(['custom_price'])}>
              <Select.Option value="role">按角色计费</Select.Option>
              <Select.Option value="custom">自定义计费</Select.Option>
            </Select>
          </Form.Item>

          {/* 显示角色定价规则 */}
          <div style={{
            padding: '12px',
            backgroundColor: '#f0f5ff',
            borderRadius: '4px',
            marginBottom: '16px',
            fontSize: '12px',
            color: '#0050b3'
          }}>
            <div style={{ marginBottom: '8px', fontWeight: 500 }}>
              按角色定价规则：
            </div>
            {Object.entries(pricingConfig.role_prices || {}).map(([role, price]) => (
              <div key={role} style={{ marginBottom: '4px' }}>
                {role === 'ADMIN' ? '管理员' : role === 'PREMIUM' ? '高级' : role === 'USER' ? '普通' : '访客'}：
                <span style={{ fontWeight: 'bold', color: '#1890ff', marginLeft: '8px' }}>
                  ¥{(price as number).toFixed(4)}/核时
                </span>
              </div>
            ))}
          </div>

          <Form.Item
            noStyle
            shouldUpdate
          >
            {({ getFieldValue }) => {
              const billingMode = getFieldValue('billing_mode');
              return billingMode === 'custom' ? (
                <Form.Item
                  name="custom_price"
                  label="自定义价格"
                  rules={[{ required: true, message: '请输入自定义价格' }]}
                >
                  <InputNumber
                    min={0.01}
                    step={0.01}
                    precision={4}
                    style={{ width: '100%' }}
                    placeholder="请输入价格（元/核时）"
                  />
                </Form.Item>
              ) : null;
            }}
          </Form.Item>

          {/* 队列设置 */}
          <Form.Item
            name="partitions"
            label={
              <Space size={4}>
                <ThunderboltOutlined style={{ color: token.colorPrimary }} />
                <span>队列权限</span>
              </Space>
            }
            tooltip="选择用户可以使用的队列（不选择则不修改）"
          >
            <Select
              mode="multiple"
              placeholder="选择队列（留空表示不修改）"
              allowClear
            >
              {partitions.length > 0 ? (
                partitions.map((p) => (
                  <Select.Option key={p.name} value={p.name}>
                    {p.name} ({p.state === 'up' ? '可用' : '不可用'})
                  </Select.Option>
                ))
              ) : (
                <>
                  <Select.Option value="cpu">cpu</Select.Option>
                  <Select.Option value="gpu">gpu</Select.Option>
                  <Select.Option value="debug">debug</Select.Option>
                </>
              )}
            </Select>
          </Form.Item>

          {/* 模块权限 */}
          <Form.Item
            name="modules"
            label={
              <Space size={4}>
                <AppstoreOutlined style={{ color: token.colorPrimary }} />
                <span>模块权限</span>
              </Space>
            }
            tooltip="选择用户可以访问的功能模块（不选择则不修改）"
          >
            <Select
              mode="multiple"
              placeholder="选择模块（留空表示不修改）"
              allowClear
            >
              <Select.Option value="electrolytes">溶液配方管理</Select.Option>
              <Select.Option value="md">溶液MD分析</Select.Option>
              <Select.Option value="analysis">溶鞘QC分析（后处理）</Select.Option>
              <Select.Option value="qc">溶元QC分析</Select.Option>
              <Select.Option value="ai-discovery">溶元AI推荐</Select.Option>
              <Select.Option value="anion-generation">溶盐FF开发</Select.Option>
            </Select>
          </Form.Item>
        </Form>
      </Modal>

      {/* 赠送核时对话框 */}
      <Modal
        title={
          <Space>
            <DollarOutlined style={{ color: token.colorSuccess }} />
            <span style={{ fontWeight: 600, color: token.colorText }}>
              赠送核时 - {grantingUser?.username}
            </span>
          </Space>
        }
        open={grantModalVisible}
        onOk={() => grantForm.submit()}
        onCancel={() => {
          setGrantModalVisible(false);
          grantForm.resetFields();
          setGrantingUser(null);
        }}
        okText="确认赠送"
        cancelText="取消"
        width={500}
      >
        <Form
          form={grantForm}
          layout="vertical"
          onFinish={handleGrantCpuHours}
          style={{ marginTop: 16 }}
        >
          {grantingUser && (
            <div style={{
              padding: 16,
              background: isDark ? 'rgba(255,255,255,0.05)' : 'rgba(0,0,0,0.02)',
              borderRadius: 8,
              marginBottom: 16
            }}>
              <Row gutter={16}>
                <Col span={12}>
                  <div style={{ marginBottom: 8 }}>
                    <Text type="secondary">当前余额</Text>
                  </div>
                  <div style={{ fontSize: 20, fontWeight: 'bold', color: grantingUser.balance_cpu_hours < 0 ? token.colorError : token.colorSuccess }}>
                    {grantingUser.balance_cpu_hours.toFixed(2)} 核时
                  </div>
                </Col>
                <Col span={12}>
                  <div style={{ marginBottom: 8 }}>
                    <Text type="secondary">管理员已赠送</Text>
                  </div>
                  <div style={{ fontSize: 20, fontWeight: 'bold', color: token.colorPrimary }}>
                    {grantingUser.admin_granted_cpu_hours?.toFixed(2) || '0.00'} 核时
                  </div>
                </Col>
              </Row>
            </div>
          )}

          <Form.Item
            name="amount"
            label="赠送核时数量"
            rules={[
              { required: true, message: '请输入赠送核时数量' },
              { type: 'number', min: 0.01, message: '赠送数量必须大于 0' },
            ]}
          >
            <InputNumber
              style={{ width: '100%' }}
              placeholder="请输入赠送核时数量"
              min={0.01}
              step={100}
              precision={2}
              addonAfter="核时"
            />
          </Form.Item>

          <Form.Item
            name="reason"
            label="赠送原因"
            rules={[{ required: true, message: '请输入赠送原因' }]}
          >
            <Input.TextArea
              placeholder="请输入赠送原因，例如：补偿用户、活动奖励等"
              rows={3}
              maxLength={200}
              showCount
            />
          </Form.Item>

          <div style={{
            padding: 12,
            background: isDark ? 'rgba(24,144,255,0.1)' : 'rgba(24,144,255,0.05)',
            borderRadius: 6,
            border: `1px solid ${isDark ? 'rgba(24,144,255,0.2)' : 'rgba(24,144,255,0.1)'}`
          }}>
            <Text type="secondary" style={{ fontSize: 12 }}>
              💡 提示：赠送核时会同时增加用户的可用余额和管理员赠送统计，所有操作都会记录到交易日志中。
            </Text>
          </div>
        </Form>
      </Modal>

    </div>
  );
};

export default UserManagement;

