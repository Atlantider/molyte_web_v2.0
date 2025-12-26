import React, { useState, useEffect } from 'react';
import {
  Table,
  Button,
  Space,
  Modal,
  Form,
  Input,
  Select,
  message,
  Spin,
  Card,
  Statistic,
  Row,
  Col,
  Drawer,
  Tabs,
  Tag,
  Progress,
  Descriptions,
  InputNumber,
  Switch,
  Tooltip,
  Badge,
  Typography,
  theme,
} from 'antd';
import {
  PlusOutlined,
  EyeOutlined,
  EditOutlined,
  DeleteOutlined,
  UserOutlined,
  TeamOutlined,
  DollarOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  ReloadOutlined,
  SearchOutlined,
  BankOutlined,
} from '@ant-design/icons';
import { formatCpuHours, formatQuota, QUOTA_PRECISION } from '../../utils/formatQuotaDisplay';
import { updateUserPrice } from '../../api/billing';
import AdminNav from '../../components/AdminNav';
import {
  getAllMasterAccounts,
  getMasterAccountDetail,
  getSubAccounts,
  createMasterAccount,
  createSubAccount,
  updateSubAccount,
  deleteSubAccount,
  deleteUser,
  deleteMasterAccount,
  MasterAccount,
  SubAccount,
  CreateMasterAccountRequest,
  CreateSubAccountRequest,
} from '../../api/admin';

const { Title, Text } = Typography;

const MasterAccountManagement: React.FC = () => {
  const { token } = theme.useToken();
  const [masterAccounts, setMasterAccounts] = useState<MasterAccount[]>([]);
  const [filteredAccounts, setFilteredAccounts] = useState<MasterAccount[]>([]);
  const [subAccounts, setSubAccounts] = useState<SubAccount[]>([]);
  const [loading, setLoading] = useState(false);
  const [createModalVisible, setCreateModalVisible] = useState(false);
  const [editModalVisible, setEditModalVisible] = useState(false);
  const [createSubModalVisible, setCreateSubModalVisible] = useState(false);
  const [editSubModalVisible, setEditSubModalVisible] = useState(false);
  const [detailDrawerVisible, setDetailDrawerVisible] = useState(false);
  const [selectedMaster, setSelectedMaster] = useState<MasterAccount | null>(null);
  const [selectedSubAccount, setSelectedSubAccount] = useState<SubAccount | null>(null);
  const [searchText, setSearchText] = useState('');
  const [form] = Form.useForm();
  const [editForm] = Form.useForm();
  const [subAccountForm] = Form.useForm();
  const [editSubForm] = Form.useForm();

  // 使用 Context 兼容的 Hook
  const [modal, contextHolder] = Modal.useModal();
  const [messageApi, contextHolderMessage] = message.useMessage();

  // 新增状态
  const [createMode, setCreateMode] = useState<'existing' | 'new'>('existing'); // 创建模式
  const [allUsers, setAllUsers] = useState<any[]>([]); // 所有用户列表
  const [loadingUsers, setLoadingUsers] = useState(false);
  const [selectedRowKeys, setSelectedRowKeys] = useState<React.Key[]>([]); // 批量选择
  const [batchModalVisible, setBatchModalVisible] = useState(false);
  const [batchForm] = Form.useForm();

  // 加载主账号列表
  const loadMasterAccounts = async () => {
    setLoading(true);
    try {
      console.log('📡 开始加载主账号列表...');
      const response = await getAllMasterAccounts();
      console.log('✅ 主账号列表加载成功:', response);
      setMasterAccounts(response || []);
      setFilteredAccounts(response || []);
    } catch (error: any) {
      console.error('❌ 加载主账号列表失败:', error);
      if (error.response) {
        console.error('响应错误:', error.response.status, error.response.data);
        messageApi.error(`加载主账号列表失败: ${error.message}`);
      }
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    loadMasterAccounts();
  }, []);

  // 搜索过滤
  useEffect(() => {
    if (!searchText) {
      setFilteredAccounts(masterAccounts);
    } else {
      const filtered = masterAccounts.filter(
        (account) =>
          account.username?.toLowerCase().includes(searchText.toLowerCase()) ||
          account.email?.toLowerCase().includes(searchText.toLowerCase()) ||
          account.organization?.toLowerCase().includes(searchText.toLowerCase())
      );
      setFilteredAccounts(filtered);
    }
  }, [searchText, masterAccounts]);

  // 加载所有用户（用于下拉选择）
  const loadAllUsers = async () => {
    setLoadingUsers(true);
    try {
      const { getAllUsers } = await import('../../api/admin');
      const users = await getAllUsers();
      // 过滤掉已经是主账号或子账号的用户
      const masterUserIds = masterAccounts.map(m => m.user_id);
      const availableUsers = users.filter(
        (u: any) => !masterUserIds.includes(u.id) && u.account_type !== 'sub_account'
      );
      setAllUsers(availableUsers);
    } catch (error) {
      console.error('加载用户列表失败:', error);
    } finally {
      setLoadingUsers(false);
    }
  };

  // 加载子账号列表
  const loadSubAccounts = async (masterId: number) => {
    try {
      const response = await getSubAccounts(masterId);
      setSubAccounts(response || []);
    } catch (error) {
      messageApi.error('加载子账号列表失败');
    }
  };

  // 创建主账号
  const handleCreateMasterAccount = async (values: any) => {
    try {
      let requestData: CreateMasterAccountRequest;

      if (createMode === 'existing') {
        // 升级已有用户
        requestData = {
          user_id: values.user_id,
          max_sub_accounts: values.max_sub_accounts || 10,
          organization: values.organization,
        };
      } else {
        // 创建新用户
        requestData = {
          username: values.username,
          email: values.email,
          password: values.password,
          organization: values.organization,
          max_sub_accounts: values.max_sub_accounts || 10,
          initial_cpu_hours: values.initial_cpu_hours || 0,
        };
      }

      await createMasterAccount(requestData);
      await createMasterAccount(requestData);
      messageApi.success('主账号创建成功');
      setCreateModalVisible(false);
      form.resetFields();
      setCreateMode('existing');
      loadMasterAccounts();
    } catch (error: any) {
      messageApi.error(error.response?.data?.detail || '创建主账号失败');
    }
  };

  // 查看详情
  const handleViewDetail = async (master: MasterAccount) => {
    setSelectedMaster(master);
    setDetailDrawerVisible(true);
    await loadSubAccounts(master.id);
  };

  // 编辑主账号
  const handleEditMasterAccount = (master: MasterAccount) => {
    setSelectedMaster(master);
    editForm.setFieldsValue({
      balance_cpu_hours: master.balance_cpu_hours,
      max_sub_accounts: master.max_sub_accounts,
      is_active: master.is_active,
    });
    setEditModalVisible(true);
  };

  // 提交编辑
  const handleSubmitEdit = async (values: any) => {
    if (!selectedMaster) return;

    try {
      // TODO: 调用更新 API
      // TODO: 调用更新 API
      messageApi.success('主账号更新成功');
      setEditModalVisible(false);
      editForm.resetFields();
      loadMasterAccounts();
    } catch (error) {
      messageApi.error('更新主账号失败');
    }
  };

  // 删除主账号
  const handleDeleteMasterAccount = (master: MasterAccount) => {
    modal.confirm({
      title: '删除主账号',
      content: `确定要删除主账号 "${master.username}" 吗？\n\n此操作将：\n1. 将主账号转换为个人账号\n2. 将所有子账号转换为个人账号\n3. 删除主账号记录\n\n用户账号和数据不会被删除。`,
      okText: '确定',
      cancelText: '取消',
      okButtonProps: { danger: true },
      onOk: async () => {
        try {
          await deleteMasterAccount(master.id);
          messageApi.success('主账号已删除，已转换为个人账号');
          loadMasterAccounts();
        } catch (error: any) {
          console.error('删除主账号失败:', error);
          messageApi.error(`删除主账号失败: ${error.response?.data?.detail || error.message}`);
        }
      },
    });
  };

  // 创建子账号
  const handleCreateSubAccount = async (values: CreateSubAccountRequest) => {
    if (!selectedMaster) return;

    try {
      await createSubAccount(selectedMaster.id, values);
      messageApi.success('子账号创建成功');
      subAccountForm.resetFields();
      setCreateSubModalVisible(false);
      await loadSubAccounts(selectedMaster.id);
      await loadMasterAccounts(); // 刷新主账号列表
    } catch (error) {
      messageApi.error('创建子账号失败');
    }
  };

  // 编辑子账号
  const handleEditSubAccount = (subAccount: SubAccount) => {
    setSelectedSubAccount(subAccount);
    editSubForm.setFieldsValue({
      allocated_quota: subAccount.allocated_quota,
      is_active: subAccount.is_active,
    });
    setEditSubModalVisible(true);
  };

  // 提交编辑子账号
  const handleSubmitEditSubAccount = async (values: any) => {
    if (!selectedMaster || !selectedSubAccount) return;

    try {
      await updateSubAccount(selectedMaster.id, selectedSubAccount.id, {
        allocated_quota: values.allocated_quota,
        is_active: values.is_active,
      });
      messageApi.success('子账号更新成功');
      setEditSubModalVisible(false);
      editSubForm.resetFields();
      setSelectedSubAccount(null);
      await loadSubAccounts(selectedMaster.id);
    } catch (error) {
      messageApi.error('更新子账号失败');
    }
  };

  // 删除子账号
  const handleDeleteSubAccount = async (subAccountId: number) => {
    // Find the sub-account record to get the user_id
    const subAccount = subAccounts.find(sub => sub.id === subAccountId);
    if (!subAccount) {
      messageApi.error('找不到子账号信息');
      return;
    }

    modal.confirm({
      title: '删除子账号',
      content: `确定要删除子账号 "${subAccount.username}" 吗？此操作不可撤销。`,
      okText: '确定',
      cancelText: '取消',
      okButtonProps: { danger: true },
      onOk: async () => {
        try {
          // Workaround: Use deleteUser instead of deleteSubAccount due to backend 405 error
          // The /admin/master-accounts/{id}/sub-accounts/{subId} DELETE endpoint is not working
          // So we use /admin/users/{userId} DELETE instead
          await deleteUser(subAccount.user_id);
          messageApi.success('子账号已删除');
          if (selectedMaster) {
            await loadSubAccounts(selectedMaster.id);
            await loadMasterAccounts(); // Refresh master account list to update counts
          }
        } catch (error: any) {
          console.error('删除子账号失败:', error);
          messageApi.error(`删除子账号失败: ${error.response?.data?.detail || error.message}`);
        }
      },
    });
  };

  // 主账号列表列定义
  const masterColumns = [
    {
      title: 'ID',
      dataIndex: 'id',
      key: 'id',
      width: 60,
    },
    {
      title: '用户信息',
      key: 'user_info',
      width: 180,
      render: (record: any) => (
        <div>
          <div style={{ fontWeight: 'bold', marginBottom: 4 }}>
            <UserOutlined style={{ marginRight: 4 }} />
            {record.username || 'Unknown'}
          </div>
          <div style={{ fontSize: '12px', color: '#666' }}>{record.email}</div>
          <div style={{ fontSize: '12px', color: '#999' }}>用户ID: {record.user_id}</div>
        </div>
      ),
    },
    {
      title: '所属组织',
      dataIndex: 'organization',
      key: 'organization',
      width: 150,
      render: (org: string) => (
        <div>
          <TeamOutlined style={{ marginRight: 4 }} />
          {org || '未知组织'}
        </div>
      ),
    },
    {
      title: '核时单价',
      key: 'custom_cpu_hour_price',
      width: 120,
      render: (record: any) => (
        <div>
          <DollarOutlined style={{ marginRight: 4 }} />
          {record.custom_cpu_hour_price ? `¥${record.custom_cpu_hour_price.toFixed(4)}/h` : '未设置'}
        </div>
      ),
    },
    {
      title: '配额使用',
      key: 'quota_usage',
      width: 220,
      render: (record: any) => {
        const balance = record.balance_cpu_hours || 0;
        const total = record.total_cpu_hours || 0;

        // 计算已使用核时: 总核时 - 当前余额
        // 如果余额为负(欠费),则已使用 = 总核时 + |欠费|
        const used = total - balance;

        // 计算使用百分比
        let usagePercent = 0;
        let status: 'success' | 'normal' | 'exception' = 'success';

        if (total > 0) {
          usagePercent = (used / total) * 100;

          // 确定状态
          if (balance < 0) {
            // 欠费状态
            status = 'exception';
          } else if (usagePercent > 90) {
            status = 'exception';
          } else if (usagePercent > 70) {
            status = 'normal';
          } else {
            status = 'success';
          }
        }

        return (
          <div>
            <div style={{ marginBottom: 8 }}>
              <Progress
                percent={Math.min(usagePercent, 100)}
                size="small"
                status={status}
                format={() => `${usagePercent.toFixed(1)}%`}
              />
            </div>
            <div style={{ fontSize: '12px' }}>
              <span style={{
                color: balance >= 0 ? '#52c41a' : '#ff4d4f',
                fontWeight: 'bold'
              }}>
                {balance >= 0 ? '可用' : '欠费'}: {formatQuota(Math.abs(balance))}h
              </span>
              {' / '}
              <span style={{ color: '#ff4d4f' }}>已用: {formatQuota(used)}h</span>
              {' / '}
              <span style={{ color: '#1890ff' }}>总计: {formatQuota(total)}h</span>
            </div>
          </div>
        );
      },
    },
    {
      title: '子账号',
      key: 'sub_accounts',
      width: 120,
      render: (record: any) => {
        const current = record.current_sub_accounts || 0;
        const max = record.max_sub_accounts || 10;
        const percent = (current / max) * 100;

        return (
          <div>
            <Badge
              count={current}
              showZero
              style={{ backgroundColor: current >= max ? '#ff4d4f' : '#52c41a' }}
            />
            <span style={{ marginLeft: 8 }}>/ {max}</span>
            <div style={{ fontSize: '12px', color: '#666', marginTop: 4 }}>
              {percent.toFixed(0)}% 已用
            </div>
          </div>
        );
      },
    },
    {
      title: '状态',
      key: 'status',
      width: 80,
      render: (record: any) => (
        <Tag color={record.is_active ? 'success' : 'error'} icon={record.is_active ? <CheckCircleOutlined /> : <CloseCircleOutlined />}>
          {record.is_active ? '活跃' : '禁用'}
        </Tag>
      ),
    },
    {
      title: '操作',
      key: 'action',
      width: 220,
      fixed: 'right' as const,
      render: (_: any, record: MasterAccount) => (
        <Space size="small">
          <Tooltip title="查看详情">
            <Button
              type="primary"
              size="small"
              icon={<EyeOutlined />}
              onClick={() => handleViewDetail(record)}
            />
          </Tooltip>
          <Tooltip title="编辑配额">
            <Button
              size="small"
              icon={<EditOutlined />}
              onClick={() => handleEditMasterAccount(record)}
            />
          </Tooltip>
          <Tooltip title="删除主账号">
            <Button
              danger
              size="small"
              icon={<DeleteOutlined />}
              onClick={() => handleDeleteMasterAccount(record)}
            />
          </Tooltip>
        </Space>
      ),
    },
  ];

  // 子账号列表列定义
  const subColumns = [
    {
      title: '用户信息',
      key: 'user_info',
      width: 180,
      render: (record: any) => (
        <div>
          <div style={{ fontWeight: 'bold', marginBottom: 4 }}>
            <UserOutlined style={{ marginRight: 4 }} />
            {record.username || `用户 ${record.user_id}`}
          </div>
          <div style={{ fontSize: '12px', color: '#666' }}>{record.email}</div>
          <div style={{ fontSize: '12px', color: '#999' }}>ID: {record.user_id}</div>
        </div>
      ),
    },
    {
      title: '配额信息',
      key: 'quota',
      width: 200,
      render: (record: any) => {
        const personalBalance = record.balance_cpu_hours || 0;
        const allocatedQuota = record.allocated_quota || 0;
        // 子账号总可用 = 个人充值余额 + 主账号分配配额（两个池的总和）
        const totalAvailable = personalBalance + allocatedQuota;

        return (
          <div>
            <div style={{ marginBottom: 8 }}>
              <Typography.Text type="secondary" style={{ fontSize: '11px' }}>个人余额</Typography.Text>
              <div style={{ fontSize: '12px', fontWeight: 'bold' }}>
                {formatQuota(personalBalance)}h
              </div>
            </div>
            <div style={{ marginBottom: 8 }}>
              <Typography.Text type="secondary" style={{ fontSize: '11px' }}>分配配额</Typography.Text>
              <div style={{ fontSize: '12px', fontWeight: 'bold' }}>
                {formatQuota(allocatedQuota)}h
              </div>
            </div>
            <div>
              <Typography.Text type="secondary" style={{ fontSize: '11px' }}>总可用</Typography.Text>
              <div style={{ fontSize: '12px', fontWeight: 'bold', color: totalAvailable > 0 ? '#52c41a' : '#ff4d4f' }}>
                {formatQuota(totalAvailable)}h
              </div>
            </div>
          </div>
        );
      },
    },
    {
      title: '状态',
      key: 'status',
      width: 80,
      render: (record: any) => (
        <Tag color={record.is_active ? 'success' : 'error'}>
          {record.is_active ? '活跃' : '禁用'}
        </Tag>
      ),
    },
    {
      title: '创建时间',
      dataIndex: 'created_at',
      key: 'created_at',
      width: 120,
      render: (date: string) => new Date(date).toLocaleDateString('zh-CN'),
    },
    {
      title: '操作',
      key: 'action',
      width: 150,
      render: (_: any, record: SubAccount) => (
        <Space size="small">
          <Tooltip title="编辑配额">
            <Button
              size="small"
              icon={<EditOutlined />}
              onClick={() => handleEditSubAccount(record)}
            />
          </Tooltip>
          <Tooltip title="删除子账号">
            <Button
              danger
              size="small"
              icon={<DeleteOutlined />}
              onClick={() => handleDeleteSubAccount(record.id)}
            />
          </Tooltip>
        </Space>
      ),
    },
  ];

  // 计算统计数据
  const stats = {
    total: masterAccounts.length,
    active: masterAccounts.filter(a => a.is_active).length,
    totalBalance: masterAccounts.reduce((sum, a) => sum + (a.balance_cpu_hours || 0), 0),
    totalFrozen: masterAccounts.reduce((sum, a) => sum + (a.frozen_cpu_hours || 0), 0),
    totalGranted: masterAccounts.reduce((sum, a) => sum + (a.admin_granted_cpu_hours || 0), 0),
    totalSubAccounts: masterAccounts.reduce((sum, a) => sum + (a.current_sub_accounts || 0), 0),
  };

  return (
    <div style={{ padding: '20px 24px', background: token.colorBgLayout, minHeight: 'calc(100vh - 64px)' }}>
      {contextHolder}
      {contextHolderMessage}
      {/* 页面标题 */}
      <div style={{ marginBottom: 16 }}>
        <Title level={3} style={{ margin: 0, marginBottom: 4 }}>
          <BankOutlined style={{ marginRight: 10, color: token.colorPrimary }} />
          主账号管理
        </Title>
        <Text type="secondary">管理主账号、子账号和配额分配</Text>
      </div>

      <AdminNav />

      <Spin spinning={loading}>
        {/* 统计卡片 - 简洁风格 */}
        <Row gutter={16} style={{ marginBottom: 20 }}>
          {[
            { label: '主账号总数', value: stats.total, suffix: '个', color: '#85a5ff', icon: <TeamOutlined /> },
            { label: '活跃账号', value: `${stats.active}/${stats.total}`, color: '#52c41a', icon: <CheckCircleOutlined />, isText: true },
            { label: '总可用余额', value: `${formatQuota(stats.totalBalance)}h`, color: '#faad14', icon: <DollarOutlined />, isText: true },
            { label: '子账号总数', value: stats.totalSubAccounts, suffix: '个', color: '#b37feb', icon: <UserOutlined /> },
          ].map((item, idx) => (
            <Col xs={12} sm={6} key={idx}>
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
                  <Text strong style={{ fontSize: item.isText ? 18 : 22, color: item.color }}>
                    {item.value}{item.suffix || ''}
                  </Text>
                </div>
              </div>
            </Col>
          ))}
        </Row>

        {/* 操作栏 */}
        <div style={{ marginBottom: 16, display: 'flex', justifyContent: 'space-between' }}>
          <Space>
            <Button
              type="primary"
              icon={<PlusOutlined />}
              onClick={() => {
                setCreateModalVisible(true);
                loadAllUsers();
              }}
            >
              创建主账号
            </Button>
            <Button
              icon={<ReloadOutlined />}
              onClick={loadMasterAccounts}
            >
              刷新
            </Button>
            {selectedRowKeys.length > 0 && (
              <>
                <Tag color="blue">已选择 {selectedRowKeys.length} 个</Tag>
                <Button
                  danger
                  onClick={() => {
                    modal.confirm({
                      title: '批量删除主账号',
                      content: `确定要删除选中的 ${selectedRowKeys.length} 个主账号吗？\n\n此操作将把这些主账号及其子账号转换为个人账号。`,
                      okText: '确定删除',
                      cancelText: '取消',
                      okButtonProps: { danger: true },
                      onOk: async () => {
                        try {
                          for (const id of selectedRowKeys) {
                            await deleteMasterAccount(id as number);
                          }
                          messageApi.success(`已删除 ${selectedRowKeys.length} 个主账号`);
                          setSelectedRowKeys([]);
                          loadMasterAccounts();
                        } catch (error: any) {
                          messageApi.error('批量删除失败: ' + (error.response?.data?.detail || error.message));
                        }
                      },
                    });
                  }}
                >
                  批量删除
                </Button>
                <Button onClick={() => setSelectedRowKeys([])}>
                  取消选择
                </Button>
              </>
            )}
          </Space>
          <Input
            placeholder="搜索用户名、邮箱或组织"
            prefix={<SearchOutlined />}
            style={{ width: 300 }}
            value={searchText}
            onChange={(e) => setSearchText(e.target.value)}
            allowClear
          />
        </div>

        <Table
          columns={masterColumns}
          dataSource={filteredAccounts}
          rowKey="id"
          rowSelection={{
            selectedRowKeys,
            onChange: (keys) => setSelectedRowKeys(keys),
          }}
          pagination={{
            pageSize: 10,
            showTotal: (total) => `共 ${total} 个主账号`,
            showSizeChanger: true,
            showQuickJumper: true,
          }}
          scroll={{ x: 1200 }}
        />
      </Spin >

      {/* 创建主账号对话框 */}
      < Modal
        title="创建主账号"
        open={createModalVisible}
        onOk={() => form.submit()}
        onCancel={() => {
          setCreateModalVisible(false);
          form.resetFields();
          setCreateMode('existing');
        }}
        width={600}
      >
        <Form
          form={form}
          layout="vertical"
          onFinish={handleCreateMasterAccount}
        >
          {/* 创建模式选择 */}
          <Form.Item label="创建方式">
            <Select
              value={createMode}
              onChange={(value) => {
                setCreateMode(value);
                form.resetFields();
                if (value === 'existing') {
                  loadAllUsers();
                }
              }}
              options={[
                { value: 'existing', label: '🔄 升级已有用户为主账号' },
                { value: 'new', label: '➕ 创建新用户并设为主账号' },
              ]}
            />
          </Form.Item>

          {createMode === 'existing' ? (
            <>
              {/* 选择已有用户 */}
              <Form.Item
                label="选择用户"
                name="user_id"
                rules={[{ required: true, message: '请选择用户' }]}
              >
                <Select
                  showSearch
                  placeholder="搜索并选择用户"
                  loading={loadingUsers}
                  filterOption={(input, option) =>
                    (option?.label as string)?.toLowerCase().includes(input.toLowerCase())
                  }
                  onFocus={loadAllUsers}
                  options={allUsers.map((user) => ({
                    value: user.id,
                    label: `${user.username} (${user.email}) - ${user.organization || '无组织'}`,
                  }))}
                />
              </Form.Item>
            </>
          ) : (
            <>
              {/* 创建新用户 */}
              <Form.Item
                label="用户名"
                name="username"
                rules={[
                  { required: true, message: '请输入用户名' },
                  { min: 3, message: '用户名至少3个字符' },
                ]}
              >
                <Input placeholder="用户名" />
              </Form.Item>
              <Form.Item
                label="邮箱"
                name="email"
                rules={[
                  { required: true, message: '请输入邮箱' },
                  { type: 'email', message: '请输入有效的邮箱' },
                ]}
              >
                <Input placeholder="邮箱" />
              </Form.Item>
              <Form.Item
                label="密码"
                name="password"
                rules={[
                  { required: true, message: '请输入密码' },
                  { min: 6, message: '密码至少6个字符' },
                ]}
              >
                <Input.Password placeholder="密码" />
              </Form.Item>
              <Form.Item
                label="初始核时"
                name="initial_cpu_hours"
                tooltip="新用户的初始核时配额"
              >
                <InputNumber style={{ width: '100%' }} placeholder="初始核时数量" min={0} step={100} />
              </Form.Item>
            </>
          )}

          {/* 共用字段 */}
          <Form.Item
            label="组织名称"
            name="organization"
            tooltip="主账号所属的组织名称"
          >
            <Input placeholder="组织名称（可选）" />
          </Form.Item>
          <Form.Item
            label="最大子账号数"
            name="max_sub_accounts"
            initialValue={10}
            tooltip="该主账号最多可创建的子账号数量"
          >
            <InputNumber style={{ width: '100%' }} placeholder="最大子账号数" min={1} max={100} />
          </Form.Item>
        </Form>
      </Modal >

      {/* 编辑主账号对话框 */}
      < Modal
        title={`编辑主账号 - ${selectedMaster?.username}`}
        open={editModalVisible}
        onOk={() => editForm.submit()}
        onCancel={() => {
          setEditModalVisible(false);
          editForm.resetFields();
          setSelectedMaster(null);
        }}
      >
        <Form
          form={editForm}
          layout="vertical"
          onFinish={handleSubmitEdit}
        >
          <Form.Item
            label="可用余额（核时）"
            name="balance_cpu_hours"
            tooltip="调整主账号的可用核时余额"
            rules={[{ required: true, message: '请输入可用余额' }]}
          >
            <InputNumber
              style={{ width: '100%' }}
              placeholder="可用余额"
              min={0}
              step={100}
              precision={QUOTA_PRECISION}
            />
          </Form.Item>
          <Form.Item
            label="最大子账号数"
            name="max_sub_accounts"
            tooltip="该主账号最多可创建的子账号数量"
            rules={[{ required: true, message: '请输入最大子账号数' }]}
          >
            <InputNumber
              style={{ width: '100%' }}
              placeholder="最大子账号数"
              min={1}
              max={100}
            />
          </Form.Item>
          <Form.Item
            label="账号状态"
            name="is_active"
            valuePropName="checked"
            tooltip="禁用后，该主账号及其子账号将无法使用"
          >
            <Switch checkedChildren="启用" unCheckedChildren="禁用" />
          </Form.Item>
        </Form>
      </Modal >

      {/* 创建子账号对话框 */}
      < Modal
        title={`创建子账号 - ${selectedMaster?.username}`}
        open={createSubModalVisible}
        onOk={() => subAccountForm.submit()}
        onCancel={() => {
          setCreateSubModalVisible(false);
          subAccountForm.resetFields();
        }}
      >
        <Form
          form={subAccountForm}
          layout="vertical"
          onFinish={handleCreateSubAccount}
        >
          <Form.Item
            label="用户名"
            name="username"
            rules={[
              { required: true, message: '请输入用户名' },
              { min: 3, message: '用户名至少3个字符' },
            ]}
          >
            <Input placeholder="用户名" />
          </Form.Item>
          <Form.Item
            label="邮箱"
            name="email"
            rules={[
              { required: true, message: '请输入邮箱' },
              { type: 'email', message: '请输入有效的邮箱地址' },
            ]}
          >
            <Input placeholder="邮箱地址" />
          </Form.Item>
          <Form.Item
            label="密码"
            name="password"
            rules={[
              { required: true, message: '请输入密码' },
              { min: 6, message: '密码至少6个字符' },
            ]}
          >
            <Input.Password placeholder="密码" />
          </Form.Item>
          <Form.Item
            label="主账号分配配额（核时）"
            name="allocated_quota"
            tooltip="主账号为该子账号分配的配额"
            rules={[{ type: 'number', min: 0, message: '配额不能为负数' }]}
          >
            <InputNumber
              style={{ width: '100%' }}
              placeholder="请输入主账号分配的配额"
              min={0}
              step={100}
              precision={QUOTA_PRECISION}
              addonAfter="核时"
            />
          </Form.Item>
        </Form>
      </Modal >

      {/* 编辑子账号对话框 */}
      < Modal
        title={`编辑子账号 - ${selectedSubAccount?.username}`}
        open={editSubModalVisible}
        onOk={() => editSubForm.submit()}
        onCancel={() => {
          setEditSubModalVisible(false);
          editSubForm.resetFields();
          setSelectedSubAccount(null);
        }}
        okText="保存"
        cancelText="取消"
      >
        <Form
          form={editSubForm}
          layout="vertical"
          onFinish={handleSubmitEditSubAccount}
        >
          {selectedSubAccount && (
            <div style={{ padding: 12, background: 'rgba(0,0,0,0.02)', borderRadius: 8, marginBottom: 16 }}>
              <div style={{ marginBottom: 8 }}>
                <Typography.Text type="secondary" style={{ fontSize: 12 }}>个人充值余额</Typography.Text>
              </div>
              <div style={{ fontSize: 16, fontWeight: 'bold', marginBottom: 12 }}>
                {formatCpuHours(selectedSubAccount.balance_cpu_hours)} 核时
              </div>
              <div style={{ marginBottom: 8 }}>
                <Typography.Text type="secondary" style={{ fontSize: 12 }}>主账号分配配额</Typography.Text>
              </div>
              <div style={{ fontSize: 16, fontWeight: 'bold' }}>
                {formatQuota(selectedSubAccount.allocated_quota)} 核时
              </div>
            </div>
          )}
          <Form.Item
            label="主账号分配配额（核时）"
            name="allocated_quota"
            tooltip="主账号可以分配给子账号的配额"
            rules={[{ type: 'number', min: 0, message: '配额不能为负数' }]}
          >
            <InputNumber
              style={{ width: '100%' }}
              placeholder="请输入主账号分配的配额"
              min={0}
              step={100}
              precision={QUOTA_PRECISION}
              addonAfter="核时"
            />
          </Form.Item>
          <Form.Item
            label="账号状态"
            name="is_active"
            valuePropName="checked"
          >
            <Switch checkedChildren="启用" unCheckedChildren="禁用" />
          </Form.Item>
        </Form>
      </Modal >

      {/* 主账号详情抽屉 */}
      < Drawer
        title={
          < div >
            <div style={{ fontSize: 18, fontWeight: 'bold' }}>
              {selectedMaster?.username}
            </div>
            <div style={{ fontSize: 12, color: '#666', marginTop: 4 }}>
              {selectedMaster?.email} | 用户ID: {selectedMaster?.user_id}
            </div>
          </div >
        }
        placement="right"
        onClose={() => {
          setDetailDrawerVisible(false);
          setSelectedMaster(null);
        }}
        open={detailDrawerVisible}
        width={720}
      >
        {selectedMaster && (
          <Tabs
            items={[
              {
                key: 'basic',
                label: '基本信息',
                children: (
                  <div>
                    {/* 基本信息 */}
                    <Card title="基本信息" style={{ marginBottom: 16 }}>
                      <Descriptions column={2} size="small">
                        <Descriptions.Item label="所属组织">
                          {selectedMaster.organization || '未知'}
                        </Descriptions.Item>
                        <Descriptions.Item label="账号状态">
                          <Tag color={selectedMaster.is_active ? 'success' : 'error'}>
                            {selectedMaster.is_active ? '活跃' : '禁用'}
                          </Tag>
                        </Descriptions.Item>
                        <Descriptions.Item label="创建时间">
                          {new Date(selectedMaster.created_at).toLocaleString('zh-CN')}
                        </Descriptions.Item>
                        <Descriptions.Item label="更新时间">
                          {new Date(selectedMaster.updated_at).toLocaleString('zh-CN')}
                        </Descriptions.Item>
                      </Descriptions>
                    </Card>

                    {/* 配额统计 */}
                    <Card title="配额统计" style={{ marginBottom: 16 }}>
                      <Row gutter={16}>
                        <Col span={8}>
                          <Statistic
                            title="可用余额"
                            value={selectedMaster.balance_cpu_hours || 0}
                            precision={2}
                            suffix="核时"
                            valueStyle={{ color: '#52c41a' }}
                          />
                        </Col>
                        <Col span={8}>
                          <Statistic
                            title="冻结核时"
                            value={selectedMaster.frozen_cpu_hours || 0}
                            precision={2}
                            suffix="核时"
                            valueStyle={{ color: '#ff4d4f' }}
                          />
                        </Col>
                        <Col span={8}>
                          <Statistic
                            title="管理员赠送"
                            value={selectedMaster.admin_granted_cpu_hours || 0}
                            precision={2}
                            suffix="核时"
                            valueStyle={{ color: '#1890ff' }}
                          />
                        </Col>
                      </Row>
                    </Card>

                    {/* 子账号管理 */}
                    <Card
                      title={
                        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                          <span>
                            子账号列表 ({selectedMaster.current_sub_accounts || 0}/{selectedMaster.max_sub_accounts || 10})
                          </span>
                          <Button
                            type="primary"
                            size="small"
                            icon={<PlusOutlined />}
                            onClick={() => setCreateSubModalVisible(true)}
                            disabled={(selectedMaster.current_sub_accounts || 0) >= (selectedMaster.max_sub_accounts || 10)}
                          >
                            创建子账号
                          </Button>
                        </div>
                      }
                    >
                      {subAccounts.length === 0 ? (
                        <div style={{ textAlign: 'center', padding: '40px 0', color: '#999' }}>
                          <UserOutlined style={{ fontSize: 48, marginBottom: 16 }} />
                          <div>暂无子账号</div>
                        </div>
                      ) : (
                        <Table
                          columns={subColumns}
                          dataSource={subAccounts}
                          rowKey="id"
                          pagination={{ pageSize: 5, size: 'small' }}
                          size="small"
                        />
                      )}
                    </Card>
                  </div>
                ),
              },
              {
                key: 'pricing',
                label: '定价设置',
                children: (
                  <Card title="核时单价设置" style={{ marginBottom: 16 }}>
                    <Form
                      layout="vertical"
                      initialValues={{
                        custom_cpu_hour_price: selectedMaster.custom_cpu_hour_price || undefined,
                      }}
                      onFinish={async (values) => {
                        try {
                          await updateUserPrice({
                            user_id: selectedMaster.user_id,
                            price: values.custom_cpu_hour_price || null,
                          });
                          message.success('定价已更新');
                          loadMasterAccounts();
                        } catch (error) {
                          message.error('更新定价失败');
                        }
                      }}
                    >
                      <Form.Item
                        label="核时单价（元/核时）"
                        name="custom_cpu_hour_price"
                        tooltip="留空表示使用系统默认价格"
                      >
                        <InputNumber
                          placeholder="请输入核时单价"
                          min={0}
                          step={0.01}
                          precision={2}
                          style={{ width: '100%' }}
                        />
                      </Form.Item>
                      <Form.Item>
                        <Button type="primary" htmlType="submit">
                          保存定价
                        </Button>
                      </Form.Item>
                    </Form>
                  </Card>
                ),
              },
            ]}
          />
        )}
      </Drawer >
    </div >
  );
};

export default MasterAccountManagement;

