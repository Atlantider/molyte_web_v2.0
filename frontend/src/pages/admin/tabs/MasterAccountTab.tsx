/**
 * 主账号管理 Tab
 * 从 MasterAccountManagement 页面提取的核心功能
 */
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
  Popconfirm,
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
} from '@ant-design/icons';
import { formatCpuHours, formatQuota, QUOTA_PRECISION } from '../../../utils/formatQuotaDisplay';
import { updateUserPrice } from '../../../api/billing';
import {
  getAllMasterAccounts,
  getMasterAccountDetail,
  getSubAccounts,
  createMasterAccount,
  createSubAccount,
  updateSubAccount,
  deleteSubAccount,
  deleteMasterAccount,
  MasterAccount,
  SubAccount,
  CreateMasterAccountRequest,
  CreateSubAccountRequest,
} from '../../../api/admin';

const { Text } = Typography;

const MasterAccountTab: React.FC = () => {
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

  // 加载主账号列表
  const loadMasterAccounts = async () => {
    setLoading(true);
    try {
      const response = await getAllMasterAccounts();
      setMasterAccounts(response || []);
      setFilteredAccounts(response || []);
    } catch (error: any) {
      message.error(`加载主账号列表失败: ${error.response?.data?.detail || error.message}`);
    } finally {
      setLoading(false);
    }
  };

  // 加载子账号列表
  const loadSubAccounts = async (masterId: number) => {
    try {
      const response = await getSubAccounts(masterId);
      setSubAccounts(response || []);
    } catch (error: any) {
      message.error(`加载子账号列表失败: ${error.response?.data?.detail || error.message}`);
    }
  };

  useEffect(() => {
    loadMasterAccounts();
  }, []);

  // 搜索过滤
  useEffect(() => {
    const filtered = masterAccounts.filter(account =>
      account.username?.toLowerCase().includes(searchText.toLowerCase()) ||
      account.email?.toLowerCase().includes(searchText.toLowerCase()) ||
      account.organization?.toLowerCase().includes(searchText.toLowerCase())
    );
    setFilteredAccounts(filtered);
  }, [searchText, masterAccounts]);

  // 创建主账号
  const handleCreateMasterAccount = async (values: CreateMasterAccountRequest) => {
    try {
      await createMasterAccount(values);
      message.success('主账号创建成功');
      setCreateModalVisible(false);
      form.resetFields();
      loadMasterAccounts();
    } catch (error) {
      message.error('创建主账号失败');
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

  // 删除主账号
  const handleDeleteMasterAccount = async (id: number) => {
    try {
      await deleteMasterAccount(id);
      message.success('主账号删除成功');
      loadMasterAccounts();
    } catch (error: any) {
      message.error(`删除主账号失败: ${error.response?.data?.detail || error.message}`);
    }
  };

  const masterColumns = [
    {
      title: '用户名',
      dataIndex: 'username',
      key: 'username',
      width: 150,
    },
    {
      title: '邮箱',
      dataIndex: 'email',
      key: 'email',
      width: 200,
    },
    {
      title: '组织',
      dataIndex: 'organization',
      key: 'organization',
      width: 150,
    },
    {
      title: '可用余额',
      dataIndex: 'balance_cpu_hours',
      key: 'balance_cpu_hours',
      width: 150,
      render: (balance: number) => (
        <span style={{ fontWeight: 'bold', color: '#1890ff' }}>
          {formatCpuHours(balance)}
        </span>
      ),
    },
    {
      title: '子账号数',
      dataIndex: 'current_sub_accounts',
      key: 'current_sub_accounts',
      width: 100,
      render: (count: number, record: MasterAccount) => (
        <span>{count} / {record.max_sub_accounts}</span>
      ),
    },
    {
      title: '状态',
      dataIndex: 'is_active',
      key: 'is_active',
      width: 100,
      render: (active: boolean) => (
        <Tag color={active ? 'green' : 'red'}>
          {active ? '活跃' : '禁用'}
        </Tag>
      ),
    },
    {
      title: '操作',
      key: 'action',
      width: 200,
      render: (_: any, record: MasterAccount) => (
        <Space size="small">
          <Button
            type="primary"
            size="small"
            icon={<EyeOutlined />}
            onClick={() => handleViewDetail(record)}
          >
            查看
          </Button>
          <Button
            size="small"
            icon={<EditOutlined />}
            onClick={() => handleEditMasterAccount(record)}
          >
            编辑
          </Button>
          <Popconfirm
            title="确定删除？"
            description="删除后无法恢复，请谨慎操作"
            onConfirm={() => handleDeleteMasterAccount(record.id)}
            okText="确定"
            cancelText="取消"
          >
            <Button
              danger
              size="small"
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
    <div>
      {/* 统计卡片 */}
      <Row gutter={[16, 16]} style={{ marginBottom: '24px' }}>
        <Col xs={24} sm={12} md={8}>
          <Card
            style={{
              borderRadius: '12px',
              boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
              border: 'none',
              background: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)',
              color: 'white',
            }}
          >
            <Statistic
              title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>主账号总数</span>}
              value={masterAccounts.length}
              valueStyle={{ color: 'white', fontSize: '28px', fontWeight: 'bold' }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={8}>
          <Card
            style={{
              borderRadius: '12px',
              boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
              border: 'none',
              background: 'linear-gradient(135deg, #f093fb 0%, #f5576c 100%)',
              color: 'white',
            }}
          >
            <Statistic
              title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>活跃账号</span>}
              value={masterAccounts.filter(a => a.is_active).length}
              suffix={`/ ${masterAccounts.length}`}
              valueStyle={{ color: 'white', fontSize: '28px', fontWeight: 'bold' }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={8}>
          <Card
            style={{
              borderRadius: '12px',
              boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
              border: 'none',
              background: 'linear-gradient(135deg, #4facfe 0%, #00f2fe 100%)',
              color: 'white',
            }}
          >
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
              <div>
                <div style={{ color: 'rgba(255,255,255,0.85)', marginBottom: '8px' }}>总子账号数</div>
                <div style={{ fontSize: '28px', fontWeight: 'bold', color: 'white' }}>
                  {masterAccounts.reduce((sum, a) => sum + (a.current_sub_accounts || 0), 0)}
                </div>
              </div>
              <Button
                type="primary"
                icon={<ReloadOutlined />}
                onClick={loadMasterAccounts}
                loading={loading}
                style={{ background: 'rgba(255,255,255,0.3)', border: 'none' }}
              >
                刷新
              </Button>
            </div>
          </Card>
        </Col>
      </Row>

      {/* 操作栏 */}
      <div style={{ marginBottom: 16, display: 'flex', justifyContent: 'space-between' }}>
        <Space>
          <Button
            type="primary"
            icon={<PlusOutlined />}
            onClick={() => setCreateModalVisible(true)}
          >
            创建主账号
          </Button>
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

      {/* 主账号表格 */}
      <Spin spinning={loading}>
        <Table
          columns={masterColumns}
          dataSource={filteredAccounts}
          rowKey="id"
          pagination={{
            pageSize: 10,
            showTotal: (total) => `共 ${total} 个主账号`,
            showSizeChanger: true,
            showQuickJumper: true,
          }}
          scroll={{ x: 1200 }}
        />
      </Spin>
    </div>
  );
};

export default MasterAccountTab;

