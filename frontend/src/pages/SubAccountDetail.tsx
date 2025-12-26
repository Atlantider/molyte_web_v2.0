/**
 * 子账号详情页面
 */

import React, { useState, useEffect } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import {
  Card,
  Row,
  Col,
  Statistic,
  Button,
  Space,
  message,
  Spin,
  Divider,
  Tag,
  Table,
  Modal,
  Form,
  InputNumber,
  Switch,
  Typography,
  Alert,
} from 'antd';
import {
  ArrowLeftOutlined,
  EditOutlined,
  DeleteOutlined,
  UserOutlined,
  WalletOutlined,
  CrownOutlined,
  HistoryOutlined,
} from '@ant-design/icons';
import { getMySubAccounts, updateMySubAccount, deleteMySubAccount, getSubAccountJobs, SubAccount, SubAccountJob, SubAccountJobsResponse } from '../api/accounts';
import { useThemeStore } from '../stores/themeStore';
import { formatCpuHours } from '../utils/formatQuotaDisplay';

const { Title, Text } = Typography;

const SubAccountDetail: React.FC = () => {
  const { id } = useParams<{ id: string }>();
  const navigate = useNavigate();
  const { mode } = useThemeStore();
  const isDark = mode === 'dark';



  const [subAccount, setSubAccount] = useState<SubAccount | null>(null);
  const [loading, setLoading] = useState(true);
  const [editModalVisible, setEditModalVisible] = useState(false);
  const [jobsModalVisible, setJobsModalVisible] = useState(false);
  const [jobs, setJobs] = useState<SubAccountJob[]>([]);
  const [jobsLoading, setJobsLoading] = useState(false);
  const [form] = Form.useForm();

  // 加载子账号详情
  const loadSubAccountDetail = async () => {
    if (!id) return;
    try {
      setLoading(true);
      const data = await getMySubAccounts();
      const account = data.find(a => a.id === parseInt(id));
      if (!account) {
        message.error('子账号不存在');
        navigate('/workspace/account-center?tab=organization');
        return;
      }
      setSubAccount(account);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '加载子账号详情失败');
      navigate('/workspace/account-center?tab=organization');
    } finally {
      setLoading(false);
    }
  };

  // 加载子账号任务
  const loadSubAccountJobs = async () => {
    if (!id) return;
    try {
      setJobsLoading(true);
      const data: SubAccountJobsResponse = await getSubAccountJobs(parseInt(id));
      setJobs(data.jobs);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '加载任务列表失败');
    } finally {
      setJobsLoading(false);
    }
  };

  useEffect(() => {
    loadSubAccountDetail();
  }, [id]);

  // 编辑子账号
  const handleEdit = () => {
    if (!subAccount) return;
    form.setFieldsValue({
      allocated_quota: subAccount.allocated_quota,
      is_active: subAccount.is_active,
    });
    setEditModalVisible(true);
  };

  // 提交编辑
  const handleUpdateSubAccount = async (values: any) => {
    if (!id) return;
    try {
      await updateMySubAccount(parseInt(id), values);
      message.success('子账号更新成功');
      setEditModalVisible(false);
      loadSubAccountDetail();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '更新失败');
    }
  };

  // 删除子账号
  const handleDelete = () => {
    if (!id) return;
    Modal.confirm({
      title: '删除子账号',
      content: '确定要删除这个子账号吗？此操作不可撤销。',
      okText: '确定',
      cancelText: '取消',
      okButtonProps: { danger: true },
      onOk: async () => {
        try {
          await deleteMySubAccount(parseInt(id));
          message.success('子账号删除成功');
          navigate('/workspace/account-center?tab=organization');
        } catch (error: any) {
          message.error(error.response?.data?.detail || '删除失败');
        }
      },
    });
  };



  if (loading) {
    return (
      <div style={{ padding: 24, textAlign: 'center' }}>
        <Spin size="large" />
      </div>
    );
  }

  if (!subAccount) {
    return (
      <div style={{ padding: 24 }}>
        <Alert message="子账号不存在" type="error" showIcon />
      </div>
    );
  }

  const jobColumns = [
    { title: '任务ID', dataIndex: 'id', key: 'id' },
    { title: '任务名称', dataIndex: 'job_name', key: 'job_name' },
    { title: '状态', dataIndex: 'status', key: 'status' },
    { title: '消耗核时', dataIndex: 'cpu_hours_used', key: 'cpu_hours_used', render: (v: number) => v?.toFixed(2) },
  ];

  return (
    <div style={{ padding: 24, background: isDark ? '#141414' : '#f5f7fa', minHeight: '100vh' }}>
      <Button icon={<ArrowLeftOutlined />} onClick={() => navigate('/workspace/account-center?tab=organization')} style={{ marginBottom: 24 }}>
        返回
      </Button>

      <Card title={`子账号详情 - ${subAccount.username}`}>
        <Row gutter={[16, 16]} style={{ marginBottom: 24 }}>
          <Col xs={24} sm={12} md={6}>
            <Card bordered={false}>
              <Statistic title="用户名" value={subAccount.username} prefix={<UserOutlined />} />
            </Card>
          </Col>
          <Col xs={24} sm={12} md={6}>
            <Card bordered={false}>
              <Statistic title="邮箱" value={subAccount.email} />
            </Card>
          </Col>
          <Col xs={24} sm={12} md={6}>
            <Card bordered={false}>
              <Statistic title="分配配额" value={formatCpuHours(subAccount.allocated_quota)} suffix="核时" />
            </Card>
          </Col>
          <Col xs={24} sm={12} md={6}>
            <Card bordered={false}>
              <Statistic title="个人余额" value={formatCpuHours(subAccount.balance_cpu_hours)} suffix="核时" />
            </Card>
          </Col>
        </Row>

        <Divider />

        <Space style={{ marginBottom: 16 }}>
          <Button type="primary" icon={<EditOutlined />} onClick={handleEdit}>
            编辑
          </Button>
          <Button icon={<HistoryOutlined />} onClick={() => { setJobsModalVisible(true); loadSubAccountJobs(); }}>
            查看任务
          </Button>
          <Button danger icon={<DeleteOutlined />} onClick={handleDelete}>
            删除
          </Button>
        </Space>

        <Tag color={subAccount.is_active ? 'green' : 'red'}>
          {subAccount.is_active ? '活跃' : '禁用'}
        </Tag>
      </Card>

      {/* 编辑对话框 */}
      <Modal
        title="编辑子账号"
        open={editModalVisible}
        onOk={() => form.submit()}
        onCancel={() => setEditModalVisible(false)}
      >
        <Form form={form} layout="vertical" onFinish={handleUpdateSubAccount}>
          <Form.Item label="分配配额" name="allocated_quota">
            <InputNumber min={0} step={0.1} precision={2} />
          </Form.Item>
          <Form.Item label="账号状态" name="is_active" valuePropName="checked">
            <Switch checkedChildren="启用" unCheckedChildren="禁用" />
          </Form.Item>
        </Form>
      </Modal>

      {/* 任务列表对话框 */}
      <Modal
        title={`${subAccount.username} 的任务列表`}
        open={jobsModalVisible}
        onCancel={() => setJobsModalVisible(false)}
        width={900}
        footer={null}
      >
        <Table columns={jobColumns} dataSource={jobs} rowKey="id" loading={jobsLoading} />
      </Modal>
    </div>
  );
};

export default SubAccountDetail;

