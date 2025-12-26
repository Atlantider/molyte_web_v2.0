/**
 * 用户定价管理 Tab
 */
import React, { useState, useEffect } from 'react';
import {
  Table,
  Button,
  Modal,
  Form,
  Input,
  InputNumber,
  Space,
  message,
  Card,
  Row,
  Col,
  Statistic,
  Tag,
  Tooltip,
  Popconfirm,
  Divider,
  Typography,
  Alert,
} from 'antd';
import {
  EditOutlined,
  DeleteOutlined,
  ReloadOutlined,
  SettingOutlined,
} from '@ant-design/icons';
import apiClient from '../../../api/client';
import type { User } from '../../../types';

const { Text } = Typography;

interface UserPricingInfo {
  user_id: number;
  username: string;
  email: string;
  custom_cpu_hour_price: number | null;
  global_price: number;
  effective_price: number;
  price_updated_at: string | null;
  price_updated_by: number | null;
}

const UserPricingTab: React.FC = () => {
  const [users, setUsers] = useState<UserPricingInfo[]>([]);
  const [loading, setLoading] = useState(false);
  const [globalPrice, setGlobalPrice] = useState<number>(0.1);
  const [isModalVisible, setIsModalVisible] = useState(false);
  const [selectedUser, setSelectedUser] = useState<UserPricingInfo | null>(null);
  const [form] = Form.useForm();

  // 批量操作状态
  const [selectedRowKeys, setSelectedRowKeys] = useState<React.Key[]>([]);
  const [batchModalVisible, setBatchModalVisible] = useState(false);
  const [batchForm] = Form.useForm();
  const [batchLoading, setBatchLoading] = useState(false);

  // 加载用户列表
  const loadUsers = async () => {
    setLoading(true);
    try {
      const response = await apiClient.get('/admin/users?skip=0&limit=1000');
      if (response.data) {
        const usersWithPricing = await Promise.all(
          response.data.map(async (user: User) => {
            try {
              const pricingRes = await apiClient.get(`/billing/admin/user-pricing/${user.id}`);
              return pricingRes.data;
            } catch (error) {
              console.error(`Failed to load pricing for user ${user.id}:`, error);
              return null;
            }
          })
        );
        setUsers(usersWithPricing.filter(Boolean));
      }
    } catch (error) {
      message.error('加载用户列表失败');
      console.error(error);
    } finally {
      setLoading(false);
    }
  };

  // 加载全局定价
  const loadGlobalPrice = async () => {
    try {
      const response = await apiClient.get('/billing/admin/pricing');
      if (response.data) {
        setGlobalPrice(response.data.cpu_hour_price);
      }
    } catch (error) {
      console.error('Failed to load global pricing:', error);
    }
  };

  useEffect(() => {
    loadUsers();
    loadGlobalPrice();
  }, []);

  // 打开编辑对话框
  const handleEdit = (user: UserPricingInfo) => {
    setSelectedUser(user);
    form.setFieldsValue({
      custom_cpu_hour_price: user.custom_cpu_hour_price,
    });
    setIsModalVisible(true);
  };

  // 删除自定义价格
  const handleDelete = async (userId: number) => {
    try {
      const response = await apiClient.put(`/billing/admin/user-pricing/${userId}`, {
        custom_cpu_hour_price: null,
      });
      if (response.data?.success) {
        message.success(response.data.message);
        loadUsers();
      }
    } catch (error) {
      message.error('删除自定义价格失败');
      console.error(error);
    }
  };

  // 保存定价
  const handleSave = async (values: any) => {
    if (!selectedUser) return;

    try {
      const response = await apiClient.put(
        `/billing/admin/user-pricing/${selectedUser.user_id}`,
        {
          custom_cpu_hour_price: values.custom_cpu_hour_price,
        }
      );

      if (response.data?.success) {
        message.success(response.data.message);
        setIsModalVisible(false);
        form.resetFields();
        loadUsers();
      }
    } catch (error) {
      message.error('保存定价失败');
      console.error(error);
    }
  };

  // 批量设置定价
  const handleBatchSave = async () => {
    if (selectedRowKeys.length === 0) {
      message.warning('请先选择用户');
      return;
    }

    try {
      const values = batchForm.getFieldsValue();
      setBatchLoading(true);

      const response = await apiClient.put('/billing/admin/user-pricing/batch', {
        user_ids: selectedRowKeys,
        custom_cpu_hour_price: values.batch_price,
      });

      if (response.data) {
        message.success(response.data.message);
        setBatchModalVisible(false);
        batchForm.resetFields();
        setSelectedRowKeys([]);
        loadUsers();
      }
    } catch (error: any) {
      message.error('批量设置定价失败: ' + (error.response?.data?.detail || '请稍后重试'));
      console.error(error);
    } finally {
      setBatchLoading(false);
    }
  };

  // 批量恢复全局定价
  const handleBatchResetToGlobal = async () => {
    if (selectedRowKeys.length === 0) {
      message.warning('请先选择用户');
      return;
    }

    try {
      setBatchLoading(true);
      const response = await apiClient.put('/billing/admin/user-pricing/batch', {
        user_ids: selectedRowKeys,
        custom_cpu_hour_price: null,
      });

      if (response.data) {
        message.success(`已将 ${response.data.success_count} 个用户恢复为全局定价`);
        setSelectedRowKeys([]);
        loadUsers();
      }
    } catch (error: any) {
      message.error('批量恢复失败: ' + (error.response?.data?.detail || '请稍后重试'));
      console.error(error);
    } finally {
      setBatchLoading(false);
    }
  };

  // 表格行选择配置
  const rowSelection = {
    selectedRowKeys,
    onChange: (keys: React.Key[]) => setSelectedRowKeys(keys),
  };

  const columns = [
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
      title: '全局定价',
      dataIndex: 'global_price',
      key: 'global_price',
      width: 120,
      render: (price: number | undefined) => price !== undefined && price !== null ? `¥${price.toFixed(4)}/核时` : '-',
    },
    {
      title: '自定义定价',
      dataIndex: 'custom_cpu_hour_price',
      key: 'custom_cpu_hour_price',
      width: 150,
      render: (price: number | null | undefined) =>
        price !== null && price !== undefined ? (
          <Tag color="blue">¥{price.toFixed(4)}/核时</Tag>
        ) : (
          <Tag>使用全局定价</Tag>
        ),
    },
    {
      title: '生效定价',
      dataIndex: 'effective_price',
      key: 'effective_price',
      width: 150,
      render: (price: number | undefined) => (
        <span style={{ fontWeight: 'bold', color: '#1890ff' }}>
          {price !== undefined && price !== null ? `¥${price.toFixed(4)}/核时` : '-'}
        </span>
      ),
    },
    {
      title: '最后更新',
      dataIndex: 'price_updated_at',
      key: 'price_updated_at',
      width: 180,
      render: (date: string | null) =>
        date ? new Date(date).toLocaleString('zh-CN') : '-',
    },
    {
      title: '操作',
      key: 'action',
      width: 150,
      render: (_: any, record: UserPricingInfo) => (
        <Space size="small">
          <Button
            type="primary"
            size="small"
            icon={<EditOutlined />}
            onClick={() => handleEdit(record)}
          >
            编辑
          </Button>
          {record.custom_cpu_hour_price !== null && (
            <Popconfirm
              title="确认删除"
              description="确定要删除该用户的自定义定价吗？将恢复使用全局定价。"
              onConfirm={() => handleDelete(record.user_id)}
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
          )}
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
              title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>全局核时单价</span>}
              value={globalPrice}
              prefix="¥"
              suffix="/核时"
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
              title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>自定义定价用户</span>}
              value={users.filter(u => u.custom_cpu_hour_price !== null).length}
              suffix={`/ ${users.length}`}
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
                <div style={{ color: 'rgba(255,255,255,0.85)', marginBottom: '8px' }}>用户总数</div>
                <div style={{ fontSize: '28px', fontWeight: 'bold', color: 'white' }}>
                  {users.length}
                </div>
              </div>
              <Button
                type="primary"
                icon={<ReloadOutlined />}
                onClick={loadUsers}
                loading={loading}
                style={{ background: 'rgba(255,255,255,0.3)', border: 'none' }}
              >
                刷新
              </Button>
            </div>
          </Card>
        </Col>
      </Row>

      {/* 批量操作按钮 */}
      {selectedRowKeys.length > 0 && (
        <Alert
          message={
            <Space>
              <span>已选择 {selectedRowKeys.length} 个用户</span>
              <Button
                type="primary"
                icon={<SettingOutlined />}
                onClick={() => setBatchModalVisible(true)}
                size="small"
              >
                批量设置单价
              </Button>
              <Popconfirm
                title="确认恢复"
                description={`确定要将 ${selectedRowKeys.length} 个用户恢复为全局定价吗？`}
                onConfirm={handleBatchResetToGlobal}
                okText="确定"
                cancelText="取消"
              >
                <Button size="small" loading={batchLoading}>
                  恢复全局定价
                </Button>
              </Popconfirm>
              <Button size="small" onClick={() => setSelectedRowKeys([])}>
                取消选择
              </Button>
            </Space>
          }
          type="info"
          style={{ marginBottom: 16 }}
        />
      )}

      {/* 用户列表 */}
      <Table
        columns={columns}
        dataSource={users}
        loading={loading}
        rowKey="user_id"
        rowSelection={rowSelection}
        pagination={{
          pageSize: 20,
          showSizeChanger: true,
          showQuickJumper: true,
          showTotal: (total) => `共 ${total} 个用户`,
        }}
        scroll={{ x: 1200 }}
      />

      {/* 编辑对话框 */}
      <Modal
        title={`编辑用户定价 - ${selectedUser?.username}`}
        open={isModalVisible}
        onOk={() => form.submit()}
        onCancel={() => {
          setIsModalVisible(false);
          form.resetFields();
        }}
        width={500}
      >
        <Form
          form={form}
          layout="vertical"
          onFinish={handleSave}
        >
          <Form.Item label="用户名">
            <Input value={selectedUser?.username} disabled />
          </Form.Item>
          <Form.Item label="邮箱">
            <Input value={selectedUser?.email} disabled />
          </Form.Item>
          <Form.Item label="全局定价">
            <Input value={`¥${(globalPrice || 0.1).toFixed(4)}/核时`} disabled />
          </Form.Item>
          <Divider />
          <Form.Item 
            label="自定义核时单价 (元/核时)" 
            tooltip="设置此用户的自定义价格，将覆盖全局定价。新价格将在下一次充值时生效。"
          >
            <InputNumber
              value={form.getFieldValue('custom_cpu_hour_price')}
              onChange={(val) => form.setFieldValue('custom_cpu_hour_price', val)}
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

      {/* 批量设置定价对话框 */}
      <Modal
        title={`批量设置定价 - 已选择 ${selectedRowKeys.length} 个用户`}
        open={batchModalVisible}
        onOk={handleBatchSave}
        onCancel={() => {
          setBatchModalVisible(false);
          batchForm.resetFields();
        }}
        confirmLoading={batchLoading}
        width={500}
      >
        <Form
          form={batchForm}
          layout="vertical"
        >
          <Form.Item label="当前全局定价">
            <Input value={`¥${(globalPrice || 0.1).toFixed(4)}/核时`} disabled />
          </Form.Item>
          <Divider />
          <Form.Item
            name="batch_price"
            label="批量设置核时单价 (元/核时)"
            tooltip="设置选中用户的自定义价格，将覆盖全局定价。留空则恢复为全局定价。"
          >
            <InputNumber
              min={0.001}
              step={0.01}
              precision={4}
              style={{ width: '100%' }}
              placeholder="输入自定义价格，留空则恢复全局定价"
            />
          </Form.Item>
          <Alert
            message="批量设置说明"
            description={
              <ul style={{ margin: 0, paddingLeft: 16 }}>
                <li>输入价格将为选中的 {selectedRowKeys.length} 个用户设置统一的自定义单价</li>
                <li>留空则将选中用户恢复为使用全局定价</li>
                <li>新价格将在用户下一次充值时生效</li>
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

export default UserPricingTab;

