/**
 * 用户定价管理页面
 * 允许管理员为不同用户设置不同的核时单价
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
  Select,
  Typography,
  theme,
} from 'antd';
import {
  EditOutlined,
  DeleteOutlined,
  DollarOutlined,
  ReloadOutlined,
  TeamOutlined,
  UserOutlined,
} from '@ant-design/icons';
import apiClient from '../../api/client';
import AdminNav from '../../components/AdminNav';
import type { User } from '../../types';

const { Title, Text } = Typography;

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

const UserPricing: React.FC = () => {
  const { token } = theme.useToken();
  const [users, setUsers] = useState<UserPricingInfo[]>([]);
  const [loading, setLoading] = useState(false);
  const [globalPrice, setGlobalPrice] = useState<number>(0.1);
  const [isModalVisible, setIsModalVisible] = useState(false);
  const [selectedUser, setSelectedUser] = useState<UserPricingInfo | null>(null);
  const [form] = Form.useForm();

  // 加载用户列表
  const loadUsers = async () => {
    setLoading(true);
    try {
      const response = await apiClient.get('/admin/users?skip=0&limit=1000');
      if (response.data) {
        // 获取每个用户的定价信息
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

  // 删除自定义价格（恢复全局定价）
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
      render: (price: number | undefined) => price ? `¥${price.toFixed(4)}/核时` : '-',
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
          {price ? `¥${price.toFixed(4)}/核时` : '-'}
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
    <div style={{ padding: '20px 24px', background: token.colorBgLayout, minHeight: 'calc(100vh - 64px)' }}>
      {/* 页面标题 */}
      <div style={{ marginBottom: 16 }}>
        <Title level={3} style={{ margin: 0, marginBottom: 4 }}>
          <DollarOutlined style={{ marginRight: 10, color: token.colorPrimary }} />
          用户定价
        </Title>
        <Text type="secondary">管理用户个性化核时单价</Text>
      </div>

      <AdminNav />

      {/* 统计卡片 - 简洁风格 */}
      <Row gutter={16} style={{ marginBottom: 20 }}>
        {[
          { label: '全局单价', value: `¥${globalPrice}/核时`, color: '#85a5ff', icon: <DollarOutlined />, isText: true },
          { label: '自定义定价', value: `${users.filter(u => u.custom_cpu_hour_price !== null).length}/${users.length}`, color: '#f5576c', icon: <UserOutlined />, isText: true },
          { label: '用户总数', value: users.length, color: '#52c41a', icon: <TeamOutlined /> },
        ].map((item, idx) => (
          <Col xs={24} sm={8} key={idx}>
            <div style={{
              padding: '16px 20px',
              background: token.colorBgContainer,
              borderRadius: 10,
              border: `1px solid ${token.colorBorder}`,
              borderLeft: `4px solid ${item.color}`,
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'space-between',
            }}>
              <div style={{ display: 'flex', alignItems: 'center', gap: 14 }}>
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
              {idx === 2 && (
                <Button icon={<ReloadOutlined />} onClick={loadUsers} loading={loading}>刷新</Button>
              )}
            </div>
          </Col>
        ))}
      </Row>

      {/* 用户列表 */}
      <Card
        style={{
          borderRadius: '12px',
          boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
          border: 'none',
        }}
      >
        <Table
          columns={columns}
          dataSource={users}
          loading={loading}
          rowKey="user_id"
          pagination={{ pageSize: 20 }}
          scroll={{ x: 1200 }}
        />
      </Card>

      {/* 编辑对话框 */}
      <Modal
        title={`编辑用户定价 - ${selectedUser?.username}`}
        visible={isModalVisible}
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
          <Form.Item
            label="全局定价"
            name="global_price"
            initialValue={globalPrice}
          >
            <InputNumber
              disabled
              prefix="¥"
              suffix="/核时"
              step={0.01}
              min={0}
            />
          </Form.Item>

          <Form.Item
            label="自定义定价（留空则使用全局定价）"
            name="custom_cpu_hour_price"
            rules={[
              {
                validator: (_, value) => {
                  if (value !== undefined && value !== null && value <= 0) {
                    return Promise.reject(new Error('定价必须大于 0'));
                  }
                  return Promise.resolve();
                },
              },
            ]}
          >
            <InputNumber
              placeholder="输入自定义定价，或留空使用全局定价"
              prefix="¥"
              suffix="/核时"
              step={0.01}
              min={0.01}
              precision={4}
            />
          </Form.Item>

          <Form.Item
            label="说明"
            name="reason"
          >
            <Input.TextArea
              placeholder="可选：说明为什么要调整该用户的定价"
              rows={3}
            />
          </Form.Item>
        </Form>
      </Modal>
    </div>
  );
};

export default UserPricing;

