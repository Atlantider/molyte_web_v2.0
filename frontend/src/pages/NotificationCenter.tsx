/**
 * 消息中心页面
 */
import React, { useState, useEffect } from 'react';
import {
  Card, List, Button, Space, Tag, Empty, Spin, message, Popconfirm, Segmented,
  Typography, Divider, theme, Badge
} from 'antd';
import {
  DeleteOutlined, CheckOutlined, ClearOutlined, BellOutlined,
  WarningOutlined, CheckCircleOutlined, InfoOutlined, ClockCircleOutlined
} from '@ant-design/icons';
import { useThemeStore } from '../stores/themeStore';
import {
  getNotifications,
  markAsRead,
  markAllAsRead,
  deleteNotification,
  deleteAllNotifications,
  Notification
} from '../api/notifications';

const { Title, Text } = Typography;

const NotificationCenter: React.FC = () => {
  const { mode } = useThemeStore();
  const { token } = theme.useToken();
  const [loading, setLoading] = useState(true);
  const [notifications, setNotifications] = useState<Notification[]>([]);
  const [filterType, setFilterType] = useState<'all' | 'unread' | 'read'>('all');
  const [actionLoading, setActionLoading] = useState<number | null>(null);

  const loadNotifications = async () => {
    setLoading(true);
    try {
      const isRead = filterType === 'read' ? true : filterType === 'unread' ? false : undefined;
      const data = await getNotifications({
        skip: 0,
        limit: 100,
        is_read: isRead
      });
      setNotifications(data);
    } catch (error) {
      console.error('加载消息失败:', error);
      message.error('加载消息失败');
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    loadNotifications();
  }, [filterType]);

  const handleMarkAsRead = async (id: number) => {
    setActionLoading(id);
    try {
      await markAsRead(id);
      setNotifications(notifications.map(n => n.id === id ? { ...n, is_read: true } : n));
      message.success('已标记为已读');
    } catch (error) {
      message.error('操作失败');
    } finally {
      setActionLoading(null);
    }
  };

  const handleMarkAllAsRead = async () => {
    try {
      await markAllAsRead();
      setNotifications(notifications.map(n => ({ ...n, is_read: true })));
      message.success('已标记所有消息为已读');
    } catch (error) {
      message.error('操作失败');
    }
  };

  const handleDelete = async (id: number) => {
    try {
      await deleteNotification(id);
      setNotifications(notifications.filter(n => n.id !== id));
      message.success('已删除');
    } catch (error) {
      message.error('删除失败');
    }
  };

  const handleDeleteAll = async () => {
    try {
      await deleteAllNotifications();
      setNotifications([]);
      message.success('已删除所有消息');
    } catch (error) {
      message.error('删除失败');
    }
  };

  const getNotificationIcon = (type: string) => {
    switch (type) {
      case 'DEBT_WARNING':
        return <WarningOutlined style={{ color: '#ff4d4f' }} />;
      case 'DEBT_CLEARED':
        return <CheckCircleOutlined style={{ color: '#52c41a' }} />;
      case 'JOB_COMPLETED':
        return <CheckCircleOutlined style={{ color: '#52c41a' }} />;
      case 'JOB_FAILED':
        return <WarningOutlined style={{ color: '#ff4d4f' }} />;
      case 'QUOTA_LOW':
        return <WarningOutlined style={{ color: '#faad14' }} />;
      case 'QUOTA_RECHARGED':
        return <CheckCircleOutlined style={{ color: '#52c41a' }} />;
      default:
        return <InfoOutlined style={{ color: '#1677ff' }} />;
    }
  };

  const getPriorityColor = (priority: string) => {
    switch (priority) {
      case 'CRITICAL':
        return 'red';
      case 'HIGH':
        return 'orange';
      case 'NORMAL':
        return 'blue';
      case 'LOW':
        return 'default';
      default:
        return 'default';
    }
  };

  const unreadCount = notifications.filter(n => !n.is_read).length;

  return (
    <div style={{
      padding: 24,
      background: token.colorBgLayout,
      minHeight: 'calc(100vh - 64px)',
      transition: 'background 0.3s',
    }}>
      {/* 页面标题 */}
      <div style={{ marginBottom: 24 }}>
        <Title level={2} style={{ margin: 0, marginBottom: 8 }}>
          <BellOutlined style={{ marginRight: 12, color: '#1677ff' }} />
          消息中心
        </Title>
        <Text type="secondary">
          {unreadCount > 0 ? (
            <>
              您有 <Badge count={unreadCount} style={{ backgroundColor: '#ff4d4f' }} /> 条未读消息
            </>
          ) : (
            '所有消息已读'
          )}
        </Text>
      </div>

      {/* 操作栏 */}
      <Card style={{
        borderRadius: 12,
        boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
        border: 'none',
        marginBottom: 24
      }}>
        <Space style={{ width: '100%', justifyContent: 'space-between' }}>
          <div>
            <Text strong>筛选：</Text>
            <Segmented
              value={filterType}
              onChange={(value) => setFilterType(value as any)}
              options={[
                { label: '全部', value: 'all' },
                { label: '未读', value: 'unread' },
                { label: '已读', value: 'read' }
              ]}
              style={{ marginLeft: 12 }}
            />
          </div>
          <Space>
            {unreadCount > 0 && (
              <Button
                type="primary"
                icon={<CheckOutlined />}
                onClick={handleMarkAllAsRead}
              >
                全部标记为已读
              </Button>
            )}
            <Popconfirm
              title="删除所有消息"
              description="确定要删除所有消息吗？此操作不可撤销。"
              onConfirm={handleDeleteAll}
              okText="确定"
              cancelText="取消"
            >
              <Button danger icon={<ClearOutlined />}>
                清空消息
              </Button>
            </Popconfirm>
          </Space>
        </Space>
      </Card>

      {/* 消息列表 */}
      <Card style={{
        borderRadius: 12,
        boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
        border: 'none'
      }}>
        <Spin spinning={loading}>
          {notifications.length === 0 ? (
            <Empty description="暂无消息" style={{ marginTop: 48, marginBottom: 48 }} />
          ) : (
            <List
              dataSource={notifications}
              renderItem={(notification) => (
                <List.Item
                  key={notification.id}
                  style={{
                    padding: '16px',
                    borderRadius: 8,
                    marginBottom: 12,
                    backgroundColor: notification.is_read ? 'transparent' : token.colorBgElevated,
                    border: `1px solid ${token.colorBorder}`,
                    transition: 'all 0.3s'
                  }}
                >
                  <List.Item.Meta
                    avatar={getNotificationIcon(notification.type)}
                    title={
                      <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                        <span>{notification.title}</span>
                        <Tag color={getPriorityColor(notification.priority)}>
                          {notification.priority}
                        </Tag>
                        {!notification.is_read && (
                          <Badge status="processing" text="未读" />
                        )}
                      </div>
                    }
                    description={
                      <div>
                        <Text>{notification.message}</Text>
                        <div style={{ marginTop: 8, fontSize: 12, color: token.colorTextSecondary }}>
                          <ClockCircleOutlined style={{ marginRight: 4 }} />
                          {new Date(notification.created_at).toLocaleString()}
                        </div>
                      </div>
                    }
                  />
                  <Space>
                    {!notification.is_read && (
                      <Button
                        type="text"
                        size="small"
                        icon={<CheckOutlined />}
                        loading={actionLoading === notification.id}
                        onClick={() => handleMarkAsRead(notification.id)}
                      >
                        标记已读
                      </Button>
                    )}
                    <Popconfirm
                      title="删除消息"
                      description="确定要删除这条消息吗？"
                      onConfirm={() => handleDelete(notification.id)}
                      okText="确定"
                      cancelText="取消"
                    >
                      <Button
                        type="text"
                        danger
                        size="small"
                        icon={<DeleteOutlined />}
                      >
                        删除
                      </Button>
                    </Popconfirm>
                  </Space>
                </List.Item>
              )}
            />
          )}
        </Spin>
      </Card>
    </div>
  );
};

export default NotificationCenter;

