/**
 * 消息/通知 API
 */
import client from './client';

export interface Notification {
  id: number;
  type: string;
  title: string;
  message: string;
  priority: string;
  is_read: boolean;
  read_at?: string;
  related_id?: string;
  related_type?: string;
  created_at: string;
  updated_at: string;
}

/**
 * 获取消息列表
 */
export const getNotifications = async (params?: {
  skip?: number;
  limit?: number;
  is_read?: boolean;
  type?: string;
}): Promise<Notification[]> => {
  const response = await client.get('/notifications/', { params });
  return response.data;
};

/**
 * 获取未读消息数
 */
export const getUnreadCount = async (): Promise<{ unread_count: number }> => {
  const response = await client.get('/notifications/unread-count');
  return response.data;
};

/**
 * 获取单个消息
 */
export const getNotification = async (notificationId: number): Promise<Notification> => {
  const response = await client.get(`/notifications/${notificationId}`);
  return response.data;
};

/**
 * 标记消息为已读
 */
export const markAsRead = async (notificationId: number): Promise<Notification> => {
  const response = await client.put(`/notifications/${notificationId}/read`);
  return response.data;
};

/**
 * 标记所有消息为已读
 */
export const markAllAsRead = async (): Promise<{ status: string }> => {
  const response = await client.put('/notifications/mark-all-as-read');
  return response.data;
};

/**
 * 删除消息
 */
export const deleteNotification = async (notificationId: number): Promise<{ status: string }> => {
  const response = await client.delete(`/notifications/${notificationId}`);
  return response.data;
};

/**
 * 删除所有消息
 */
export const deleteAllNotifications = async (): Promise<{ status: string }> => {
  const response = await client.delete('/notifications/');
  return response.data;
};

