/**
 * useQuota Hook
 * 
 * 用于获取和管理用户配额信息
 * 支持多种配额来源：个人配额、组织配额、主账号配额
 */

import { useState, useEffect, useCallback } from 'react';
import { message } from 'antd';
import client from '../api/client';

export interface QuotaSource {
  personal?: number;
  organization?: number;
  master_account?: number;
  personal_quota?: number;
  organization_quota?: number;
  master_account_quota?: number;
  free_granted?: number;  // 初始赠送
  recharge?: number;  // 充值获得
  admin_granted?: number;  // 管理员赠送
}

export interface QuotaData {
  available_quota: number;
  quota_sources: QuotaSource;
  account_type: string;
  balance_cpu_hours?: number;
  frozen_cpu_hours?: number;
  debt_cpu_hours?: number;
  total_cpu_hours?: number;
  used_cpu_hours?: number;
  account_details?: {
    [key: string]: any;
  };
}

interface UseQuotaReturn {
  quota: QuotaData | null;
  loading: boolean;
  error: string | null;
  refetch: () => Promise<void>;
  getAvailableQuota: () => number;
  getQuotaPercentage: () => number;
  isQuotaSufficient: (requiredHours: number) => boolean;
  getQuotaSourceBreakdown: () => QuotaSource;
}

/**
 * 获取用户配额信息
 */
export const useQuota = (): UseQuotaReturn => {
  const [quota, setQuota] = useState<QuotaData | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  const fetchQuota = useCallback(async () => {
    try {
      setLoading(true);
      setError(null);

      const response = await client.get('/users/me/quota');
      setQuota(response.data);
    } catch (err: any) {
      const errorMsg = err.response?.data?.detail || '获取配额信息失败';
      setError(errorMsg);
      message.error(errorMsg);
    } finally {
      setLoading(false);
    }
  }, []);

  useEffect(() => {
    fetchQuota();
  }, [fetchQuota]);

  const getAvailableQuota = useCallback((): number => {
    return quota?.available_quota ?? 0;
  }, [quota]);

  const getQuotaPercentage = useCallback((): number => {
    if (!quota) return 0;

    const available = quota.available_quota;
    const used = quota.account_details?.used_cpu_hours ?? 0;
    const total = available + used;

    if (total === 0) return 0;
    return Math.round((used / total) * 100);
  }, [quota]);

  const isQuotaSufficient = useCallback((requiredHours: number = 1): boolean => {
    return (quota?.available_quota ?? 0) >= requiredHours;
  }, [quota]);

  const getQuotaSourceBreakdown = useCallback((): QuotaSource => {
    return quota?.quota_sources ?? {};
  }, [quota]);

  return {
    quota,
    loading,
    error,
    refetch: fetchQuota,
    getAvailableQuota,
    getQuotaPercentage,
    isQuotaSufficient,
    getQuotaSourceBreakdown,
  };
};

/**
 * 格式化配额显示 - 最多两位小数
 */
export const formatQuota = (hours: number | null | undefined, decimals: number = 2): string => {
  if (hours === null || hours === undefined) {
    return `${'0.' + '0'.repeat(decimals)} h`;
  }
  return `${hours.toFixed(decimals)} h`;
};

/**
 * 获取配额状态颜色
 */
export const getQuotaStatusColor = (percentage: number): string => {
  if (percentage >= 90) return '#ff4d4f'; // 红色 - 即将用尽
  if (percentage >= 70) return '#faad14'; // 橙色 - 使用较多
  if (percentage >= 50) return '#1890ff'; // 蓝色 - 使用中等
  return '#52c41a'; // 绿色 - 充足
};

/**
 * 获取配额状态文本
 */
export const getQuotaStatusText = (percentage: number): string => {
  if (percentage >= 90) return '即将用尽';
  if (percentage >= 70) return '使用较多';
  if (percentage >= 50) return '使用中等';
  return '充足';
};

/**
 * 获取配额来源标签
 */
export const getQuotaSourceLabel = (source: string): string => {
  const labels: Record<string, string> = {
    personal_quota: '个人配额',
    organization_quota: '组织配额',
    master_account_quota: '主账号配额',
  };
  return labels[source] || source;
};

export default useQuota;

