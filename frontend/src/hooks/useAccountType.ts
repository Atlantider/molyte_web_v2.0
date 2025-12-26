/**
 * useAccountType Hook
 *
 * 用于获取和管理用户账号类型信息
 * 支持 3 种账号类型：PERSONAL, MASTER_ACCOUNT, SUB_ACCOUNT
 */

import { useState, useEffect } from 'react';
import { message } from 'antd';
import client from '../api/client';

export interface AccountInfo {
  account_type: 'PERSONAL' | 'MASTER_ACCOUNT' | 'SUB_ACCOUNT' | 'personal' | 'master_account' | 'sub_account';
  account_info: {
    account_type?: string;
    user_id?: number;
    username?: string;
    email?: string;
    total_cpu_hours?: number;
    balance_cpu_hours?: number;
    used_cpu_hours?: number;
    frozen_cpu_hours?: number;
    max_sub_accounts?: number;
    current_sub_accounts?: number;
    personal_quota?: number;
    personal_used?: number;
    quota_source?: string;
    [key: string]: any;
  };
  permissions: {
    can_create_jobs?: boolean;
    can_manage_master_account?: boolean;
    can_manage_sub_accounts?: boolean;
    can_access_admin?: boolean;
    [key: string]: any;
  };
}

export interface QuotaInfo {
  available_quota: number;
  quota_sources: {
    personal_quota?: number;
    organization_quota?: number;
    master_account_quota?: number;
  };
  account_type: string;
  account_details: {
    [key: string]: any;
  };
}

interface UseAccountTypeReturn {
  accountInfo: AccountInfo | null;
  quotaInfo: QuotaInfo | null;
  loading: boolean;
  error: string | null;
  refetch: () => Promise<void>;
}

/**
 * 获取用户账号类型和信息
 */
export const useAccountType = (): UseAccountTypeReturn => {
  const [accountInfo, setAccountInfo] = useState<AccountInfo | null>(null);
  const [quotaInfo, setQuotaInfo] = useState<QuotaInfo | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  const fetchAccountInfo = async () => {
    try {
      setLoading(true);
      setError(null);

      // 并行获取账号信息和配额信息
      const [accountResponse, quotaResponse] = await Promise.all([
        client.get('/users/me/account-info'),
        client.get('/users/me/quota'),
      ]);

      setAccountInfo(accountResponse.data);
      setQuotaInfo(quotaResponse.data);
    } catch (err: any) {
      const errorMsg = err.response?.data?.detail || '获取账号信息失败';
      setError(errorMsg);
      message.error(errorMsg);
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    fetchAccountInfo();
  }, []);

  return {
    accountInfo,
    quotaInfo,
    loading,
    error,
    refetch: fetchAccountInfo,
  };
};

/**
 * 获取账号类型的显示名称
 */
export const getAccountTypeLabel = (accountType: string): string => {
  const labels: Record<string, string> = {
    PERSONAL: '个人账号',
    personal: '个人账号',
    MASTER_ACCOUNT: '主账号',
    master_account: '主账号',
    SUB_ACCOUNT: '子账号',
    sub_account: '子账号',
  };
  return labels[accountType] || accountType;
};

/**
 * 获取账号类型的颜色
 */
export const getAccountTypeColor = (accountType: string): string => {
  const colors: Record<string, string> = {
    PERSONAL: 'blue',
    personal: 'blue',
    MASTER_ACCOUNT: 'gold',
    master_account: 'gold',
    SUB_ACCOUNT: 'purple',
    sub_account: 'purple',
  };
  return colors[accountType] || 'default';
};

/**
 * 检查用户是否有特定权限
 */
export const hasPermission = (accountInfo: AccountInfo | null, permission: string): boolean => {
  if (!accountInfo) return false;
  return accountInfo.permissions[permission] ?? false;
};

/**
 * 获取可用配额（小时）
 */
export const getAvailableQuota = (quotaInfo: QuotaInfo | null): number => {
  if (!quotaInfo) return 0;
  return quotaInfo.available_quota;
};

/**
 * 获取配额使用百分比
 */
export const getQuotaUsagePercentage = (quotaInfo: QuotaInfo | null, accountInfo: AccountInfo | null): number => {
  if (!quotaInfo || !accountInfo) return 0;

  const available = quotaInfo.available_quota;
  const used = accountInfo.account_info.used_cpu_hours || 0;
  const total = available + used;

  if (total === 0) return 0;
  return Math.round((used / total) * 100);
};

/**
 * 检查配额是否充足
 */
export const isQuotaSufficient = (quotaInfo: QuotaInfo | null, requiredHours: number = 1): boolean => {
  if (!quotaInfo) return false;
  return quotaInfo.available_quota >= requiredHours;
};

export default useAccountType;

