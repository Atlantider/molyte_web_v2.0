/**
 * 账户管理 API
 * 
 * 处理主账号和子账号的相关操作
 */

import client from './client';

// 子账号类型定义
export interface SubAccount {
  id: number;
  master_account_id: number;
  user_id: number;
  username?: string;
  email?: string;
  master_username?: string;
  master_email?: string;
  personal_quota?: number;
  personal_used?: number;
  balance_cpu_hours: number;           // 子账号自己充值的余额
  frozen_cpu_hours?: number;           // 冻结核时
  allocated_quota: number;              // 主账号分配的配额
  is_active: boolean;
  created_at: string;
  updated_at: string;
}

// 创建子账号请求
export interface CreateSubAccountRequest {
  username: string;
  email: string;
  password: string;
  personal_quota?: number;
}

// 更新子账号请求
export interface UpdateSubAccountRequest {
  personal_quota?: number;
  is_active?: boolean;
}

// 主账号信息
export interface MasterAccountInfo {
  id: number;
  user_id: number;
  total_cpu_hours: number;
  balance_cpu_hours: number;
  frozen_cpu_hours: number;
  used_cpu_hours: number;
  current_sub_accounts: number;
  max_sub_accounts: number;
  is_active: boolean;
  created_at: string;
  updated_at: string;
  username?: string;
  email?: string;
  organization?: string;
}

/**
 * 获取我的子账号列表
 */
export const getMySubAccounts = async (): Promise<SubAccount[]> => {
  const response = await client.get('/accounts/my-sub-accounts');
  return response.data;
};

/**
 * 创建子账号
 */
export const createMySubAccount = async (data: CreateSubAccountRequest): Promise<SubAccount> => {
  const response = await client.post('/accounts/my-sub-accounts', data);
  return response.data;
};

/**
 * 更新子账号
 */
export const updateMySubAccount = async (subAccountId: number, data: UpdateSubAccountRequest): Promise<SubAccount> => {
  const response = await client.put(`/accounts/my-sub-accounts/${subAccountId}`, data);
  return response.data;
};

/**
 * 删除子账号
 */
export const deleteMySubAccount = async (subAccountId: number): Promise<void> => {
  await client.delete(`/accounts/my-sub-accounts/${subAccountId}`);
};

/**
 * 获取主账号信息
 */
export const getMyMasterAccountInfo = async (): Promise<MasterAccountInfo> => {
  const response = await client.get('/accounts/my-master-account');
  return response.data;
};

/**
 * 获取子账号信息
 */
export const getMySubAccountInfo = async (): Promise<SubAccount> => {
  const response = await client.get('/accounts/my-sub-account');
  return response.data;
};

/**
 * 添加现有用户为子账号
 */
export interface AddExistingUserRequest {
  username_or_email: string;
  allocated_quota?: number;
}

export const addExistingUserAsSubAccount = async (data: AddExistingUserRequest): Promise<SubAccount> => {
  const response = await client.post('/accounts/my-sub-accounts/add-existing', data);
  return response.data;
};

/**
 * 获取子账号的任务列表
 */
export interface SubAccountJob {
  id: number;
  system_id: number;
  status: string;
  username?: string;
  user_id?: number;
  created_at: string;
  started_at?: string;
  finished_at?: string;
  actual_cpu_hours?: number;
  estimated_cpu_hours?: number;
}

export interface SubAccountJobsResponse {
  sub_account_id?: number;
  username?: string;
  total: number;
  skip: number;
  limit: number;
  jobs: SubAccountJob[];
}

export const getSubAccountJobs = async (subAccountId: number, skip = 0, limit = 20): Promise<SubAccountJobsResponse> => {
  const response = await client.get(`/accounts/my-sub-accounts/${subAccountId}/jobs`, {
    params: { skip, limit }
  });
  return response.data;
};

/**
 * 获取所有子账号的任务汇总
 */
export const getAllSubAccountsJobs = async (skip = 0, limit = 50): Promise<SubAccountJobsResponse> => {
  const response = await client.get('/accounts/my-sub-accounts/all-jobs', {
    params: { skip, limit }
  });
  return response.data;
};

