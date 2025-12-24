/**
 * Admin API client
 */
import axios from 'axios';

const API_BASE_URL = '/api/v1';

// ============ Types ============

export interface UserListItem {
  id: number;
  username: string;
  email: string;
  role: 'ADMIN' | 'PREMIUM' | 'USER' | 'GUEST';
  user_type: 'STUDENT' | 'RESEARCHER' | 'COMPANY';
  organization: string | null;
  department: string | null;
  is_active: boolean;
  balance_cpu_hours: number;
  frozen_cpu_hours: number;
  free_cpu_hours_granted: number;
  recharge_cpu_hours: number;
  admin_granted_cpu_hours: number;
  daily_job_limit: number;
  concurrent_job_limit: number;
  storage_quota_gb: number;
  allowed_partitions: string[] | null;
  allowed_modules: string[] | null;
  last_login_at: string | null;
  created_at: string;
  custom_cpu_hour_price?: number | null;
}

export interface UserDetail extends UserListItem {
  updated_at: string;
  used_cpu_hours: number;
  today_jobs: number;
  running_jobs: number;
  total_jobs: number;
  completed_jobs: number;
  failed_jobs: number;
}

export interface UserUpdate {
  email?: string;
  role?: 'ADMIN' | 'PREMIUM' | 'USER' | 'GUEST';
  is_active?: boolean;
  balance_cpu_hours?: number;
  free_cpu_hours_granted?: number;
  recharge_cpu_hours?: number;
  admin_granted_cpu_hours?: number;
  daily_job_limit?: number;
  concurrent_job_limit?: number;
  storage_quota_gb?: number;
  allowed_partitions?: string[] | null;
  allowed_modules?: string[] | null;
  custom_cpu_hour_price?: number | null;
}

export interface UserCreate {
  username: string;
  email: string;
  password: string;
  role?: 'ADMIN' | 'PREMIUM' | 'USER' | 'GUEST';
  balance_cpu_hours?: number;
  free_cpu_hours_granted?: number;
  daily_job_limit?: number;
  concurrent_job_limit?: number;
  storage_quota_gb?: number;
  allowed_partitions?: string[] | null;
  allowed_modules?: string[] | null;
}

export interface GlobalStats {
  total_users: number;
  active_users: number;
  total_jobs: number;
  running_jobs: number;
  queued_jobs: number;
  completed_jobs: number;
  failed_jobs: number;
  total_cpu_hours_used: number;
  total_cpu_hours_allocated: number;
  total_storage_used_gb: number;
  total_storage_allocated_gb: number;
}

export interface UserUsageStatsItem {
  user_id: number;
  username: string;
  email: string;
  role: string;
  used_cpu_hours: number;
  total_cpu_hours: number;
  usage_percentage: number;
  total_jobs: number;
  running_jobs: number;
  completed_jobs: number;
  failed_jobs: number;
  last_job_at: string | null;
}

export interface UserRanking {
  user_id: number;
  username: string;
  email: string;
  metric_value: number;
  rank: number;
}

export interface AuditLogItem {
  id: number;
  user_id: number | null;
  username: string | null;
  action: string;
  resource_type: string | null;
  resource_id: number | null;
  details: any;
  ip_address: string | null;
  created_at: string;
}

export interface QuotaCheckResponse {
  allowed: boolean;
  reason?: string;
  details: any;
}

// ============ API Functions ============

const getAuthHeaders = () => {
  const token = localStorage.getItem('access_token');
  return {
    Authorization: `Bearer ${token}`,
  };
};

// User Management
export const getAllUsers = async (params?: {
  skip?: number;
  limit?: number;
  role?: string;
  is_active?: boolean;
}): Promise<UserListItem[]> => {
  const response = await axios.get(`${API_BASE_URL}/admin/users`, {
    headers: getAuthHeaders(),
    params,
  });
  return response.data;
};

export const getUserDetail = async (userId: number): Promise<UserDetail> => {
  const response = await axios.get(`${API_BASE_URL}/admin/users/${userId}`, {
    headers: getAuthHeaders(),
  });
  return response.data;
};

export const updateUser = async (userId: number, data: UserUpdate): Promise<UserDetail> => {
  const response = await axios.put(`${API_BASE_URL}/admin/users/${userId}`, data, {
    headers: getAuthHeaders(),
  });
  return response.data;
};

export const createUser = async (data: UserCreate): Promise<UserDetail> => {
  const response = await axios.post(`${API_BASE_URL}/admin/users`, data, {
    headers: getAuthHeaders(),
  });
  return response.data;
};

export const deleteUser = async (userId: number): Promise<void> => {
  await axios.delete(`${API_BASE_URL}/admin/users/${userId}`, {
    headers: getAuthHeaders(),
  });
};

export const updateUserStatus = async (userId: number, isActive: boolean): Promise<any> => {
  const response = await axios.put(
    `${API_BASE_URL}/admin/users/${userId}/status`,
    null,
    {
      headers: getAuthHeaders(),
      params: { is_active: isActive },
    }
  );
  return response.data;
};

export const checkUserQuota = async (userId: number): Promise<QuotaCheckResponse> => {
  const response = await axios.get(`${API_BASE_URL}/admin/users/${userId}/quota`, {
    headers: getAuthHeaders(),
  });
  return response.data;
};

// Statistics and Monitoring
export const getGlobalStats = async (): Promise<GlobalStats> => {
  const response = await axios.get(`${API_BASE_URL}/admin/stats/overview`, {
    headers: getAuthHeaders(),
  });
  return response.data;
};

export const getUserUsageStats = async (params?: {
  skip?: number;
  limit?: number;
  sort_by?: string;
}): Promise<UserUsageStatsItem[]> => {
  const response = await axios.get(`${API_BASE_URL}/admin/stats/users`, {
    headers: getAuthHeaders(),
    params,
  });
  return response.data;
};

export const getCpuUsageRanking = async (limit: number = 10): Promise<UserRanking[]> => {
  const response = await axios.get(`${API_BASE_URL}/admin/stats/ranking/cpu`, {
    headers: getAuthHeaders(),
    params: { limit },
  });
  return response.data;
};

export const getJobCountRanking = async (limit: number = 10): Promise<UserRanking[]> => {
  const response = await axios.get(`${API_BASE_URL}/admin/stats/ranking/jobs`, {
    headers: getAuthHeaders(),
    params: { limit },
  });
  return response.data;
};

// Job Management
export const getAllJobs = async (params?: {
  skip?: number;
  limit?: number;
  status?: string;
  user_id?: number;
}): Promise<any[]> => {
  const response = await axios.get(`${API_BASE_URL}/admin/jobs/all`, {
    headers: getAuthHeaders(),
    params,
  });
  return response.data;
};

export const adminCancelJob = async (jobId: number, reason?: string): Promise<any> => {
  const response = await axios.post(
    `${API_BASE_URL}/admin/jobs/${jobId}/cancel`,
    null,
    {
      headers: getAuthHeaders(),
      params: { reason },
    }
  );
  return response.data;
};

// Audit Logs
export const getAuditLogs = async (params?: {
  skip?: number;
  limit?: number;
  action?: string;
  user_id?: number;
  resource_type?: string;
  start_date?: string;
  end_date?: string;
}): Promise<AuditLogItem[]> => {
  const response = await axios.get(`${API_BASE_URL}/admin/logs`, {
    headers: getAuthHeaders(),
    params,
  });
  return response.data;
};

// Statistics Summary
export const getStatisticsSummary = async (period: string = '7days'): Promise<any> => {
  const response = await axios.get(`${API_BASE_URL}/admin/statistics/summary`, {
    headers: getAuthHeaders(),
    params: { period },
  });
  return response.data;
};

// ============ Partition Management ============

export interface PartitionInfo {
  name: string;
  state: string;
  total_nodes: number;
  available_nodes: number;
  total_cpus: number;
  available_cpus: number;
  max_time: string | null;
}

/**
 * Get all available partitions (admin only)
 * Returns all partitions without filtering for user assignment
 */
export const getAllPartitions = async (): Promise<PartitionInfo[]> => {
  const response = await axios.get(`${API_BASE_URL}/admin/partitions`, {
    headers: getAuthHeaders(),
  });
  return response.data;
};

// ============ Master Account Management ============

export interface MasterAccount {
  id: number;
  user_id: number;
  username?: string;
  email?: string;
  organization?: string;
  custom_cpu_hour_price?: number | null;
  // 主账号的配额来自 User 表
  balance_cpu_hours: number;
  frozen_cpu_hours: number;
  free_cpu_hours_granted: number;
  recharge_cpu_hours: number;
  admin_granted_cpu_hours: number;
  current_sub_accounts: number;
  max_sub_accounts: number;
  is_active: boolean;
  created_at: string;
  updated_at: string;
}

export interface SubAccount {
  id: number;
  master_account_id: number;
  user_id: number;
  username?: string;
  email?: string;
  // 子账号的配额来自两个来源：
  // 1. User.balance_cpu_hours - 子账号自己充值的余额
  // 2. SubAccount.allocated_quota - 主账号分配的配额
  // 实际可用 = min(User.balance_cpu_hours, SubAccount.allocated_quota)
  balance_cpu_hours?: number;  // 来自 User 表
  frozen_cpu_hours?: number;   // 来自 User 表
  allocated_quota: number;     // 主账号分配的配额
  is_active: boolean;
  created_at: string;
  updated_at: string;
}

export interface CreateMasterAccountRequest {
  // 方式1: 升级已有用户
  user_id?: number;
  // 方式2: 创建新用户
  username?: string;
  email?: string;
  password?: string;
  // 共用字段
  organization?: string;
  max_sub_accounts?: number;
  initial_cpu_hours?: number;
}

export interface CreateSubAccountRequest {
  username: string;
  email: string;
  password: string;
  allocated_quota?: number;  // 主账号分配的配额
}

/**
 * Get all master accounts (admin only)
 */
export const getAllMasterAccounts = async (params?: {
  skip?: number;
  limit?: number;
}): Promise<MasterAccount[]> => {
  const response = await axios.get(`${API_BASE_URL}/admin/master-accounts`, {
    headers: getAuthHeaders(),
    params,
  });
  return response.data;
};

/**
 * Get master account detail (admin only)
 */
export const getMasterAccountDetail = async (masterId: number): Promise<MasterAccount> => {
  const response = await axios.get(`${API_BASE_URL}/admin/master-accounts/${masterId}`, {
    headers: getAuthHeaders(),
  });
  return response.data;
};

/**
 * Get sub accounts of a master account (admin only)
 */
export const getSubAccounts = async (masterId: number): Promise<SubAccount[]> => {
  const response = await axios.get(`${API_BASE_URL}/admin/master-accounts/${masterId}/sub-accounts`, {
    headers: getAuthHeaders(),
  });
  return response.data;
};

/**
 * Create master account (admin only)
 */
export const createMasterAccount = async (data: CreateMasterAccountRequest): Promise<MasterAccount> => {
  const response = await axios.post(`${API_BASE_URL}/admin/master-accounts`, data, {
    headers: getAuthHeaders(),
  });
  return response.data;
};

/**
 * Create sub account (admin only)
 */
export const createSubAccount = async (masterId: number, data: CreateSubAccountRequest): Promise<SubAccount> => {
  const response = await axios.post(`${API_BASE_URL}/admin/master-accounts/${masterId}/sub-accounts`, data, {
    headers: getAuthHeaders(),
  });
  return response.data;
};

/**
 * Update sub account (admin only)
 */
export const updateSubAccount = async (
  masterId: number,
  subAccountId: number,
  data: { allocated_quota?: number; is_active?: boolean }
): Promise<SubAccount> => {
  const response = await axios.put(
    `${API_BASE_URL}/admin/master-accounts/${masterId}/sub-accounts/${subAccountId}`,
    null,
    {
      headers: getAuthHeaders(),
      params: data,
    }
  );
  return response.data;
};

/**
 * Delete sub account (admin only)
 */
export const deleteSubAccount = async (masterId: number, subAccountId: number): Promise<void> => {
  await axios.delete(`${API_BASE_URL}/admin/master-accounts/${masterId}/sub-accounts/${subAccountId}`, {
    headers: getAuthHeaders(),
  });
};

/**
 * Delete master account (admin only)
 * 将主账号转换为个人账号，所有子账号也转换为个人账号
 */
export const deleteMasterAccount = async (masterId: number): Promise<void> => {
  await axios.delete(`${API_BASE_URL}/admin/master-accounts/${masterId}`, {
    headers: getAuthHeaders(),
  });
};

