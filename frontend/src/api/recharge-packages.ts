/**
 * 充值套餐 API
 */
import client from './client';

export interface RechargePackage {
  id: number;
  name: string;
  description?: string;
  user_type: string;
  price: number;
  cpu_hours: number;
  display_order: number;
  badge?: string;
  color: string;
  icon?: string;
  is_active: boolean;
}

/**
 * 获取充值套餐列表
 * @param userType 用户类型（可选）
 */
export const getRechargePackages = async (userType?: string): Promise<RechargePackage[]> => {
  const params = new URLSearchParams();
  if (userType) {
    params.append('user_type', userType);
  }
  
  const response = await client.get(`/recharge-packages/?${params.toString()}`);
  return response.data;
};

/**
 * 获取单个充值套餐
 * @param packageId 套餐ID
 */
export const getRechargePackage = async (packageId: number): Promise<RechargePackage> => {
  const response = await client.get(`/recharge-packages/${packageId}`);
  return response.data;
};

/**
 * 创建充值套餐（管理员）
 */
export const createRechargePackage = async (data: Omit<RechargePackage, 'id' | 'is_active'>): Promise<RechargePackage> => {
  const response = await client.post('/recharge-packages/', data);
  return response.data;
};

/**
 * 更新充值套餐（管理员）
 */
export const updateRechargePackage = async (packageId: number, data: Partial<RechargePackage>): Promise<RechargePackage> => {
  const response = await client.put(`/recharge-packages/${packageId}`, data);
  return response.data;
};

/**
 * 删除充值套餐（管理员）
 */
export const deleteRechargePackage = async (packageId: number): Promise<void> => {
  await client.delete(`/recharge-packages/${packageId}`);
};

/**
 * 获取所有充值套餐（管理员）
 */
export const getAllRechargePackages = async (userType?: string): Promise<RechargePackage[]> => {
  const params = new URLSearchParams();
  if (userType) {
    params.append('user_type', userType);
  }
  
  const response = await client.get(`/recharge-packages/admin/all?${params.toString()}`);
  return response.data;
};

