import apiClient from './client';

const API_BASE_URL = '/pricing';

export interface TaskTypePrice {
    task_type: string;
    price_per_hour: number;
}

export interface UserTypePrice {
    user_type: string;
    core_hour_price: number;
}

export interface UserDiscountInfo {
    discount_rate: number;
    discount_source: string;
    user_type?: string;
}

export interface PriceCalculationRequest {
    task_type: string;
    cpu_hours: number;
}

export interface PriceCalculationResponse {
    task_type: string;
    cpu_hours: number;
    task_price: number;
    discount_rate: number;
    final_price: number;
}

// 用户端点
export const getTaskTypePrices = async (): Promise<TaskTypePrice[]> => {
    const response = await apiClient.get(`${API_BASE_URL}/task-type-prices`);
    return response.data;
};

export const getMyDiscount = async (): Promise<UserDiscountInfo> => {
    const response = await apiClient.get(`${API_BASE_URL}/my-discount`);
    return response.data;
};

export const calculatePrice = async (
    request: PriceCalculationRequest
): Promise<PriceCalculationResponse> => {
    const response = await apiClient.post(`${API_BASE_URL}/calculate`, request);
    return response.data;
};

// 管理员端点 - 任务类型单价
export const adminGetTaskTypePrices = async (): Promise<TaskTypePrice[]> => {
    const response = await apiClient.get(`${API_BASE_URL}/admin/task-types`);
    return response.data;
};

export const adminUpdateTaskTypePrice = async (
    taskType: string,
    pricePerHour: number
): Promise<{ success: boolean; message: string }> => {
    const response = await apiClient.put(
        `${API_BASE_URL}/admin/task-types/${taskType}`,
        {
            price_per_hour: pricePerHour,
        }
    );
    return response.data;
};

// 管理员端点 - 用户类型核时单价
export const adminGetUserTypePrices = async (): Promise<UserTypePrice[]> => {
    const response = await apiClient.get(`${API_BASE_URL}/admin/user-type-prices`);
    return response.data;
};

export const adminUpdateUserTypePrice = async (
    userType: string,
    coreHourPrice: number
): Promise<{ success: boolean; message: string }> => {
    const response = await apiClient.put(
        `${API_BASE_URL}/admin/user-type-prices/${userType}`,
        {
            core_hour_price: coreHourPrice,
        }
    );
    return response.data;
};

export const adminSetUserCustomDiscount = async (
    userId: number,
    discountRate: number | null
): Promise<{ success: boolean; message: string }> => {
    const response = await apiClient.put(
        `${API_BASE_URL}/admin/users/${userId}/custom-discount`,
        {
            discount_rate: discountRate,
        }
    );
    return response.data;
};

// 计费配置相关
export interface BillingConfig {
    pricing_mode: 'CORE_HOUR' | 'TASK_TYPE';
}

export const getBillingConfig = async (): Promise<BillingConfig> => {
    const response = await apiClient.get(`${API_BASE_URL}/config`);
    return response.data;
};

export const adminGetBillingConfig = async (): Promise<BillingConfig> => {
    const response = await apiClient.get(`${API_BASE_URL}/admin/config`);
    return response.data;
};

export const adminUpdateBillingConfig = async (config: Partial<BillingConfig>): Promise<{ success: boolean; message: string }> => {
    const response = await apiClient.put(`${API_BASE_URL}/admin/config`, config);
    return response.data;
};
