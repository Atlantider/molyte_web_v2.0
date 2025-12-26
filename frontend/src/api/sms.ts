import apiClient from './client';

export interface SendCodeRequest {
  phone: string;
  purpose?: 'register' | 'bind' | 'reset';
}

export interface VerifyCodeRequest {
  phone: string;
  code: string;
  purpose?: 'register' | 'bind' | 'reset';
}

export interface BindPhoneRequest {
  phone: string;
  code: string;
}

/**
 * 发送短信验证码
 */
export const sendSmsCode = async (data: SendCodeRequest) => {
  const response = await apiClient.post('/sms/send-code', data);
  return response.data;
};

/**
 * 验证短信验证码
 */
export const verifySmsCode = async (data: VerifyCodeRequest) => {
  const response = await apiClient.post('/sms/verify-code', data);
  return response.data;
};

/**
 * 绑定手机号到当前用户
 */
export const bindPhone = async (data: BindPhoneRequest) => {
  const response = await apiClient.post('/sms/bind-phone', data);
  return response.data;
};

