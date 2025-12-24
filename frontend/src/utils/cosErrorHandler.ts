/**
 * COS（腾讯云对象存储）错误处理工具
 * 
 * 处理COS相关的错误，特别是账户欠费等问题
 */

import type { AxiosError } from 'axios';

export interface COSError {
  code: string;
  message: string;
  resource?: string;
  requestid?: string;
  traceid?: string;
}

/**
 * 检查是否是COS错误
 */
export function isCOSError(error: any): error is COSError {
  return error?.code && error?.message;
}

/**
 * 检查是否是账户欠费错误
 */
export function isAccountArrearsError(error: any): boolean {
  if (isCOSError(error)) {
    return error.code === 'UnavailableForLegalReasons' && 
           error.message?.includes('arrears');
  }
  return false;
}

/**
 * 获取用户友好的错误提示
 */
export function getCOSErrorMessage(error: any): string {
  if (!isCOSError(error)) {
    return '文件访问失败，请稍后重试';
  }

  switch (error.code) {
    case 'UnavailableForLegalReasons':
      if (error.message?.includes('arrears')) {
        return '账户欠费，请充值后重试';
      }
      return '资源暂时不可用，请稍后重试';
    
    case 'NoSuchKey':
      return '文件不存在';
    
    case 'AccessDenied':
      return '没有权限访问此文件';
    
    case 'InvalidBucketName':
      return '存储桶配置错误';
    
    case 'ServiceUnavailable':
      return '服务暂时不可用，请稍后重试';
    
    case 'RequestTimeout':
      return '请求超时，请检查网络连接';
    
    default:
      return `文件访问失败: ${error.message || '未知错误'}`;
  }
}

/**
 * 处理COS错误
 */
export function handleCOSError(error: any): {
  message: string;
  isArrearsError: boolean;
  isRetryable: boolean;
} {
  const message = getCOSErrorMessage(error);
  const isArrearsError = isAccountArrearsError(error);
  
  // 判断是否可重试
  const isRetryable = !isArrearsError && 
                      error?.code !== 'NoSuchKey' && 
                      error?.code !== 'AccessDenied';

  return {
    message,
    isArrearsError,
    isRetryable,
  };
}

/**
 * 从Axios错误中提取COS错误信息
 */
export function extractCOSErrorFromAxios(axiosError: AxiosError): COSError | null {
  const data = axiosError.response?.data;

  if (isCOSError(data)) {
    return data as COSError;
  }

  return null;
}

/**
 * 清除浏览器缓存
 * 用于解决充值后仍显示欠费错误的问题
 */
export function clearBrowserCache(): void {
  // 清除 localStorage 中的缓存数据
  const cacheKeys = [
    'job_detail_cache',
    'md_jobs_cache',
    'qc_jobs_cache',
    'electrolyte_cache',
    'cos_error_cache',
  ];

  cacheKeys.forEach(key => {
    localStorage.removeItem(key);
  });

  // 清除 sessionStorage
  sessionStorage.clear();

  console.log('✅ 浏览器缓存已清除');
}

/**
 * 强制刷新页面（清除所有缓存）
 * 用于解决充值后仍显示欠费错误的问题
 */
export function forceRefreshPage(): void {
  // 清除缓存
  clearBrowserCache();

  // 强制刷新页面（不使用缓存）
  window.location.reload();
}

