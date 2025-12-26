/**
 * 前端统一错误处理工具
 * 
 * 提供错误处理、错误消息转换和错误展示功能
 */

/**
 * API错误接口
 */
export interface APIError {
    code: string;
    message: string;
    details?: Record<string, any>;
}

/**
 * 错误码枚举
 */
export enum ErrorCode {
    // 客户端错误
    INVALID_INPUT = 'INVALID_INPUT',
    INVALID_SMILES = 'INVALID_SMILES',
    NOT_FOUND = 'NOT_FOUND',
    JOB_NOT_FOUND = 'JOB_NOT_FOUND',
    UNAUTHORIZED = 'UNAUTHORIZED',
    FORBIDDEN = 'FORBIDDEN',
    QUOTA_EXCEEDED = 'QUOTA_EXCEEDED',
    INSUFFICIENT_BALANCE = 'INSUFFICIENT_BALANCE',
    JOB_ALREADY_SUBMITTED = 'JOB_ALREADY_SUBMITTED',

    // 服务端错误
    INTERNAL_ERROR = 'INTERNAL_ERROR',
    DATABASE_ERROR = 'DATABASE_ERROR',
    SLURM_ERROR = 'SLURM_ERROR',

    // 未知错误
    UNKNOWN_ERROR = 'UNKNOWN_ERROR',
}

/**
 * 用户友好的错误消息映射
 */
const ERROR_MESSAGES: Record<string, string> = {
    // 客户端错误
    [ErrorCode.INVALID_INPUT]: '输入参数不正确，请检查后重试',
    [ErrorCode.INVALID_SMILES]: 'SMILES格式不正确，请检查分子结构',
    [ErrorCode.NOT_FOUND]: '请求的资源不存在',
    [ErrorCode.JOB_NOT_FOUND]: '任务不存在或已被删除',
    [ErrorCode.UNAUTHORIZED]: '请先登录',
    [ErrorCode.FORBIDDEN]: '您没有权限执行此操作',
    [ErrorCode.QUOTA_EXCEEDED]: '配额不足，请充值或联系管理员',
    [ErrorCode.INSUFFICIENT_BALANCE]: '余额不足，请充值',
    [ErrorCode.JOB_ALREADY_SUBMITTED]: '任务已提交，无法重复提交',

    // 服务端错误
    [ErrorCode.INTERNAL_ERROR]: '服务器错误，请稍后重试',
    [ErrorCode.DATABASE_ERROR]: '数据库操作失败，请稍后重试',
    [ErrorCode.SLURM_ERROR]: '集群调度失败，请联系管理员',

    // 未知错误
    [ErrorCode.UNKNOWN_ERROR]: '操作失败，请稍后重试',
};

/**
 * 从axios错误中提取API错误
 */
export function extractAPIError(error: any): APIError {
    // 检查是否是axios错误响应
    if (error.response?.data) {
        const data = error.response.data;

        // 检查是否是标准API错误格式
        if (data.code && data.message) {
            return {
                code: data.code,
                message: data.message,
                details: data.details || {},
            };
        }

        // 兼容旧的错误格式
        if (data.detail) {
            // Pydantic验证错误
            if (Array.isArray(data.detail)) {
                return {
                    code: ErrorCode.INVALID_INPUT,
                    message: '输入参数验证失败',
                    details: { validation_errors: data.detail },
                };
            }

            // 字符串详情
            if (typeof data.detail === 'string') {
                return {
                    code: ErrorCode.UNKNOWN_ERROR,
                    message: data.detail,
                };
            }
        }
    }

    // 网络错误
    if (error.message === 'Network Error') {
        return {
            code: ErrorCode.INTERNAL_ERROR,
            message: '网络连接失败，请检查网络后重试',
        };
    }

    // 超时错误
    if (error.code === 'ECONNABORTED') {
        return {
            code: ErrorCode.INTERNAL_ERROR,
            message: '请求超时，请稍后重试',
        };
    }

    // 未知错误
    return {
        code: ErrorCode.UNKNOWN_ERROR,
        message: error.message || '操作失败，请稍后重试',
    };
}

/**
 * 获取用户友好的错误消息
 */
export function getUserFriendlyMessage(error: APIError): string {
    return ERROR_MESSAGES[error.code] || error.message || '操作失败';
}

/**
 * 显示错误通知
 * 
 * 注意: 需要先导入 antd 的 message 或 notification
 */
export function showErrorNotification(error: any, customMessage?: string) {
    const apiError = extractAPIError(error);
    const message = customMessage || getUserFriendlyMessage(apiError);

    // 使用 antd 的 message 组件
    // 调用方需要先导入: import { message } from 'antd';
    return {
        message,
        description: apiError.details?.validation_errors
            ? '请检查输入参数是否正确'
            : undefined,
        error: apiError,
    };
}

/**
 * 错误处理装饰器（用于async函数）
 */
export function handleErrors<T extends (...args: any[]) => Promise<any>>(
    fn: T,
    errorHandler?: (error: APIError) => void
): T {
    return (async (...args: any[]) => {
        try {
            return await fn(...args);
        } catch (error) {
            const apiError = extractAPIError(error);

            if (errorHandler) {
                errorHandler(apiError);
            } else {
                console.error('API Error:', apiError);
            }

            throw apiError;
        }
    }) as T;
}

/**
 * React Hook: 错误处理
 */
export function useErrorHandler() {
    const handleError = (error: any, customMessage?: string) => {
        const apiError = extractAPIError(error);
        const message = customMessage || getUserFriendlyMessage(apiError);

        return {
            apiError,
            message,
            showNotification: (type: 'error' | 'warning' = 'error') => {
                // 返回通知配置，由调用方决定如何显示
                return {
                    type,
                    message,
                    description: apiError.details?.validation_errors
                        ? '请检查输入参数'
                        : undefined,
                };
            },
        };
    };

    return { handleError };
}

/**
 * 判断是否是认证错误
 */
export function isAuthError(error: APIError): boolean {
    return error.code === ErrorCode.UNAUTHORIZED || error.code === ErrorCode.FORBIDDEN;
}

/**
 * 判断是否是配额错误
 */
export function isQuotaError(error: APIError): boolean {
    return error.code === ErrorCode.QUOTA_EXCEEDED || error.code === ErrorCode.INSUFFICIENT_BALANCE;
}

/**
 * 判断是否是服务器错误
 */
export function isServerError(error: APIError): boolean {
    return [
        ErrorCode.INTERNAL_ERROR,
        ErrorCode.DATABASE_ERROR,
        ErrorCode.SLURM_ERROR,
    ].includes(error.code as ErrorCode);
}
