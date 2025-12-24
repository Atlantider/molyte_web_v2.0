/**
 * COS错误处理组件
 * 
 * 用于处理和显示COS相关的错误，特别是账户欠费错误
 * 提供清除缓存和重试的选项
 */

import { Alert, Button, Space, message } from 'antd';
import { ReloadOutlined, DeleteOutlined } from '@ant-design/icons';
import { handleCOSError, clearBrowserCache, forceRefreshPage } from '../utils/cosErrorHandler';

interface COSErrorHandlerProps {
  error: any;
  onRetry?: () => void;
  showClearCache?: boolean;
}

export default function COSErrorHandler({
  error,
  onRetry,
  showClearCache = true,
}: COSErrorHandlerProps) {
  const { message: errorMessage, isArrearsError, isRetryable } = handleCOSError(error);

  const handleClearCache = () => {
    clearBrowserCache();
    message.success('缓存已清除，请重新加载页面');
    
    // 2秒后自动刷新
    setTimeout(() => {
      window.location.reload();
    }, 2000);
  };

  const handleForceRefresh = () => {
    message.loading('正在清除缓存并刷新页面...');
    setTimeout(() => {
      forceRefreshPage();
    }, 1000);
  };

  return (
    <Alert
      message={isArrearsError ? '账户欠费' : '文件访问错误'}
      description={
        <div>
          <p>{errorMessage}</p>
          {isArrearsError && (
            <p style={{ marginTop: 8, fontSize: 12, color: 'rgba(0,0,0,0.65)' }}>
              💡 提示：如果您已经充值，请清除浏览器缓存后重试
            </p>
          )}
        </div>
      }
      type={isArrearsError ? 'error' : 'warning'}
      showIcon
      action={
        <Space size="small">
          {isRetryable && onRetry && (
            <Button
              size="small"
              type="primary"
              icon={<ReloadOutlined />}
              onClick={onRetry}
            >
              重试
            </Button>
          )}
          {showClearCache && (
            <>
              <Button
                size="small"
                icon={<DeleteOutlined />}
                onClick={handleClearCache}
              >
                清除缓存
              </Button>
              <Button
                size="small"
                type="primary"
                danger
                onClick={handleForceRefresh}
              >
                强制刷新
              </Button>
            </>
          )}
        </Space>
      }
      style={{ marginBottom: 16 }}
    />
  );
}

