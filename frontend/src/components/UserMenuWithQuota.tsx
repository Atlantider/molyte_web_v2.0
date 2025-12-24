/**
 * UserMenuWithQuota 组件
 * 
 * 在用户菜单中显示账号类型和配额信息
 */

import React from 'react';
import { Space, Tag, Progress, Tooltip, Divider, Typography } from 'antd';
import { ThunderboltOutlined, CrownOutlined, TeamOutlined, UserOutlined, BankOutlined } from '@ant-design/icons';
import { useAccountType, getAccountTypeLabel, getAccountTypeColor } from '../hooks/useAccountType';
import { useQuota, formatQuota, getQuotaStatusColor, getQuotaStatusText } from '../hooks/useQuota';

const { Text } = Typography;

interface UserMenuWithQuotaProps {
  compact?: boolean;
}

/**
 * 获取账号类型图标
 */
const getAccountTypeIcon = (accountType: string) => {
  const icons: Record<string, React.ReactNode> = {
    PERSONAL: <UserOutlined />,
    MASTER_ACCOUNT: <CrownOutlined />,
    SUB_ACCOUNT: <BankOutlined />,
  };
  return icons[accountType] || <UserOutlined />;
};

const UserMenuWithQuota: React.FC<UserMenuWithQuotaProps> = ({ compact = false }) => {
  const { accountInfo, loading: accountLoading } = useAccountType();
  const { quota, loading: quotaLoading } = useQuota();

  if (accountLoading || quotaLoading) {
    return <Text type="secondary">加载中...</Text>;
  }

  if (!accountInfo || !quota) {
    return null;
  }

  const { account_type } = accountInfo;
  const label = getAccountTypeLabel(account_type);
  const color = getAccountTypeColor(account_type);
  const icon = getAccountTypeIcon(account_type);

  const available = quota.available_quota;
  const used = quota.account_details?.used_cpu_hours ?? 0;
  const total = available + used;
  const percentage = total > 0 ? Math.round((used / total) * 100) : 0;
  const statusColor = getQuotaStatusColor(percentage);
  const statusText = getQuotaStatusText(percentage);

  if (compact) {
    // 紧凑模式 - 用于导航栏
    return (
      <Tooltip
        title={
          <div style={{ fontSize: '12px' }}>
            <div style={{ marginBottom: '8px' }}>
              <strong>账号类型：</strong> {label}
            </div>
            <div style={{ marginBottom: '8px' }}>
              <strong>可用配额：</strong> {formatQuota(available)}
            </div>
            <div style={{ marginBottom: '8px' }}>
              <strong>已用配额：</strong> {formatQuota(used)}
            </div>
            <div>
              <strong>使用情况：</strong> {percentage}% ({statusText})
            </div>
          </div>
        }
      >
        <Space size="small" style={{ cursor: 'pointer' }}>
          <Tag color={color} icon={icon} style={{ margin: 0 }}>
            {label}
          </Tag>
          <Tag
            color={statusColor === '#ff4d4f' ? 'red' : statusColor === '#faad14' ? 'orange' : statusColor === '#1890ff' ? 'blue' : 'green'}
            icon={<ThunderboltOutlined />}
            style={{ margin: 0 }}
          >
            {formatQuota(available)}
          </Tag>
        </Space>
      </Tooltip>
    );
  }

  // 完整模式 - 用于下拉菜单
  return (
    <div style={{ padding: '12px 0' }}>
      <div style={{ marginBottom: '12px' }}>
        <div style={{ marginBottom: '8px', fontSize: '12px', color: '#666' }}>
          <strong>账号类型</strong>
        </div>
        <Tag color={color} icon={icon}>
          {label}
        </Tag>
      </div>

      <Divider style={{ margin: '8px 0' }} />

      <div style={{ marginBottom: '12px' }}>
        <div style={{ marginBottom: '8px', fontSize: '12px', color: '#666' }}>
          <strong>配额信息</strong>
        </div>
        <div style={{ marginBottom: '8px', fontSize: '12px' }}>
          <Space>
            <span>可用：{formatQuota(available)}</span>
            <span>已用：{formatQuota(used)}</span>
          </Space>
        </div>
        <Progress
          percent={percentage}
          strokeColor={statusColor}
          status={percentage >= 90 ? 'exception' : 'active'}
          size="small"
        />
        <div style={{ marginTop: '4px', fontSize: '12px', color: '#666' }}>
          {statusText} ({percentage}%)
        </div>
      </div>

      {/* 配额来源 */}
      {quota.quota_sources && Object.keys(quota.quota_sources).length > 0 && (
        <>
          <Divider style={{ margin: '8px 0' }} />
          <div style={{ fontSize: '12px', color: '#666' }}>
            <strong>配额来源</strong>
            {Object.entries(quota.quota_sources).map(([key, value]) => (
              value !== undefined && value !== null && (
                <div key={key} style={{ marginTop: '4px' }}>
                  <span>{key === 'personal_quota' ? '个人' : key === 'organization_quota' ? '组织' : '主账号'}：</span>
                  <span style={{ fontWeight: 'bold' }}>{formatQuota(value)}</span>
                </div>
              )
            ))}
          </div>
        </>
      )}

      {/* 配额不足警告 */}
      {available < 1 && (
        <>
          <Divider style={{ margin: '8px 0' }} />
          <div style={{ fontSize: '12px', color: '#ff4d4f', fontWeight: 'bold' }}>
            ⚠️ 配额即将用尽，请及时充值
          </div>
        </>
      )}
    </div>
  );
};

export default UserMenuWithQuota;

