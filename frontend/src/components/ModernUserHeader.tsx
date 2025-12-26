/**
 * 现代化用户头部组件
 * 
 * 采用现代 SaaS 应用的设计风格：
 * - 简洁的配额显示
 * - 优雅的用户下拉菜单
 * - 清晰的信息层次
 */

import React, { useState, useEffect } from 'react';
import { Space, Dropdown, Avatar, Typography, Tag, Button, Tooltip, Badge } from 'antd';
import {
  UserOutlined,
  LogoutOutlined,
  SettingOutlined,
  WalletOutlined,
  ThunderboltOutlined,
  CrownOutlined,
  TeamOutlined,
  BankOutlined,
  DashboardOutlined,
  BellOutlined,
  WarningOutlined,
} from '@ant-design/icons';
import { useNavigate, useLocation, useSearchParams } from 'react-router-dom';
import { useAuthStore } from '../stores/authStore';
import { useAccountType } from '../hooks/useAccountType';
import { useQuota } from '../hooks/useQuota';
import { useThemeStore } from '../stores/themeStore';
import { getUnreadCount } from '../api/notifications';
import type { MenuProps } from 'antd';

const { Text } = Typography;

interface ModernUserHeaderProps {
  onLogout: () => void;
}

/**
 * 获取账号类型图标和颜色
 */
const getAccountTypeConfig = (accountType: string) => {
  const configs = {
    personal: { icon: <UserOutlined />, color: '#1890ff', label: '个人' },
    master_account: { icon: <CrownOutlined />, color: '#faad14', label: '主账号' },
    sub_account: { icon: <BankOutlined />, color: '#722ed1', label: '子账号' },
  };
  return configs[accountType as keyof typeof configs] || configs.personal;
};

/**
 * 格式化配额显示 - 最多两位小数
 */
const formatQuotaDisplay = (quota: number | undefined): string => {
  if (!quota && quota !== 0) {
    return '0.00';
  }
  if (quota >= 1000) {
    return `${(quota / 1000).toFixed(2)}k`;
  }
  return quota.toFixed(2);
};

const ModernUserHeader: React.FC<ModernUserHeaderProps> = ({ onLogout }) => {
  const navigate = useNavigate();
  const location = useLocation();
  const [searchParams, setSearchParams] = useSearchParams();
  const { user } = useAuthStore();
  const { accountInfo, loading: accountLoading } = useAccountType();
  const { quota, loading: quotaLoading } = useQuota();
  const { mode } = useThemeStore();
  const isDark = mode === 'dark';
  const [unreadCount, setUnreadCount] = useState(0);

  if (!user) return null;

  // 判断是否在账户中心页面
  const isInAccountCenter = location.pathname === '/workspace/account-center';

  // 导航到账户中心的特定 tab
  const navigateToAccountTab = (tab?: string) => {
    if (isInAccountCenter) {
      // 已经在账户中心，直接更新 URL 参数
      if (tab) {
        setSearchParams({ tab });
      } else {
        setSearchParams({});
      }
    } else {
      // 不在账户中心，使用 navigate
      if (tab) {
        navigate(`/workspace/account-center?tab=${tab}`);
      } else {
        navigate('/workspace/account-center');
      }
    }
  };

  // 定期获取未读消息数
  useEffect(() => {
    const fetchUnreadCount = async () => {
      try {
        const result = await getUnreadCount();
        setUnreadCount(result.unread_count);
      } catch (error) {
        console.error('获取未读消息数失败:', error);
      }
    };

    fetchUnreadCount();
    const interval = setInterval(fetchUnreadCount, 30000); // 每30秒更新一次
    return () => clearInterval(interval);
  }, []);

  const loading = accountLoading || quotaLoading;
  const available = quota?.available_quota || 0;
  const isLowQuota = available < 10;
  const isCriticalQuota = available < 1;

  const accountConfig = accountInfo ? getAccountTypeConfig(accountInfo.account_type) : null;

  // 用户下拉菜单
  const userMenuItems: MenuProps['items'] = [
    // 用户信息头部
    {
      key: 'user-info',
      label: (
        <div style={{ padding: '12px 0 8px 0', borderBottom: '1px solid #f0f0f0' }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: 12, marginBottom: 8 }}>
            <Avatar size={40} icon={<UserOutlined />} style={{ 
              background: 'linear-gradient(135deg, #5B8DEF 0%, #7C6EAF 100%)' 
            }} />
            <div>
              <div style={{ fontWeight: 600, fontSize: 14, marginBottom: 2 }}>
                {user.username}
              </div>
              <div style={{ fontSize: 12, color: '#666' }}>
                {user.email}
              </div>
            </div>
          </div>
          
          {/* 账号类型 */}
          {accountConfig && (
            <div style={{ marginBottom: 8 }}>
              <Tag 
                icon={accountConfig.icon} 
                color={accountConfig.color}
                style={{ fontSize: 12 }}
              >
                {accountConfig.label}
              </Tag>
            </div>
          )}

          {/* 配额信息 */}
          <div style={{
            background: isCriticalQuota
              ? (isDark ? '#2a1215' : '#fff2f0')
              : isLowQuota
                ? (isDark ? '#2b2611' : '#fffbe6')
                : (isDark ? '#162312' : '#f6ffed'),
            padding: '8px 12px',
            borderRadius: 6,
            border: `1px solid ${isCriticalQuota
              ? (isDark ? '#5c2c2c' : '#ffccc7')
              : isLowQuota
                ? (isDark ? '#594214' : '#ffe58f')
                : (isDark ? '#274916' : '#d9f7be')}`,
          }}>
            <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
              <span style={{ fontSize: 12, color: '#666' }}>可用配额</span>
              {(isLowQuota || isCriticalQuota) && (
                <WarningOutlined style={{ 
                  color: isCriticalQuota ? '#ff4d4f' : '#faad14',
                  fontSize: 12 
                }} />
              )}
            </div>
            <div style={{ 
              fontSize: 16, 
              fontWeight: 600, 
              color: isCriticalQuota ? '#ff4d4f' : isLowQuota ? '#faad14' : '#52c41a',
              marginTop: 2
            }}>
              {loading ? '...' : `${formatQuotaDisplay(available)} 核时`}
            </div>
          </div>
        </div>
      ),
      disabled: true,
    },
    
    // 功能菜单
    {
      key: 'account-dashboard',
      icon: <DashboardOutlined />,
      label: '账户总览',
    },
    {
      key: 'recharge',
      icon: <WalletOutlined />,
      label: '充值中心',
    },
    {
      key: 'settings',
      icon: <SettingOutlined />,
      label: '账户设置',
    },
    {
      type: 'divider',
    },
    {
      key: 'logout',
      icon: <LogoutOutlined />,
      label: '退出登录',
      danger: true,
    },
  ];

  const handleMenuClick = ({ key }: { key: string }) => {
    switch (key) {
      case 'account-dashboard':
        navigateToAccountTab('overview');
        break;
      case 'recharge':
        navigateToAccountTab('recharge');
        break;
      case 'settings':
        navigateToAccountTab('settings');
        break;
      case 'logout':
        onLogout();
        break;
    }
  };

  return (
    <Space size={16}>
      {/* 配额快速显示 */}
      <Tooltip title={isCriticalQuota ? `欠费 ${Math.abs(available || 0).toFixed(2)} 核时，请立即充值` : `可用配额：${(available || 0).toFixed(2)} 核时`}>
        <div style={{
          display: 'flex',
          alignItems: 'center',
          gap: 6,
          padding: '4px 12px',
          background: isCriticalQuota
            ? (isDark ? '#2a1215' : '#fff2f0')
            : isLowQuota
              ? (isDark ? '#2b2611' : '#fffbe6')
              : (isDark ? '#162312' : '#f6ffed'),
          borderRadius: 16,
          border: `1px solid ${isCriticalQuota
            ? (isDark ? '#5c2c2c' : '#ffccc7')
            : isLowQuota
              ? (isDark ? '#594214' : '#ffe58f')
              : (isDark ? '#274916' : '#d9f7be')}`,
          cursor: 'pointer',
        }}
        onClick={() => navigateToAccountTab('recharge')}
        >
          <ThunderboltOutlined style={{
            color: isCriticalQuota ? '#ff4d4f' : isLowQuota ? '#faad14' : '#52c41a',
            fontSize: 14
          }} />
          <Text style={{
            fontSize: 13,
            fontWeight: 500,
            color: isCriticalQuota ? '#ff4d4f' : isLowQuota ? '#faad14' : '#52c41a',
          }}>
            {loading ? '...' : isCriticalQuota ? `欠费 ${Math.abs(available).toFixed(2)}` : `${formatQuotaDisplay(available)}`}
          </Text>
        </div>
      </Tooltip>

      {/* 通知按钮 */}
      <Tooltip title={unreadCount > 0 ? `${unreadCount} 条未读消息` : '消息中心'}>
        <Badge
          count={unreadCount > 0 ? unreadCount : 0}
          size="small"
          style={{
            backgroundColor: unreadCount > 0 ? '#faad14' : 'transparent',
            color: '#fff',
            boxShadow: unreadCount > 0 ? '0 0 0 1px #fff' : 'none',
          }}
        >
          <Button
            type="text"
            icon={<BellOutlined />}
            style={{
              width: 36,
              height: 36,
              borderRadius: 8,
              color: unreadCount > 0 ? '#faad14' : 'inherit'
            }}
            onClick={() => navigate('/workspace/notifications')}
          />
        </Badge>
      </Tooltip>

      {/* 用户下拉菜单 */}
      <Dropdown 
        menu={{ items: userMenuItems, onClick: handleMenuClick }} 
        placement="bottomRight"
        trigger={['click']}
      >
        <div style={{
          display: 'flex',
          alignItems: 'center',
          gap: 8,
          padding: '4px 4px 4px 8px',
          borderRadius: 20,
          background: isDark ? '#2a2a2a' : '#f5f7fa',
          cursor: 'pointer',
          transition: 'all 0.2s',
        }}
        onMouseEnter={(e) => {
          e.currentTarget.style.background = isDark ? '#3a3a3a' : '#e8ecf3';
        }}
        onMouseLeave={(e) => {
          e.currentTarget.style.background = isDark ? '#2a2a2a' : '#f5f7fa';
        }}
        >
          <Text style={{ fontSize: 13, fontWeight: 500, maxWidth: 80, overflow: 'hidden', textOverflow: 'ellipsis' }}>
            {user.username}
          </Text>
          <Avatar size={28} icon={<UserOutlined />} style={{ 
            background: 'linear-gradient(135deg, #5B8DEF 0%, #7C6EAF 100%)' 
          }} />
        </div>
      </Dropdown>
    </Space>
  );
};

export default ModernUserHeader;
