/**
 * 主布局组件
 */
import { useState, useEffect } from 'react';
import { Outlet, useNavigate, useLocation } from 'react-router-dom';
import { Layout as AntLayout, Menu, Dropdown, Avatar, Space, Button, Badge, Typography, Tooltip, Divider, Tag } from 'antd';
import {
  DashboardOutlined,
  ProjectOutlined,
  ExperimentOutlined,
  RocketOutlined,
  UserOutlined,
  LogoutOutlined,
  SettingOutlined,
  ControlOutlined,
  DatabaseOutlined,
  WalletOutlined,
  BellOutlined,
  ThunderboltOutlined,
  EyeOutlined,
  SunOutlined,
  MoonOutlined,
  AppstoreOutlined,
  LineChartOutlined,
  ClusterOutlined,
  TeamOutlined,
  CrownOutlined,
  LinkOutlined,
  BankOutlined,
  MenuFoldOutlined,
  MenuUnfoldOutlined,
} from '@ant-design/icons';
import { useAuthStore } from '../stores/authStore';
import { useThemeStore, themeTokens } from '../stores/themeStore';
import UserMenuWithQuota from './UserMenuWithQuota';
import ModernUserHeader from './ModernUserHeader';
import type { MenuProps } from 'antd';

const { Header, Sider, Content, Footer } = AntLayout;
const { Text } = Typography;

export default function Layout() {
  const navigate = useNavigate();
  const location = useLocation();
  const { user, logout } = useAuthStore();
  const { mode, toggleTheme } = useThemeStore();
  const [collapsed, setCollapsed] = useState(false);
  const [userOrgInfo, setUserOrgInfo] = useState<any>(null);
  const [loadingOrgInfo, setLoadingOrgInfo] = useState(false);

  // 获取当前主题的颜色
  const colors = themeTokens[mode];
  const isDark = mode === 'dark';

  // 根据当前路径确定选中的菜单项
  const getSelectedKey = () => {
    const path = location.pathname;

    // 系统管理子菜单 - 精确匹配子页面
    if (path === '/workspace/admin/users') return '/workspace/admin/users';
    if (path.startsWith('/workspace/admin/users/')) return '/workspace/admin/users';
    if (path === '/workspace/admin/organizations') return '/workspace/admin/organizations';
    if (path.startsWith('/workspace/admin/organizations/')) return '/workspace/admin/organizations';
    if (path === '/workspace/admin/billing') return '/workspace/admin/billing';
    if (path === '/workspace/admin/visibility') return '/workspace/admin/visibility';
    if (path === '/workspace/admin/permissions-billing') return '/workspace/admin/permissions-billing';
    if (path === '/workspace/admin/logs') return '/workspace/admin/logs';
    if (path === '/workspace/admin' || path === '/workspace/admin/') return '/workspace/admin';
    // 主账号管理和用户定价高亮各自的菜单项
    if (path === '/workspace/admin/master-accounts') return '/workspace/admin/master-accounts';
    if (path === '/workspace/admin/pricing') return '/workspace/admin/pricing';
    if (path.startsWith('/workspace/admin/master-accounts/')) return '/workspace/admin/master-accounts';
    if (path.startsWith('/workspace/admin/pricing/')) return '/workspace/admin/pricing';

    // 溶液电解质子菜单
    if (path.startsWith('/workspace/liquid-electrolyte/electrolytes')) return '/workspace/liquid-electrolyte/electrolytes';
    if (path.startsWith('/workspace/liquid-electrolyte/analysis')) return '/workspace/liquid-electrolyte/analysis';
    if (path.startsWith('/workspace/liquid-electrolyte/qc')) return '/workspace/liquid-electrolyte/qc';
    if (path.startsWith('/workspace/liquid-electrolyte/ai-discovery')) return '/workspace/liquid-electrolyte/ai-discovery';
    if (path.startsWith('/workspace/liquid-electrolyte/anion-generation')) return '/workspace/liquid-electrolyte/anion-generation';
    if (path.startsWith('/workspace/liquid-electrolyte/md') || path.startsWith('/workspace/liquid-electrolyte')) return '/workspace/liquid-electrolyte/md';

    // 旧路由兼容
    if (path.startsWith('/workspace/qc-jobs')) return '/workspace/liquid-electrolyte/qc';
    if (path.startsWith('/workspace/jobs')) return '/workspace/liquid-electrolyte/md';

    // 账户中心（包含所有账户和组织相关功能）
    if (path.startsWith('/workspace/account-center')) return '/workspace/account-center';
    if (path.startsWith('/workspace/account-dashboard')) return '/workspace/account-center';
    if (path.startsWith('/workspace/account-settings')) return '/workspace/account-center';
    if (path.startsWith('/workspace/quota-recharge')) return '/workspace/account-center';
    if (path.startsWith('/workspace/my-organizations')) return '/workspace/account-center';
    if (path.startsWith('/workspace/sub-accounts')) return '/workspace/account-center';
    // 兼容旧路由
    if (path.startsWith('/workspace/my-account')) return '/workspace/account-dashboard';
    if (path.startsWith('/workspace/account-management')) return '/workspace/account-dashboard';

    // 其他页面
    if (path.startsWith('/workspace/electrolytes')) return '/workspace/electrolytes';
    if (path.startsWith('/workspace/projects')) return '/workspace/projects';
    if (path.startsWith('/workspace/research')) return '/workspace/research';
    if (path.startsWith('/workspace/data-visibility')) return '/workspace/data-visibility';
    if (path.startsWith('/workspace/dashboard')) return '/workspace/dashboard';

    return '/workspace/dashboard';
  };

  // 获取打开的子菜单
  const getOpenKeys = () => {
    const path = location.pathname;
    const keys: string[] = [];

    if (path.startsWith('/workspace/admin')) {
      keys.push('/workspace/admin');
    }
    if (path.startsWith('/workspace/liquid-electrolyte') ||
        path.startsWith('/workspace/jobs') ||
        path.startsWith('/workspace/qc-jobs')) {
      keys.push('/workspace/liquid-electrolyte');
    }
    // 账户中心不需要子菜单展开（单一页面）
    return keys;
  };

  // 检查用户是否有权限访问某个模块
  const hasModuleAccess = (moduleName: string): boolean => {
    // 管理员可以访问所有模块
    if (user?.role === 'ADMIN') return true;

    // 如果用户没有 allowed_modules 字段，默认不允许访问任何模块（严格模式）
    // 这确保了模块权限的强制执行
    if (!user?.allowed_modules || user.allowed_modules.length === 0) return false;

    // 检查用户是否有权限访问该模块
    return user.allowed_modules.includes(moduleName);
  };

  // 侧边栏菜单项
  const menuItems: MenuProps['items'] = [
    {
      key: '/workspace/dashboard',
      icon: <DashboardOutlined />,
      label: '个人中心',
    },
    {
      key: '/workspace/projects',
      icon: <ProjectOutlined />,
      label: '项目管理',
    },
    // 溶液电解质模块（二级菜单）
    {
      key: '/workspace/liquid-electrolyte',
      icon: <AppstoreOutlined />,
      label: '溶元调配',
      children: [
        ...(hasModuleAccess('electrolytes') ? [{
          key: '/workspace/liquid-electrolyte/electrolytes',
          icon: <ExperimentOutlined />,
          label: '溶液配方管理',
        }] : []),
        ...(hasModuleAccess('md') ? [{
          key: '/workspace/liquid-electrolyte/md',
          icon: <RocketOutlined />,
          label: '溶液MD分析',
        }] : []),
        ...(hasModuleAccess('analysis') ? [{
          key: '/workspace/liquid-electrolyte/analysis',
          icon: <LineChartOutlined />,
          label: '溶鞘QC分析',
        }] : []),
        ...(hasModuleAccess('qc') ? [{
          key: '/workspace/liquid-electrolyte/qc',
          icon: <ThunderboltOutlined />,
          label: '溶元QC分析',
        }] : []),
        ...(hasModuleAccess('anion-generation') ? [{
          key: '/workspace/liquid-electrolyte/anion-generation',
          icon: <ExperimentOutlined />,
          label: '溶盐FF开发',
        }] : []),
        ...(hasModuleAccess('ai-discovery') ? [{
          key: '/workspace/liquid-electrolyte/ai-discovery',
          icon: <ClusterOutlined />,
          label: '溶元AI推荐',
        }] : []),
      ],
    },
    {
      key: '/workspace/research',
      icon: <DatabaseOutlined />,
      label: '数据管理',
    },
    // Admin menu - only show for admin users
    ...(user?.role === 'ADMIN' ? [{
      key: '/workspace/admin',
      icon: <ControlOutlined />,
      label: '系统管理',
      children: [
        {
          key: '/workspace/admin',
          label: '管理面板',
        },
        {
          key: '/workspace/admin/users',
          label: '用户管理',
        },
        {
          key: '/workspace/admin/master-accounts',
          label: '主账号管理',
        },
        {
          key: '/workspace/admin/billing',
          label: '计费管理',
        },
        {
          key: '/workspace/admin/pricing',
          label: '用户定价',
        },
        {
          key: '/workspace/admin/visibility',
          label: '数据公开',
        },
        {
          key: '/workspace/admin/logs',
          label: '审计日志',
        },
      ],
    }] : []),
    // 账户中心 - 统一的账户管理入口（包含组织管理）
    {
      key: '/workspace/account-center',
      icon: <UserOutlined />,
      label: '账户中心',
    },
  ];

  // 用户下拉菜单
  console.log('🎨 渲染用户菜单，账号信息状态:', {
    userOrgInfo,
    loadingOrgInfo,
    hasMasterAccount: !!userOrgInfo?.master_account,
    hasSubAccount: !!userOrgInfo?.sub_account
  });

  const userMenuItems: MenuProps['items'] = [
    // Phase 3: 账号类型和配额信息
    {
      key: 'account-quota',
      label: <UserMenuWithQuota compact={false} />,
      disabled: true,
    },
    {
      type: 'divider',
    },
    // 主账号信息
    ...(userOrgInfo && userOrgInfo.master_account ? [
      {
        key: 'master-account-info',
        label: (
          <div style={{ padding: '8px 0' }}>
            <div style={{
              fontSize: '12px',
              color: colors.colorTextSecondary,
              marginBottom: '8px',
              fontWeight: 500
            }}>
              主账号信息
            </div>
            <div style={{ padding: '4px 0' }}>
              <div style={{
                display: 'flex',
                alignItems: 'center',
                gap: '6px',
                marginBottom: '2px'
              }}>
                <CrownOutlined style={{ fontSize: '12px', color: '#faad14' }} />
                <span style={{
                  fontSize: '13px',
                  fontWeight: 500,
                  color: colors.colorText
                }}>
                  {userOrgInfo.master_account.organization_name}
                </span>
                <Tag color="gold" style={{ fontSize: '10px', lineHeight: '16px' }}>
                  主账号
                </Tag>
              </div>
              <div style={{
                fontSize: '11px',
                color: colors.colorTextSecondary,
                marginLeft: '18px'
              }}>
                总配额: {userOrgInfo.master_account.total_cpu_hours !== undefined && userOrgInfo.master_account.total_cpu_hours !== null ? userOrgInfo.master_account.total_cpu_hours.toFixed(1) : '0.0'}h /
                可用: {userOrgInfo.master_account.balance_cpu_hours !== undefined && userOrgInfo.master_account.balance_cpu_hours !== null ? userOrgInfo.master_account.balance_cpu_hours.toFixed(1) : '0.0'}h
              </div>
              <div style={{
                fontSize: '11px',
                color: colors.colorTextSecondary,
                marginLeft: '18px'
              }}>
                子账号: {userOrgInfo.master_account.current_sub_accounts}/{userOrgInfo.master_account.max_sub_accounts}
              </div>
            </div>
          </div>
        ),
        disabled: true,
      },
    ] : []),

    // 子账号信息
    ...(userOrgInfo && userOrgInfo.sub_account ? [
      {
        key: 'sub-account-info',
        label: (
          <div style={{ padding: '8px 0' }}>
            <div style={{
              fontSize: '12px',
              color: colors.colorTextSecondary,
              marginBottom: '8px',
              fontWeight: 500
            }}>
              子账号信息
            </div>
            <div style={{ padding: '4px 0' }}>
              <div style={{
                display: 'flex',
                alignItems: 'center',
                gap: '6px',
                marginBottom: '2px'
              }}>
                <LinkOutlined style={{ fontSize: '12px', color: colors.colorPrimary }} />
                <span style={{
                  fontSize: '13px',
                  fontWeight: 500,
                  color: colors.colorText
                }}>
                  关联到: {userOrgInfo.sub_account.master_username}
                </span>
                <Tag color="blue" style={{ fontSize: '10px', lineHeight: '16px' }}>
                  子账号
                </Tag>
              </div>
              <div style={{
                fontSize: '11px',
                color: colors.colorTextSecondary,
                marginLeft: '18px'
              }}>
                个人配额: {userOrgInfo.sub_account.personal_quota !== undefined && userOrgInfo.sub_account.personal_quota !== null ? userOrgInfo.sub_account.personal_quota.toFixed(1) : '0.0'}h /
                已用: {userOrgInfo.sub_account.personal_used !== undefined && userOrgInfo.sub_account.personal_used !== null ? userOrgInfo.sub_account.personal_used.toFixed(1) : '0.0'}h
              </div>
            </div>
          </div>
        ),
        disabled: true,
      },
    ] : []),

    // 分隔线（如果有账号信息）
    ...(userOrgInfo && (userOrgInfo.master_account || userOrgInfo.sub_account) ? [
      { type: 'divider' as const }
    ] : []),

    {
      key: 'profile',
      icon: <UserOutlined />,
      label: '个人信息',
    },
    {
      key: 'recharge',
      icon: <WalletOutlined />,
      label: '充值中心',
    },
    {
      key: 'settings',
      icon: <SettingOutlined />,
      label: '修改密码',
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
    // 直接导航到对应路由
    navigate(key);
  };

  const handleUserMenuClick = ({ key }: { key: string }) => {
    if (key === 'logout') {
      logout();
      navigate('/login');
    } else if (key === 'settings') {
      navigate('/workspace/change-password');
    } else if (key === 'profile') {
      navigate('/workspace/profile');
    } else if (key === 'recharge') {
      // Phase 3: 导航到新的配额充值页面
      navigate('/workspace/quota-recharge');
    }
  };

  return (
    <AntLayout style={{ minHeight: '100vh' }}>
      <Sider
        collapsible
        collapsed={collapsed}
        onCollapse={setCollapsed}
        theme="dark"
        width={220}
        collapsedWidth={72}
        trigger={null}
        style={{
          height: '100vh',
          position: 'fixed',
          left: 0,
          top: 0,
          bottom: 0,
          background: '#0f172a',
          borderRight: '1px solid rgba(255, 255, 255, 0.06)',
          zIndex: 200,
          overflow: 'hidden',
        }}
      >
        <div style={{
          height: '100%',
          display: 'flex',
          flexDirection: 'column',
        }}>
          {/* Logo 区域 */}
          <div
            onClick={() => navigate('/')}
            style={{
              height: 64,
              minHeight: 64,
              display: 'flex',
              alignItems: 'center',
              padding: collapsed ? '0 16px' : '0 20px',
              cursor: 'pointer',
              borderBottom: '1px solid rgba(255, 255, 255, 0.06)',
              flexShrink: 0,
            }}
          >
            <div style={{
              width: 36,
              height: 36,
              borderRadius: 10,
              background: 'linear-gradient(135deg, #3b82f6 0%, #8b5cf6 100%)',
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
              flexShrink: 0,
            }}>
              <ThunderboltOutlined style={{ color: '#fff', fontSize: 18 }} />
            </div>
            {!collapsed && (
              <span style={{
                marginLeft: 12,
                color: '#fff',
                fontSize: 18,
                fontWeight: 600,
                letterSpacing: 0.5,
              }}>
                Molyte
              </span>
            )}
          </div>

          {/* 菜单区域 - 可滚动 */}
          <div style={{
            flex: 1,
            overflowY: 'auto',
            overflowX: 'hidden',
            padding: '8px 0',
          }} className="sidebar-menu-scroll">
            <Menu
              theme="dark"
              mode="inline"
              selectedKeys={[getSelectedKey()]}
              defaultOpenKeys={collapsed ? [] : ['/workspace/liquid-electrolyte', ...getOpenKeys()]}
              items={menuItems}
              onClick={handleMenuClick}
              style={{
                background: 'transparent',
                border: 'none',
              }}
            />
          </div>

          {/* 底部区域：折叠按钮 + 版权 */}
          <div style={{
            borderTop: '1px solid rgba(255, 255, 255, 0.06)',
            flexShrink: 0,
          }}>
            {/* 折叠按钮 */}
            <div
              onClick={() => setCollapsed(!collapsed)}
              style={{
                height: 44,
                display: 'flex',
                alignItems: 'center',
                justifyContent: collapsed ? 'center' : 'flex-start',
                padding: collapsed ? 0 : '0 20px',
                cursor: 'pointer',
                color: 'rgba(255, 255, 255, 0.5)',
                transition: 'all 0.2s',
              }}
            >
              {collapsed ? (
                <MenuUnfoldOutlined style={{ fontSize: 16 }} />
              ) : (
                <>
                  <MenuFoldOutlined style={{ fontSize: 16 }} />
                  <span style={{ marginLeft: 12, fontSize: 13 }}>收起菜单</span>
                </>
              )}
            </div>

            {/* 版权信息 */}
            <div style={{
              padding: collapsed ? '12px 8px' : '12px 20px',
              color: 'rgba(255, 255, 255, 0.3)',
              fontSize: 11,
              textAlign: collapsed ? 'center' : 'left',
            }}>
              {collapsed ? '©' : '© 2025 Molyte'}
            </div>
          </div>
        </div>
      </Sider>

      <AntLayout style={{ marginLeft: collapsed ? 72 : 220, transition: 'margin-left 0.2s' }}>
        <Header
          style={{
            padding: '0 24px',
            background: colors.colorBgContainer,
            display: 'flex',
            justifyContent: 'space-between',
            alignItems: 'center',
            boxShadow: isDark ? '0 1px 4px rgba(0, 0, 0, 0.3)' : '0 1px 4px rgba(0, 21, 41, 0.08)',
            position: 'sticky',
            top: 0,
            zIndex: 100,
            height: 64,
            transition: 'background 0.3s',
          }}
        >
          {/* 左侧：页面标识 */}
          <div style={{ display: 'flex', alignItems: 'center' }}>
            <Text style={{
              fontSize: 15,
              color: colors.colorTextSecondary,
              fontWeight: 400,
            }}>
              工作台
            </Text>
          </div>

          {/* 右侧：操作区 */}
          <Space size={12}>
            {/* 主题切换按钮 */}
            <Tooltip title={isDark ? '切换到浅色模式' : '切换到深色模式'}>
              <Button
                type="text"
                icon={isDark ? <SunOutlined style={{ fontSize: 18, color: '#FFC53D' }} /> : <MoonOutlined style={{ fontSize: 18, color: '#5B8DEF' }} />}
                onClick={toggleTheme}
                style={{
                  width: 36,
                  height: 36,
                  borderRadius: 8,
                  display: 'flex',
                  alignItems: 'center',
                  justifyContent: 'center',
                  background: 'transparent',
                }}
              />
            </Tooltip>

            {/* 现代化用户头部组件 */}
            <ModernUserHeader onLogout={() => {
              logout();
              navigate('/login');
            }} />
          </Space>
        </Header>

        <Content style={{
          margin: 0,
          padding: 0,
          background: colors.colorBgLayout,
          minHeight: 'calc(100vh - 64px - 48px)', // 减去 header 和 footer 高度
          overflow: 'auto',
          transition: 'background 0.3s',
        }}>
          <Outlet />
        </Content>

        {/* Footer */}
        <Footer style={{
          textAlign: 'center',
          padding: '12px 24px',
          background: colors.colorBgContainer,
          borderTop: `1px solid ${colors.colorBorder}`,
          fontSize: 12,
          color: colors.colorTextSecondary,
        }}>
          © 2025 Molyte Platform
        </Footer>
      </AntLayout>
    </AntLayout>
  );
}

