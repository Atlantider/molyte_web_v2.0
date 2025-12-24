/**
 * AccountTypeCard 组件
 * 
 * 显示用户账号类型和基本信息
 */

import React from 'react';
import { Card, Row, Col, Statistic, Tag, Space, Empty, Spin } from 'antd';
import {
  UserOutlined,
  TeamOutlined,
  CrownOutlined,
  BankOutlined,
} from '@ant-design/icons';
import { AccountInfo, getAccountTypeLabel, getAccountTypeColor } from '../hooks/useAccountType';

interface AccountTypeCardProps {
  accountInfo: AccountInfo | null;
  loading?: boolean;
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

/**
 * 获取账号类型描述
 */
const getAccountTypeDescription = (accountType: string): string => {
  const descriptions: Record<string, string> = {
    PERSONAL: '个人账号，拥有个人配额',
    MASTER_ACCOUNT: '主账号，可管理子账号和配额分配',
    SUB_ACCOUNT: '子账号，可使用个人配额和主账号配额',
  };
  return descriptions[accountType] || '';
};

const AccountTypeCard: React.FC<AccountTypeCardProps> = ({ accountInfo, loading = false }) => {
  if (loading) {
    return (
      <Card style={{ marginBottom: '24px' }}>
        <Spin />
      </Card>
    );
  }

  if (!accountInfo) {
    return (
      <Card style={{ marginBottom: '24px' }}>
        <Empty description="无法加载账号信息" />
      </Card>
    );
  }

  const { account_type, account_info, permissions } = accountInfo;
  const label = getAccountTypeLabel(account_type);
  const color = getAccountTypeColor(account_type);
  const icon = getAccountTypeIcon(account_type);
  const description = getAccountTypeDescription(account_type);

  return (
    <Card
      style={{ marginBottom: '24px' }}
      title={
        <Space>
          {icon}
          <span>账号信息</span>
        </Space>
      }
    >
      {/* 账号类型 */}
      <Row gutter={16} style={{ marginBottom: '24px' }}>
        <Col span={24}>
          <div style={{ display: 'flex', alignItems: 'center', gap: '12px' }}>
            <span style={{ fontWeight: 'bold' }}>账号类型：</span>
            <Tag color={color}>{label}</Tag>
            <span style={{ color: '#666', fontSize: '12px' }}>{description}</span>
          </div>
        </Col>
      </Row>

      {/* 配额统计 */}
      <Row gutter={16} style={{ marginBottom: '24px' }}>
        {account_info.total_cpu_hours !== undefined && (
          <Col xs={24} sm={12} md={8}>
            <Statistic
              title="总核时"
              value={account_info.total_cpu_hours}
              suffix="h"
              precision={1}
            />
          </Col>
        )}

        {account_info.balance_cpu_hours !== undefined && (
          <Col xs={24} sm={12} md={8}>
            <Statistic
              title="可用核时"
              value={account_info.balance_cpu_hours}
              suffix="h"
              precision={1}
              valueStyle={{ color: '#52c41a' }}
            />
          </Col>
        )}

        {account_info.used_cpu_hours !== undefined && (
          <Col xs={24} sm={12} md={8}>
            <Statistic
              title="已用核时"
              value={account_info.used_cpu_hours}
              suffix="h"
              precision={1}
            />
          </Col>
        )}
      </Row>

      {/* 主账号特定信息 */}
      {account_type === 'MASTER_ACCOUNT' && (
        <Row gutter={16} style={{ marginBottom: '24px' }}>
          {account_info.max_sub_accounts !== undefined && (
            <Col xs={24} sm={12} md={8}>
              <Statistic
                title="子账号"
                value={`${account_info.current_sub_accounts || 0}/${account_info.max_sub_accounts}`}
              />
            </Col>
          )}

          {account_info.frozen_cpu_hours !== undefined && (
            <Col xs={24} sm={12} md={8}>
              <Statistic
                title="冻结核时"
                value={account_info.frozen_cpu_hours}
                suffix="h"
                precision={1}
              />
            </Col>
          )}
        </Row>
      )}

      {/* 权限信息 */}
      <Row gutter={16}>
        <Col span={24}>
          <div style={{ fontWeight: 'bold', marginBottom: '8px' }}>权限：</div>
          <Space wrap>
            {permissions.can_create_jobs && <Tag color="blue">创建任务</Tag>}
            {permissions.can_manage_master_account && <Tag color="purple">管理主账号</Tag>}
            {permissions.can_manage_sub_accounts && <Tag color="orange">管理子账号</Tag>}
            {permissions.can_access_admin && <Tag color="red">管理员</Tag>}
          </Space>
        </Col>
      </Row>
    </Card>
  );
};

export default AccountTypeCard;

