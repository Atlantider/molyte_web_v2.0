/**
 * 统一的账户管理页面
 * 包含用户定价和主账号管理两个 Tab
 */
import React, { useState, useEffect } from 'react';
import {
  Tabs,
  Card,
  Spin,
  message,
} from 'antd';
import {
  DollarOutlined,
  UserOutlined,
} from '@ant-design/icons';
import AdminNav from '../../components/AdminNav';
import UserPricingTab from './tabs/UserPricingTab';
import MasterAccountTab from './tabs/MasterAccountTab';

const AccountManagement: React.FC = () => {
  const [activeTab, setActiveTab] = useState('pricing');
  const [loading, setLoading] = useState(false);

  const tabItems = [
    {
      key: 'pricing',
      label: (
        <span>
          <DollarOutlined style={{ marginRight: 8 }} />
          用户定价管理
        </span>
      ),
      children: <UserPricingTab />,
    },
    {
      key: 'master-accounts',
      label: (
        <span>
          <UserOutlined style={{ marginRight: 8 }} />
          主账号管理
        </span>
      ),
      children: <MasterAccountTab />,
    },
  ];

  return (
    <div style={{ padding: '24px', background: '#f5f7fa', minHeight: '100vh' }}>
      <AdminNav />
      
      <Card
        style={{
          borderRadius: '12px',
          boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
          border: 'none',
          marginTop: '24px',
        }}
      >
        <Spin spinning={loading}>
          <Tabs
            activeKey={activeTab}
            onChange={setActiveTab}
            items={tabItems}
            size="large"
          />
        </Spin>
      </Card>
    </div>
  );
};

export default AccountManagement;

