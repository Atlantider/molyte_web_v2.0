/**
 * 充值套餐卡片组件
 */
import React from 'react';
import { Card, Button, Tag, Typography, Space } from 'antd';
import { CheckOutlined } from '@ant-design/icons';
import { RechargePackage } from '../api/recharge-packages';

const { Title, Text } = Typography;

interface RechargePackageCardProps {
  package: RechargePackage;
  onSelect: (pkg: RechargePackage) => void;
  loading?: boolean;
  pricePerHour?: number;
}

const RechargePackageCard: React.FC<RechargePackageCardProps> = ({
  package: pkg,
  onSelect,
  loading = false,
  pricePerHour = 0.1
}) => {
  const cpuHoursPerYuan = 1 / pricePerHour;

  return (
    <Card
      hoverable
      style={{
        borderRadius: 12,
        border: 'none',
        boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
        height: '100%',
        display: 'flex',
        flexDirection: 'column',
        position: 'relative',
        overflow: 'hidden'
      }}
      bodyStyle={{ padding: '24px', flex: 1, display: 'flex', flexDirection: 'column' }}
    >
      {/* 背景色条 */}
      <div
        style={{
          position: 'absolute',
          top: 0,
          left: 0,
          right: 0,
          height: 4,
          backgroundColor: pkg.color
        }}
      />

      {/* 徽章 */}
      {pkg.badge && (
        <div style={{ marginBottom: 12 }}>
          <Tag color={pkg.color} style={{ borderRadius: 4 }}>
            {pkg.badge}
          </Tag>
        </div>
      )}

      {/* 套餐名称 */}
      <Title level={4} style={{ margin: '0 0 8px 0' }}>
        {pkg.icon && <span style={{ marginRight: 8 }}>{pkg.icon}</span>}
        {pkg.name}
      </Title>

      {/* 描述 */}
      {pkg.description && (
        <Text type="secondary" style={{ fontSize: 12, marginBottom: 16, display: 'block' }}>
          {pkg.description}
        </Text>
      )}

      {/* 价格 */}
      <div style={{ margin: '16px 0', flex: 1 }}>
        <div style={{ marginBottom: 8 }}>
          <Text strong style={{ fontSize: 24, color: pkg.color }}>
            ¥{pkg.price.toFixed(0)}
          </Text>
        </div>
        <div>
          <Text type="secondary">
            获得 <Text strong style={{ color: '#52c41a' }}>{pkg.cpu_hours.toFixed(0)}</Text> 核时
          </Text>
        </div>
        <div style={{ marginTop: 4 }}>
          <Text type="secondary" style={{ fontSize: 12 }}>
            相当于 ¥{(pkg.price / pkg.cpu_hours).toFixed(4)}/核时
          </Text>
        </div>
      </div>

      {/* 按钮 */}
      <Button
        type="primary"
        block
        size="large"
        loading={loading}
        onClick={() => onSelect(pkg)}
        style={{
          borderRadius: 8,
          backgroundColor: pkg.color,
          borderColor: pkg.color,
          marginTop: 'auto'
        }}
      >
        立即充值
      </Button>
    </Card>
  );
};

export default RechargePackageCard;

