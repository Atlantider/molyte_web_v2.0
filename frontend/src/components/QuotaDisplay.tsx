/**
 * QuotaDisplay 组件
 * 
 * 显示用户配额信息和使用情况
 */

import React from 'react';
import { Card, Row, Col, Statistic, Progress, Space, Empty, Spin, Divider, Tag } from 'antd';
import { ThunderboltOutlined, AlertOutlined } from '@ant-design/icons';
import {
  formatQuota,
  getQuotaStatusColor,
  getQuotaStatusText,
  getQuotaSourceLabel,
  QuotaData,
} from '../hooks/useQuota';

interface QuotaDisplayProps {
  quota: QuotaData | null;
  loading?: boolean;
  showSourceBreakdown?: boolean;
}

const QuotaDisplay: React.FC<QuotaDisplayProps> = ({
  quota,
  loading = false,
  showSourceBreakdown = true,
}) => {
  if (loading) {
    return (
      <Card style={{ marginBottom: '24px' }}>
        <Spin />
      </Card>
    );
  }

  if (!quota) {
    return (
      <Card style={{ marginBottom: '24px' }}>
        <Empty description="无法加载配额信息" />
      </Card>
    );
  }

  const available = quota.available_quota;
  const used = quota.account_details?.used_cpu_hours ?? 0;
  const total = available + used;
  const percentage = total > 0 ? Math.round((used / total) * 100) : 0;
  const statusColor = getQuotaStatusColor(percentage);
  const statusText = getQuotaStatusText(percentage);

  return (
    <Card
      style={{ marginBottom: '24px' }}
      title={
        <Space>
          <ThunderboltOutlined />
          <span>配额信息</span>
        </Space>
      }
    >
      {/* 配额概览 */}
      <Row gutter={16} style={{ marginBottom: '24px' }}>
        <Col xs={24} sm={12} md={8}>
          <Statistic
            title="可用配额"
            value={available}
            suffix="h"
            precision={1}
            valueStyle={{ color: '#52c41a' }}
          />
        </Col>

        <Col xs={24} sm={12} md={8}>
          <Statistic
            title="已用配额"
            value={used}
            suffix="h"
            precision={1}
          />
        </Col>

        <Col xs={24} sm={12} md={8}>
          <Statistic
            title="总配额"
            value={total}
            suffix="h"
            precision={1}
          />
        </Col>
      </Row>

      {/* 配额使用进度 */}
      <Row gutter={16} style={{ marginBottom: '24px' }}>
        <Col span={24}>
          <div style={{ marginBottom: '8px' }}>
            <Space>
              <span style={{ fontWeight: 'bold' }}>配额使用情况：</span>
              <Tag color={statusColor === '#ff4d4f' ? 'red' : statusColor === '#faad14' ? 'orange' : statusColor === '#1890ff' ? 'blue' : 'green'}>
                {statusText}
              </Tag>
              <span style={{ color: '#666' }}>
                {percentage}% ({formatQuota(used)} / {formatQuota(total)})
              </span>
            </Space>
          </div>
          <Progress
            percent={percentage}
            strokeColor={statusColor}
            status={percentage >= 90 ? 'exception' : 'active'}
          />
        </Col>
      </Row>

      {/* 配额来源分解 */}
      {showSourceBreakdown && quota.quota_sources && Object.keys(quota.quota_sources).length > 0 && (
        <>
          <Divider />
          <Row gutter={16}>
            <Col span={24}>
              <div style={{ fontWeight: 'bold', marginBottom: '12px' }}>配额来源：</div>
            </Col>
          </Row>

          <Row gutter={16}>
            {Object.entries(quota.quota_sources).map(([key, value]) => (
              value !== undefined && value !== null && (
                <Col xs={24} sm={12} md={8} key={key}>
                  <Card size="small" style={{ backgroundColor: '#fafafa' }}>
                    <Statistic
                      title={getQuotaSourceLabel(key)}
                      value={value}
                      suffix="h"
                      precision={1}
                    />
                  </Card>
                </Col>
              )
            ))}
          </Row>
        </>
      )}

      {/* 配额不足警告 */}
      {available < 1 && (
        <>
          <Divider />
          <Row gutter={16}>
            <Col span={24}>
              <div
                style={{
                  padding: '12px',
                  backgroundColor: '#fff7e6',
                  border: '1px solid #ffc069',
                  borderRadius: '4px',
                  display: 'flex',
                  alignItems: 'center',
                  gap: '8px',
                }}
              >
                <AlertOutlined style={{ color: '#faad14', fontSize: '16px' }} />
                <span style={{ color: '#ad6800' }}>
                  配额即将用尽，请及时充值以继续使用服务
                </span>
              </div>
            </Col>
          </Row>
        </>
      )}

      {/* 账号类型提示 */}
      {quota.account_type && (
        <>
          <Divider />
          <Row gutter={16}>
            <Col span={24}>
              <div style={{ color: '#666', fontSize: '12px' }}>
                <strong>账号类型：</strong> {quota.account_type}
              </div>
            </Col>
          </Row>
        </>
      )}
    </Card>
  );
};

export default QuotaDisplay;

