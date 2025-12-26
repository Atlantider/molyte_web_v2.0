/**
 * User Daily Statistics Page
 */
import React, { useState, useEffect } from 'react';
import { Card, Row, Col, Table, Statistic, Spin, message, Typography, Space, theme, Button, Select, Tooltip } from 'antd';
import {
  ReloadOutlined,
  ThunderboltOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  RocketOutlined,
  DatabaseOutlined,
} from '@ant-design/icons';
import { Column } from '@ant-design/plots';
import dayjs from 'dayjs';
import { getMyDailyStats, DailyStats } from '../api/userStats';
import { useAuthStore } from '../stores/authStore';
import { useThemeStore } from '../stores/themeStore';

const { Title, Text } = Typography;

const UserStats: React.FC = () => {
  const { token } = theme.useToken();
  const { mode } = useThemeStore();
  const isDark = mode === 'dark';
  const { user } = useAuthStore();
  const [loading, setLoading] = useState(false);
  const [stats, setStats] = useState<DailyStats[]>([]);
  const [days, setDays] = useState(7);

  useEffect(() => {
    loadStats();
  }, [days]);

  const loadStats = async () => {
    setLoading(true);
    try {
      const data = await getMyDailyStats(days);
      setStats(data);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '加载统计数据失败');
    } finally {
      setLoading(false);
    }
  };

  // Calculate totals
  const totals = {
    jobs_submitted: stats.reduce((sum, s) => sum + s.jobs_submitted, 0),
    jobs_completed: stats.reduce((sum, s) => sum + s.jobs_completed, 0),
    jobs_failed: stats.reduce((sum, s) => sum + s.jobs_failed, 0),
    cpu_hours_used: stats.reduce((sum, s) => sum + s.cpu_hours_used, 0),
    cluster_analysis_cpu_hours: stats.reduce((sum, s) => sum + s.cluster_analysis_cpu_hours, 0),
    cluster_analysis_task_count: stats.reduce((sum, s) => sum + s.cluster_analysis_task_count, 0),
  };

  // Prepare chart data
  const chartData = stats.map(s => ({
    date: dayjs(s.date).format('MM-DD'),
    'CPU核时': s.cpu_hours_used,
    '后处理核时': s.cluster_analysis_cpu_hours,
  }));

  const columns = [
    {
      title: '日期',
      dataIndex: 'date',
      key: 'date',
      render: (text: string) => dayjs(text).format('YYYY-MM-DD'),
    },
    {
      title: '提交任务',
      dataIndex: 'jobs_submitted',
      key: 'jobs_submitted',
      align: 'center' as const,
    },
    {
      title: '完成任务',
      dataIndex: 'jobs_completed',
      key: 'jobs_completed',
      align: 'center' as const,
    },
    {
      title: '失败任务',
      dataIndex: 'jobs_failed',
      key: 'jobs_failed',
      align: 'center' as const,
    },
    {
      title: 'CPU核时',
      dataIndex: 'cpu_hours_used',
      key: 'cpu_hours_used',
      align: 'right' as const,
      render: (value: number) => value.toFixed(2),
    },
    {
      title: '后处理核时',
      dataIndex: 'cluster_analysis_cpu_hours',
      key: 'cluster_analysis_cpu_hours',
      align: 'right' as const,
      render: (value: number) => value.toFixed(2),
    },
    {
      title: '后处理任务数',
      dataIndex: 'cluster_analysis_task_count',
      key: 'cluster_analysis_task_count',
      align: 'center' as const,
    },
    {
      title: '最大并发',
      dataIndex: 'max_concurrent_jobs',
      key: 'max_concurrent_jobs',
      align: 'center' as const,
    },
  ];

  return (
    <div style={{
      background: token.colorBgLayout,
      minHeight: 'calc(100vh - 64px)',
      padding: '24px',
    }}>
      <div style={{ maxWidth: 1600, margin: '0 auto' }}>
        {/* Header */}
        <div style={{ marginBottom: 24 }}>
          <Title level={2} style={{ margin: 0, marginBottom: 8 }}>
            使用统计
          </Title>
          <Text type="secondary">
            查看您的任务和资源使用情况
          </Text>
        </div>

        <Spin spinning={loading}>
          {/* KPI Cards */}
          <Row gutter={[16, 16]} style={{ marginBottom: 24 }}>
            <Col xs={24} sm={12} lg={6}>
              <Card
                bordered={false}
                style={{
                  borderRadius: 12,
                  background: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)',
                  boxShadow: '0 4px 12px rgba(102, 126, 234, 0.3)',
                }}
              >
                <Statistic
                  title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>提交任务</span>}
                  value={totals.jobs_submitted}
                  valueStyle={{ color: '#fff', fontSize: 28 }}
                  prefix={<RocketOutlined />}
                />
              </Card>
            </Col>
            <Col xs={24} sm={12} lg={6}>
              <Card
                bordered={false}
                style={{
                  borderRadius: 12,
                  background: 'linear-gradient(135deg, #11998e 0%, #38ef7d 100%)',
                  boxShadow: '0 4px 12px rgba(17, 153, 142, 0.3)',
                }}
              >
                <Statistic
                  title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>完成任务</span>}
                  value={totals.jobs_completed}
                  valueStyle={{ color: '#fff', fontSize: 28 }}
                  prefix={<CheckCircleOutlined />}
                />
              </Card>
            </Col>
            <Col xs={24} sm={12} lg={6}>
              <Card
                bordered={false}
                style={{
                  borderRadius: 12,
                  background: 'linear-gradient(135deg, #fa709a 0%, #fee140 100%)',
                  boxShadow: '0 4px 12px rgba(250, 112, 154, 0.3)',
                }}
              >
                <Statistic
                  title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>CPU核时</span>}
                  value={totals.cpu_hours_used.toFixed(2)}
                  valueStyle={{ color: '#fff', fontSize: 28 }}
                  prefix={<ThunderboltOutlined />}
                  suffix={<span style={{ fontSize: 16, color: 'rgba(255,255,255,0.7)' }}>h</span>}
                />
              </Card>
            </Col>
            <Col xs={24} sm={12} lg={6}>
              <Card
                bordered={false}
                style={{
                  borderRadius: 12,
                  background: 'linear-gradient(135deg, #f093fb 0%, #f5576c 100%)',
                  boxShadow: '0 4px 12px rgba(245, 87, 108, 0.3)',
                }}
              >
                <Statistic
                  title={<span style={{ color: 'rgba(255,255,255,0.85)' }}>后处理任务</span>}
                  value={totals.cluster_analysis_task_count}
                  valueStyle={{ color: '#fff', fontSize: 28 }}
                  prefix={<DatabaseOutlined />}
                />
              </Card>
            </Col>
          </Row>

          {/* Chart */}
          <Card
            title="CPU核时趋势"
            bordered={false}
            extra={
              <Space>
                <Select
                  value={days}
                  onChange={setDays}
                  style={{ width: 120 }}
                  options={[
                    { label: '最近7天', value: 7 },
                    { label: '最近14天', value: 14 },
                    { label: '最近30天', value: 30 },
                  ]}
                />
                <Button
                  icon={<ReloadOutlined />}
                  onClick={loadStats}
                  loading={loading}
                >
                  刷新
                </Button>
              </Space>
            }
            style={{
              borderRadius: 12,
              boxShadow: isDark ? '0 2px 8px rgba(0,0,0,0.3)' : '0 10px 30px rgba(15, 100, 255, 0.08)',
              marginBottom: 24,
            }}
          >
            {chartData.length > 0 ? (
              <Column
                data={chartData}
                xField="date"
                yField="value"
                seriesField="type"
                isGroup={true}
                columnStyle={{
                  radius: [4, 4, 0, 0],
                }}
                color={[token.colorPrimary, token.colorSuccess]}
                legend={{
                  position: 'top-right',
                }}
                yAxis={{
                  label: {
                    formatter: (v: string) => `${v}h`,
                  },
                  title: {
                    text: 'CPU核时 (h)',
                  },
                }}
                theme={isDark ? 'dark' : 'light'}
                smooth={true}
                height={300}
              />
            ) : null}
          </Card>

          {/* Table */}
          <Card
            title="每日统计详情"
            bordered={false}
            style={{
              borderRadius: 12,
              boxShadow: isDark ? '0 2px 8px rgba(0,0,0,0.3)' : '0 10px 30px rgba(15, 100, 255, 0.08)',
            }}
          >
            <Table
              columns={columns}
              dataSource={stats}
              rowKey="date"
              pagination={{ pageSize: 10 }}
              scroll={{ x: 1200 }}
            />
          </Card>
        </Spin>
      </div>
    </div>
  );
};

export default UserStats;

