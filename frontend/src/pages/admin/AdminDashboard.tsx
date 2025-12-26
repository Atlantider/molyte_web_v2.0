/**
 * Admin Dashboard Page
 */
import React, { useState, useEffect } from 'react';
import { Card, Row, Col, Statistic, Table, Tag, Progress, Spin, message, Typography, Space, theme, Button, Tooltip } from 'antd';
import {
  UserOutlined,
  RocketOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  ThunderboltOutlined,
  TeamOutlined,
  ControlOutlined,
  DatabaseOutlined,
  ReloadOutlined,
  ExperimentOutlined,
} from '@ant-design/icons';
import { Column, Pie } from '@ant-design/plots';
import { useNavigate, useLocation } from 'react-router-dom';
import AdminNav from '../../components/AdminNav';
import {
  getGlobalStats,
  getUserUsageStats,
  getCpuUsageRanking,
  getJobCountRanking,
  GlobalStats,
  UserUsageStatsItem,
  UserRanking,
} from '../../api/admin';
import { refreshIons } from '../../api/electrolytes';
import { useThemeStore } from '../../stores/themeStore';

const { Title, Text } = Typography;

const AdminDashboard: React.FC = () => {
  const navigate = useNavigate();
  const location = useLocation();
  const { token } = theme.useToken();
  const { mode } = useThemeStore();
  const isDark = mode === 'dark';
  const [loading, setLoading] = useState(false);
  const [stats, setStats] = useState<GlobalStats | null>(null);
  const [userStats, setUserStats] = useState<UserUsageStatsItem[]>([]);
  const [cpuRanking, setCpuRanking] = useState<UserRanking[]>([]);
  const [jobRanking, setJobRanking] = useState<UserRanking[]>([]);
  const [refreshingIons, setRefreshingIons] = useState(false);

  // 每次导航到此页面时刷新数据（使用 key 强制刷新）
  useEffect(() => {
    loadData();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [location.key, location.pathname]);

  const loadData = async () => {
    setLoading(true);
    try {
      const [globalStats, usageStats, cpuRank, jobRank] = await Promise.all([
        getGlobalStats(),
        getUserUsageStats({ limit: 10, sort_by: 'usage_percentage' }),
        getCpuUsageRanking(5),
        getJobCountRanking(5),
      ]);

      setStats(globalStats);
      setUserStats(usageStats);
      setCpuRanking(cpuRank);
      setJobRanking(jobRank);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '加载数据失败');
    } finally {
      setLoading(false);
    }
  };

  const handleRefreshIons = async () => {
    setRefreshingIons(true);
    try {
      const result = await refreshIons();
      message.success(`离子缓存刷新成功！阳离子: ${result.cations.length} 个，阴离子: ${result.anions.length} 个`);
    } catch (error: any) {
      message.error('刷新离子缓存失败');
    } finally {
      setRefreshingIons(false);
    }
  };

  if (loading || !stats) {
    return (
      <div style={{ textAlign: 'center', padding: '100px 0' }}>
        <Spin size="large" />
      </div>
    );
  }

  // Calculate usage percentage
  const cpuUsagePercentage = stats.total_cpu_hours_allocated > 0
    ? (stats.total_cpu_hours_used / stats.total_cpu_hours_allocated) * 100
    : 0;

  const storageUsagePercentage = stats.total_storage_allocated_gb > 0
    ? (stats.total_storage_used_gb / stats.total_storage_allocated_gb) * 100
    : 0;

  // Prepare job status distribution data
  const jobStatusData = [
    { type: '已完成', value: stats.completed_jobs },
    { type: '运行中', value: stats.running_jobs },
    { type: '排队中', value: stats.queued_jobs },
    { type: '失败', value: stats.failed_jobs },
  ].filter(item => item.value > 0);

  // Prepare CPU ranking chart data (格式化为1位小数避免显示过长)
  const cpuRankingData = cpuRanking.map(item => ({
    username: item.username,
    value: Math.round(item.metric_value * 10) / 10,
  }));

  return (
    <div style={{ padding: '20px 24px', background: token.colorBgLayout, minHeight: 'calc(100vh - 64px)' }}>
      {/* 页面标题 - 统一风格 */}
      <div style={{ marginBottom: 16 }}>
        <Title level={3} className="admin-dashboard-title" style={{ margin: 0, marginBottom: 4 }}>
          <ControlOutlined style={{ marginRight: 10, color: token.colorPrimary }} />
          管理面板
        </Title>
        <Text className="admin-dashboard-desc">系统资源监控与用户管理</Text>
      </div>

      {/* Admin Navigation Menu */}
      <AdminNav />

      {/* KPI Cards - 统一简洁风格 */}
      <Row gutter={16} style={{ marginBottom: 20 }}>
        {[
          { label: '总用户数', value: stats.total_users, color: '#85a5ff', icon: <TeamOutlined /> },
          { label: '总任务数', value: stats.total_jobs, color: '#b37feb', icon: <RocketOutlined /> },
          { label: '已完成任务', value: stats.completed_jobs, color: '#52c41a', icon: <CheckCircleOutlined /> },
          { label: 'CPU使用率', value: `${cpuUsagePercentage.toFixed(1)}%`, color: '#faad14', icon: <ThunderboltOutlined />, isText: true },
        ].map((item, idx) => (
          <Col xs={12} sm={6} key={idx}>
            <div style={{
              padding: '16px 20px',
              background: isDark ? 'rgba(255,255,255,0.03)' : 'white',
              borderRadius: 10,
              border: `1px solid ${token.colorBorder}`,
              borderLeft: `4px solid ${item.color}`,
              display: 'flex',
              alignItems: 'center',
              gap: 14,
            }}>
              <div style={{
                width: 42,
                height: 42,
                borderRadius: 10,
                background: `${item.color}15`,
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                fontSize: 20,
                color: item.color,
              }}>
                {item.icon}
              </div>
              <div>
                <Text type="secondary" style={{ fontSize: 13, display: 'block' }}>{item.label}</Text>
                <Text strong style={{ fontSize: item.isText ? 18 : 22, color: item.color }}>{item.value}</Text>
              </div>
            </div>
          </Col>
        ))}
      </Row>

      {/* Charts Row */}
      <Row gutter={[16, 16]} style={{ marginBottom: 24 }}>
        {/* CPU Usage Ranking */}
        <Col xs={24} lg={12}>
          <Card
            title={
              <Space>
                <ThunderboltOutlined style={{ color: token.colorPrimary }} />
                <span>CPU 核时使用排行 Top 5</span>
              </Space>
            }
            bordered={false}
            style={{
              borderRadius: 12,
              border: 'none',
              boxShadow: isDark ? '0 2px 8px rgba(0,0,0,0.3)' : '0 2px 8px rgba(0,0,0,0.06)',
              background: token.colorBgContainer,
            }}
          >
            {cpuRankingData.length > 0 ? (
              <Column
                data={cpuRankingData}
                xField="username"
                yField="value"
                theme={isDark ? 'dark' : undefined}
                label={{
                  position: 'top',
                  style: {
                    fill: token.colorText,
                    fontSize: 12,
                  },
                  content: (originData: any) => {
                    const val = originData?.value;
                    return val !== undefined ? `${val.toFixed(1)}h` : '';
                  },
                }}
                tooltip={{
                  formatter: (datum: any) => ({
                    name: 'CPU核时',
                    value: `${datum.value?.toFixed(1)} h`,
                  }),
                }}
                color="l(270) 0:#1677ff 1:#13c2c2"
                columnStyle={{
                  radius: [8, 8, 0, 0],
                }}
                xAxis={{
                  label: {
                    autoRotate: false,
                  },
                }}
                yAxis={{
                  title: {
                    text: 'CPU 核时 (h)',
                  },
                }}
                height={300}
              />
            ) : (
              <div style={{ textAlign: 'center', padding: '50px 0', color: token.colorTextSecondary }}>
                暂无数据
              </div>
            )}
          </Card>
        </Col>

        {/* Job Status Distribution */}
        <Col xs={24} lg={12}>
          <Card
            title={
              <Space>
                <DatabaseOutlined style={{ color: token.colorPrimary }} />
                <span>任务状态分布</span>
              </Space>
            }
            bordered={false}
            style={{
              borderRadius: 12,
              border: 'none',
              boxShadow: isDark ? '0 2px 8px rgba(0,0,0,0.3)' : '0 2px 8px rgba(0,0,0,0.06)',
              background: token.colorBgContainer,
            }}
          >
            {jobStatusData.length > 0 ? (
              <Pie
                data={jobStatusData}
                angleField="value"
                colorField="type"
                radius={0.8}
                innerRadius={0.6}
                theme={isDark ? 'dark' : undefined}
                label={{
                  type: 'inner',
                  offset: '-30%',
                  content: (data: any) => data?.value ?? '',
                  style: {
                    fontSize: 14,
                    textAlign: 'center',
                  },
                }}
                legend={{
                  position: 'bottom',
                }}
                statistic={{
                  title: {
                    content: '总任务',
                  },
                  content: {
                    content: stats.total_jobs.toString(),
                  },
                }}
                height={300}
              />
            ) : (
              <div style={{ textAlign: 'center', padding: '50px 0', color: token.colorTextSecondary }}>
                暂无数据
              </div>
            )}
          </Card>
        </Col>
      </Row>

      {/* User Usage Stats Table */}
      <Card
        title={
          <Space>
            <TeamOutlined style={{ color: token.colorPrimary }} />
            用户使用情况概览
          </Space>
        }
        bordered={false}
        extra={
          <Space>
            <Button
              icon={<ReloadOutlined />}
              onClick={loadData}
              loading={loading}
              style={{ borderRadius: 8 }}
            >
              刷新
            </Button>
            <Button
              type="primary"
              onClick={() => navigate('/workspace/admin/users')}
              style={{ borderRadius: 8 }}
            >
              查看全部用户 →
            </Button>
          </Space>
        }
        style={{
          borderRadius: '12px',
          boxShadow: isDark ? '0 2px 8px rgba(0,0,0,0.3)' : '0 10px 30px rgba(15, 100, 255, 0.08)',
          background: token.colorBgContainer,
        }}
      >
        <Table
          dataSource={userStats}
          rowKey="user_id"
          pagination={false}
          columns={[
            {
              title: '用户名',
              dataIndex: 'username',
              key: 'username',
              render: (text: string, record: UserUsageStatsItem) => (
                <a onClick={() => navigate(`/workspace/admin/users/${record.user_id}`)}>
                  {text}
                </a>
              ),
            },
            {
              title: '角色',
              dataIndex: 'role',
              key: 'role',
              render: (role: string) => {
                const roleMap: any = {
                  ADMIN: { color: 'red', text: '管理员' },
                  PREMIUM: { color: 'gold', text: '高级用户' },
                  USER: { color: 'blue', text: '普通用户' },
                  GUEST: { color: 'default', text: '访客' },
                };
                const config = roleMap[role] || { color: 'default', text: role };
                return <Tag color={config.color}>{config.text}</Tag>;
              },
            },
            {
              title: 'CPU 核时使用',
              key: 'cpu_usage',
              render: (_: any, record: UserUsageStatsItem) => {
                const percent = Math.round(record.usage_percentage * 100) / 100;
                return (
                  <div>
                    <div style={{ marginBottom: '4px' }}>
                      {record.used_cpu_hours.toFixed(1)} / {record.total_cpu_hours.toFixed(1)} h
                    </div>
                    <Progress
                      percent={percent}
                      size="small"
                      format={(p) => `${p?.toFixed(1)}%`}
                      strokeColor={
                        percent > 80
                          ? '#f5222d'
                          : percent > 50
                            ? '#fa8c16'
                            : '#52c41a'
                      }
                    />
                  </div>
                );
              },
            },
            {
              title: '任务统计',
              key: 'jobs',
              render: (_: any, record: UserUsageStatsItem) => (
                <div>
                  总计: {record.total_jobs} | 运行中: {record.running_jobs} |
                  完成: {record.completed_jobs} | 失败: {record.failed_jobs}
                </div>
              ),
            },
            {
              title: '最后任务时间',
              dataIndex: 'last_job_at',
              key: 'last_job_at',
              render: (time: string | null) =>
                time ? new Date(time).toLocaleString('zh-CN') : '-',
            },
          ]}
        />
      </Card>
    </div>
  );
};

export default AdminDashboard;
