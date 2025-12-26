/**
 * 仪表盘页面 - 重新设计版
 */
import { useEffect, useState } from 'react';
import {
  Card, Row, Col, Button, Space, Typography, Spin,
  Progress, Table, Tag, Badge, Segmented, Empty, theme
} from 'antd';
import {
  ProjectOutlined,
  ExperimentOutlined,
  RocketOutlined,
  ClockCircleOutlined,
  ThunderboltOutlined,
  SyncOutlined,
  BarChartOutlined,
  PieChartOutlined,
  LineChartOutlined,
} from '@ant-design/icons';
import { Column, Pie } from '@ant-design/plots';
import { useNavigate } from 'react-router-dom';
import { useAuthStore } from '../stores/authStore';
import { useThemeStore } from '../stores/themeStore';
import { getProjects } from '../api/projects';
import { getElectrolytes } from '../api/electrolytes';
import { getMDJobs } from '../api/jobs';
import { JobStatus } from '../types';
import type { MDJob } from '../types';
import StatusTag from '../components/StatusTag';
import dayjs from 'dayjs';

const { Text, Title } = Typography;

// 计算 CPU 核时
// 使用从 Slurm 获取的实际核时（CPUTimeRAW），而不是时间差
// 包括 MD 核时和 RESP 核时的总和
const calculateCPUHours = (job: MDJob): number => {
  // 只有已完成或失败的任务才有实际核时
  if (job.status !== 'COMPLETED' && job.status !== 'FAILED') {
    return 0;
  }

  // 计算 MD 核时
  const mdHours = (job.actual_cpu_hours && job.actual_cpu_hours > 0) ? job.actual_cpu_hours : 0;

  // 计算 RESP 核时
  const respHours = (job.resp_cpu_hours && job.resp_cpu_hours > 0) ? job.resp_cpu_hours : 0;

  // 返回总核时 = MD 核时 + RESP 核时
  return mdHours + respHours;
};

export default function Dashboard() {
  const { user } = useAuthStore();
  const { mode } = useThemeStore();
  const { token } = theme.useToken();
  const isDark = mode === 'dark';
  const navigate = useNavigate();
  const [loading, setLoading] = useState(false);

  //动态 KPI 卡片样式
  const kpiCardStyle: React.CSSProperties = {
    borderRadius: 16,
    background: isDark
      ? 'rgba(30, 30, 30, 0.8)'
      : 'rgba(255, 255, 255, 0.9)',
    backdropFilter: 'blur(12px)',
    WebkitBackdropFilter: 'blur(12px)',
    border: isDark ? '1px solid rgba(255, 255, 255, 0.12)' : '1px solid rgba(0, 0, 0, 0.06)',
    boxShadow: isDark
      ? '0 4px 12px rgba(0, 0, 0, 0.3)'
      : '0 4px 12px rgba(0, 0, 0, 0.08)',
    height: '100%',
    transition: 'all 0.3s cubic-bezier(0.4, 0, 0.2, 1)',
  };
  const [timeRange, setTimeRange] = useState<'today' | '7days' | '30days'>('7days');
  const [chartType, setChartType] = useState<'count' | 'hours'>('count');
  const [taskFilter, setTaskFilter] = useState<'all' | 'running' | 'completed' | 'failed'>('all');
  const [stats, setStats] = useState({
    projectCount: 0,
    electrolyteCount: 0,
    runningJobCount: 0,
    queuedJobCount: 0,
    failedJobCount: 0,
    completedJobCount: 0,
    todayJobCount: 0,
    totalCPUHours: 0,
    usedCPUHours: 0,
    weekNewElectrolytes: 0,
  });
  const [allJobs, setAllJobs] = useState<MDJob[]>([]);
  const [runningJobs, setRunningJobs] = useState<MDJob[]>([]);
  const [trendData, setTrendData] = useState<any[]>([]);
  const [statusDistribution, setStatusDistribution] = useState<any[]>([]);
  const [projectDistribution, setProjectDistribution] = useState<any[]>([]);

  // 加载统计数据
  const loadStats = async () => {
    setLoading(true);
    try {
      const [projects, electrolytes, jobs] = await Promise.all([
        getProjects(),
        getElectrolytes(),
        getMDJobs(),
      ]);

      console.log('Loaded jobs:', jobs);

      // 计算各种统计数据
      const runningCount = jobs.filter(j => j.status === JobStatus.RUNNING).length;
      const queuedCount = jobs.filter(j => j.status === JobStatus.QUEUED).length;
      const failedCount = jobs.filter(j => j.status === JobStatus.FAILED).length;
      const completedCount = jobs.filter(j => j.status === JobStatus.COMPLETED).length;

      // 今日提交的任务数
      const today = dayjs().startOf('day');
      const todayCount = jobs.filter(j => dayjs(j.created_at).isAfter(today)).length;

      // 本周新增配方数
      const weekAgo = dayjs().subtract(7, 'day');
      const weekNewCount = electrolytes.filter(e =>
        dayjs(e.created_at).isAfter(weekAgo)
      ).length;

      // 使用API返回的准确核时数据，fallback到前端计算
      const usedHours = (user?.used_cpu_hours !== undefined && user?.used_cpu_hours !== null)
        ? user.used_cpu_hours
        : jobs
          .filter(j => j.status === JobStatus.COMPLETED || j.status === JobStatus.FAILED)
          .reduce((sum, job) => sum + calculateCPUHours(job), 0);

      // 从用户配置中获取总机时配额
      const totalHours = user?.total_cpu_hours || 200;

      setStats({
        projectCount: projects.length,
        electrolyteCount: electrolytes.length,
        runningJobCount: runningCount,
        queuedJobCount: queuedCount,
        failedJobCount: failedCount,
        completedJobCount: completedCount,
        todayJobCount: todayCount,
        totalCPUHours: totalHours,
        usedCPUHours: usedHours,
        weekNewElectrolytes: weekNewCount,
      });

      // 设置所有任务和运行中的任务
      const sortedJobs = jobs.sort((a, b) =>
        new Date(b.created_at).getTime() - new Date(a.created_at).getTime()
      );
      setAllJobs(sortedJobs);

      setRunningJobs(jobs.filter(j =>
        j.status === JobStatus.RUNNING || j.status === JobStatus.QUEUED
      ).slice(0, 5));

      // 生成趋势数据（最近7天）
      generateTrendData(jobs, chartType);

      // 生成状态分布数据
      generateStatusDistribution(jobs);

      // 生成项目分布数据
      generateProjectDistribution(jobs, projects);

    } catch (error) {
      console.error('加载统计数据失败:', error);
    } finally {
      setLoading(false);
    }
  };

  // 生成趋势数据
  const generateTrendData = (jobs: MDJob[], type: 'count' | 'hours' = 'count') => {
    const days = 7;
    const data: any[] = [];

    for (let i = days - 1; i >= 0; i--) {
      const date = dayjs().subtract(i, 'day');
      const dateStr = date.format('MM-DD');

      const dayJobs = jobs.filter(j =>
        dayjs(j.created_at).format('YYYY-MM-DD') === date.format('YYYY-MM-DD')
      );

      if (type === 'count') {
        const completed = dayJobs.filter(j => j.status === JobStatus.COMPLETED).length;
        const failed = dayJobs.filter(j => j.status === JobStatus.FAILED).length;
        const running = dayJobs.filter(j => j.status === JobStatus.RUNNING).length;

        data.push({ date: dateStr, type: '完成', value: completed });
        data.push({ date: dateStr, type: '失败', value: failed });
        data.push({ date: dateStr, type: '运行中', value: running });
      } else {
        // 按机时
        const completedHours = dayJobs
          .filter(j => j.status === JobStatus.COMPLETED)
          .reduce((sum, j) => sum + calculateCPUHours(j), 0);
        const failedHours = dayJobs
          .filter(j => j.status === JobStatus.FAILED)
          .reduce((sum, j) => sum + calculateCPUHours(j), 0);
        const runningHours = dayJobs
          .filter(j => j.status === JobStatus.RUNNING)
          .reduce((sum, j) => sum + calculateCPUHours(j), 0);

        data.push({ date: dateStr, type: '完成', value: Number(completedHours.toFixed(2)) });
        data.push({ date: dateStr, type: '失败', value: Number(failedHours.toFixed(2)) });
        data.push({ date: dateStr, type: '运行中', value: Number(runningHours.toFixed(2)) });
      }
    }

    console.log('Trend Data:', data);
    setTrendData(data);
  };

  // 生成状态分布数据
  const generateStatusDistribution = (jobs: MDJob[]) => {
    const distribution = [
      { type: '已完成', value: jobs.filter(j => j.status === JobStatus.COMPLETED).length },
      { type: '运行中', value: jobs.filter(j => j.status === JobStatus.RUNNING).length },
      { type: '排队中', value: jobs.filter(j => j.status === JobStatus.QUEUED).length },
      { type: '失败', value: jobs.filter(j => j.status === JobStatus.FAILED).length },
      { type: '已取消', value: jobs.filter(j => j.status === JobStatus.CANCELLED).length },
    ].filter(item => item.value > 0);

    console.log('Status Distribution:', distribution);
    setStatusDistribution(distribution);
  };

  // 生成项目分布数据（按任务ID显示核时占用Top 5）
  const generateProjectDistribution = (jobs: MDJob[], projects: any[]) => {
    // 只统计已完成或运行中的任务，并按核时排序
    const jobsWithHours = jobs
      .filter(j => j.status === JobStatus.COMPLETED || j.status === JobStatus.RUNNING)
      .map(job => ({
        id: job.id,
        hours: calculateCPUHours(job),
      }))
      .filter(item => item.hours > 0)
      .sort((a, b) => b.hours - a.hours)
      .slice(0, 5); // 只显示前5个

    const distribution = jobsWithHours.map(item => ({
      type: `任务 #${item.id}`,
      value: Number(item.hours.toFixed(2)),
    }));

    console.log('Project Distribution:', distribution);
    setProjectDistribution(distribution);
  };

  useEffect(() => {
    loadStats();
  }, []);

  // 当图表类型切换时，重新生成趋势数据
  useEffect(() => {
    if (allJobs.length > 0) {
      generateTrendData(allJobs, chartType);
    }
  }, [chartType]);

  return (
    <div style={{
      background: token.colorBgLayout,
      minHeight: 'calc(100vh - 64px)',
      transition: 'background 0.3s',
    }}>
      {/* 顶部渐变背景条 */}
      <div style={{
        background: isDark
          ? 'linear-gradient(135deg, rgba(107, 154, 255, 0.08) 0%, rgba(124, 110, 175, 0.08) 100%)'
          : 'linear-gradient(135deg, rgba(91, 141, 239, 0.05) 0%, rgba(124, 110, 175, 0.05) 100%)',
        padding: '32px 24px 24px',
        marginBottom: 24,
      }}>
        <div style={{ maxWidth: 1600, margin: '0 auto' }}>
          <Title level={2} style={{ margin: 0, marginBottom: 8, color: token.colorPrimary }}>
            欢迎回来，{user?.username}！
          </Title>
          <Space size={16}>
            <Text type="secondary">
              分子动力学云平台 · 您的计算工作台
            </Text>
            <Text type="secondary">|</Text>
            <Text type="secondary">
              今日提交 <Text strong style={{ color: token.colorPrimary }}>{stats.todayJobCount}</Text> 个任务
            </Text>
            <Text type="secondary">|</Text>
            <Text type="secondary">
              机时使用率 <Text strong style={{ color: stats.usedCPUHours / stats.totalCPUHours > 0.8 ? token.colorError : token.colorSuccess }}>
                {((stats.usedCPUHours / stats.totalCPUHours) * 100).toFixed(1)}%
              </Text>
            </Text>
          </Space>
        </div>
      </div>

      <div style={{ padding: '0 24px 24px', maxWidth: 1600, margin: '0 auto' }}>
        <Spin spinning={loading}>
          {/* 第一层：全局概览 KPI 卡片 */}
          <Row gutter={[16, 16]} style={{ marginBottom: 24 }}>
            {/* 运行中任务 */}
            <Col xs={24} sm={12} lg={6}>
              <Card
                style={kpiCardStyle}
                hoverable
                onClick={() => navigate('/workspace/liquid-electrolyte/md')}
                onMouseEnter={(e) => {
                  e.currentTarget.style.transform = 'translateY(-4px)';
                  e.currentTarget.style.boxShadow = isDark
                    ? '0 12px 24px rgba(0, 0, 0, 0.4), 0 0 40px rgba(91, 141, 239, 0.2)'
                    : '0 12px 24px rgba(0, 0, 0, 0.12), 0 0 40px rgba(91, 141, 239, 0.25)';
                }}
                onMouseLeave={(e) => {
                  e.currentTarget.style.transform = 'translateY(0)';
                  e.currentTarget.style.boxShadow = kpiCardStyle.boxShadow as string;
                }}
              >
                <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
                  <div style={{ flex: 1 }}>
                    <Text type="secondary" style={{ fontSize: 12, display: 'block', marginBottom: 8, fontWeight: 600 }}>
                      RUNNING JOBS
                    </Text>
                    <div style={{ fontSize: 36, fontWeight: 800, color: token.colorPrimary, marginBottom: 8, lineHeight: 1 }}>
                      {stats.runningJobCount}
                    </div>
                    <Space size={4} style={{ fontSize: 12 }}>
                      <Text type="secondary">队列中</Text>
                      <Text strong>{stats.queuedJobCount}</Text>
                      <Text type="secondary">|</Text>
                      <Text type="secondary">失败</Text>
                      <Text strong style={{ color: token.colorError }}>{stats.failedJobCount}</Text>
                    </Space>
                  </div>
                  <div style={{
                    width: 56,
                    height: 56,
                    borderRadius: 14,
                    background: 'linear-gradient(135deg, #5B8DEF 0%, #7C6EAF 100%)',
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'center',
                    boxShadow: '0 8px 20px rgba(91, 141, 239, 0.4), 0 0 30px rgba(91, 141, 239, 0.2)',
                    transition: 'transform 0.3s cubic-bezier(0.4, 0, 0.2, 1)',
                  }}
                    onMouseEnter={(e) => e.currentTarget.style.transform = 'scale(1.1) rotate(5deg)'}
                    onMouseLeave={(e) => e.currentTarget.style.transform = 'scale(1) rotate(0deg)'}
                  >
                    <RocketOutlined style={{ fontSize: 28, color: '#fff' }} />
                  </div>
                </div>
              </Card>
            </Col>

            {/* 已创建配方 */}
            <Col xs={24} sm={12} lg={6}>
              <Card
                style={kpiCardStyle}
                hoverable
                onClick={() => navigate('/workspace/electrolytes')}
                onMouseEnter={(e) => {
                  e.currentTarget.style.transform = 'translateY(-4px)';
                  e.currentTarget.style.boxShadow = isDark
                    ? '0 12px 24px rgba(0, 0, 0, 0.4), 0 0 40px rgba(240, 147, 251, 0.2)'
                    : '0 12px 24px rgba(0, 0, 0, 0.12), 0 0 40px rgba(240, 147, 251, 0.25)';
                }}
                onMouseLeave={(e) => {
                  e.currentTarget.style.transform = 'translateY(0)';
                  e.currentTarget.style.boxShadow = kpiCardStyle.boxShadow as string;
                }}
              >
                <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
                  <div style={{ flex: 1 }}>
                    <Text type="secondary" style={{ fontSize: 12, display: 'block', marginBottom: 8, fontWeight: 600 }}>
                      ELECTROLYTE RECIPES
                    </Text>
                    <div style={{ fontSize: 36, fontWeight: 800, color: '#13c2c2', marginBottom: 8, lineHeight: 1 }}>
                      {stats.electrolyteCount}
                    </div>
                    <Space size={4} style={{ fontSize: 12 }}>
                      <Text type="secondary">本周新增</Text>
                      <Text strong style={{ color: token.colorSuccess }}>+{stats.weekNewElectrolytes}</Text>
                    </Space>
                  </div>
                  <div style={{
                    width: 56,
                    height: 56,
                    borderRadius: 14,
                    background: 'linear-gradient(135deg, #f093fb 0%, #f5576c 100%)',
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'center',
                    boxShadow: '0 8px 20px rgba(240, 147, 251, 0.4), 0 0 30px rgba(240, 147, 251, 0.2)',
                    transition: 'transform 0.3s cubic-bezier(0.4, 0, 0.2, 1)',
                  }}
                    onMouseEnter={(e) => e.currentTarget.style.transform = 'scale(1.1) rotate(5deg)'}
                    onMouseLeave={(e) => e.currentTarget.style.transform = 'scale(1) rotate(0deg)'}
                  >
                    <ExperimentOutlined style={{ fontSize: 28, color: '#fff' }} />
                  </div>
                </div>
              </Card>
            </Col>

            {/* 项目数量 */}
            <Col xs={24} sm={12} lg={6}>
              <Card
                style={kpiCardStyle}
                hoverable
                onClick={() => navigate('/workspace/projects')}
                onMouseEnter={(e) => {
                  e.currentTarget.style.transform = 'translateY(-4px)';
                  e.currentTarget.style.boxShadow = isDark
                    ? '0 12px 24px rgba(0, 0, 0, 0.4), 0 0 40px rgba(79, 172, 254, 0.2)'
                    : '0 12px 24px rgba(0, 0, 0, 0.12), 0 0 40px rgba(79, 172, 254, 0.25)';
                }}
                onMouseLeave={(e) => {
                  e.currentTarget.style.transform = 'translateY(0)';
                  e.currentTarget.style.boxShadow = kpiCardStyle.boxShadow as string;
                }}
              >
                <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
                  <div style={{ flex: 1 }}>
                    <Text type="secondary" style={{ fontSize: 12, display: 'block', marginBottom: 8, fontWeight: 600 }}>
                      PROJECTS
                    </Text>
                    <div style={{ fontSize: 36, fontWeight: 800, color: '#9254de', marginBottom: 8, lineHeight: 1 }}>
                      {stats.projectCount}
                    </div>
                    <Space size={4} style={{ fontSize: 12 }}>
                      <Text type="secondary">活跃项目</Text>
                      <Text strong>{stats.projectCount}</Text>
                    </Space>
                  </div>
                  <div style={{
                    width: 56,
                    height: 56,
                    borderRadius: 14,
                    background: 'linear-gradient(135deg, #4facfe 0%, #00f2fe 100%)',
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'center',
                    boxShadow: '0 8px 20px rgba(79, 172, 254, 0.4), 0 0 30px rgba(79, 172, 254, 0.2)',
                    transition: 'transform 0.3s cubic-bezier(0.4, 0, 0.2, 1)',
                  }}
                    onMouseEnter={(e) => e.currentTarget.style.transform = 'scale(1.1) rotate(5deg)'}
                    onMouseLeave={(e) => e.currentTarget.style.transform = 'scale(1) rotate(0deg)'}
                  >
                    <ProjectOutlined style={{ fontSize: 28, color: '#fff' }} />
                  </div>
                </div>
              </Card>
            </Col>

            {/* 机时使用情况 */}
            <Col xs={24} sm={12} lg={6}>
              <Card
                style={kpiCardStyle}
                hoverable
                onMouseEnter={(e) => {
                  e.currentTarget.style.transform = 'translateY(-4px)';
                  e.currentTarget.style.boxShadow = isDark
                    ? '0 12px 24px rgba(0, 0, 0, 0.4), 0 0 40px rgba(250, 112, 154, 0.2)'
                    : '0 12px 24px rgba(0, 0, 0, 0.12), 0 0 40px rgba(250, 112, 154, 0.25)';
                }}
                onMouseLeave={(e) => {
                  e.currentTarget.style.transform = 'translateY(0)';
                  e.currentTarget.style.boxShadow = kpiCardStyle.boxShadow as string;
                }}
              >
                <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start' }}>
                  <div style={{ flex: 1 }}>
                    <Text type="secondary" style={{ fontSize: 12, display: 'block', marginBottom: 8, fontWeight: 600 }}>
                      CPU HOURS USED
                    </Text>
                    <div style={{ fontSize: 36, fontWeight: 800, color: '#fa8c16', marginBottom: 8, lineHeight: 1 }}>
                      {stats.usedCPUHours.toFixed(1)}
                    </div>
                    <Space size={4} style={{ fontSize: 12 }}>
                      <Text type="secondary">/ {stats.totalCPUHours} h</Text>
                    </Space>
                    <Progress
                      percent={Math.round((stats.usedCPUHours / stats.totalCPUHours) * 100)}
                      strokeColor={{
                        '0%': '#fa8c16',
                        '100%': '#ff4d4f',
                      }}
                      size="small"
                      style={{ marginTop: 8 }}
                    />
                  </div>
                  <div style={{
                    width: 56,
                    height: 56,
                    borderRadius: 14,
                    background: 'linear-gradient(135deg, #fa709a 0%, #fee140 100%)',
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'center',
                    boxShadow: '0 8px 20px rgba(250, 112, 154, 0.4), 0 0 30px rgba(250, 112, 154, 0.2)',
                    transition: 'transform 0.3s cubic-bezier(0.4, 0, 0.2, 1)',
                  }}
                    onMouseEnter={(e) => e.currentTarget.style.transform = 'scale(1.1) rotate(5deg)'}
                    onMouseLeave={(e) => e.currentTarget.style.transform = 'scale(1) rotate(0deg)'}
                  >
                    <ThunderboltOutlined style={{ fontSize: 28, color: '#fff' }} />
                  </div>
                </div>
              </Card>
            </Col>
          </Row>

          {/* 第二层：运行中任务摘要 */}
          {runningJobs.length > 0 && (
            <Card
              title={
                <Space>
                  <SyncOutlined spin style={{ color: token.colorPrimary }} />
                  <Text strong>运行中的任务</Text>
                  <Badge count={runningJobs.length} style={{ backgroundColor: token.colorPrimary }} />
                </Space>
              }
              style={{ marginBottom: 24, borderRadius: 12, boxShadow: isDark ? '0 4px 12px rgba(0, 0, 0, 0.3)' : '0 4px 12px rgba(0, 0, 0, 0.06)' }}
            >
              <Space direction="vertical" size={16} style={{ width: '100%' }}>
                {runningJobs.map(job => (
                  <div
                    key={job.id}
                    style={{
                      padding: 16,
                      background: token.colorBgContainer,
                      borderRadius: 8,
                      border: `1px solid ${token.colorBorder}`,
                      cursor: 'pointer',
                      transition: 'all 0.3s',
                    }}
                    onClick={() => navigate(`/workspace/liquid-electrolyte/md/${job.id}`)}
                    onMouseEnter={(e) => {
                      e.currentTarget.style.borderColor = token.colorPrimary;
                      e.currentTarget.style.boxShadow = `0 4px 12px ${isDark ? 'rgba(107, 154, 255, 0.2)' : 'rgba(91, 141, 239, 0.15)'}`;
                    }}
                    onMouseLeave={(e) => {
                      e.currentTarget.style.borderColor = token.colorBorder;
                      e.currentTarget.style.boxShadow = 'none';
                    }}
                  >
                    <Row gutter={16} align="middle">
                      <Col flex="auto">
                        <Space direction="vertical" size={4} style={{ width: '100%' }}>
                          <Space>
                            <Text strong style={{ fontSize: 14 }}>
                              {job.config?.job_name || `任务 #${job.id}`}
                            </Text>
                            <StatusTag status={job.status} />
                            {job.slurm_job_id && (
                              <Tag color="blue">Slurm: {job.slurm_job_id}</Tag>
                            )}
                          </Space>
                          <Text type="secondary" style={{ fontSize: 12 }}>
                            {job.config?.user_note && `备注: ${job.config.user_note} | `}配方 ID: {job.system_id} | 创建于 {dayjs(job.created_at).format('YYYY-MM-DD HH:mm')}
                          </Text>
                        </Space>
                      </Col>
                      <Col flex="200px">
                        <div style={{ textAlign: 'right' }}>
                          <Text type="secondary" style={{ fontSize: 12, display: 'block', marginBottom: 4 }}>
                            进度
                          </Text>
                          <Progress
                            percent={job.progress}
                            status={job.status === JobStatus.RUNNING ? 'active' : 'normal'}
                            strokeColor={{
                              '0%': token.colorPrimary,
                              '100%': token.colorSuccess,
                            }}
                            size="small"
                          />
                        </div>
                      </Col>
                    </Row>
                  </div>
                ))}
              </Space>
            </Card>
          )}

          {/* 第三层：趋势与结构图表 */}
          <Row gutter={[16, 16]} style={{ marginBottom: 24 }}>
            {/* 左侧：任务趋势图 */}
            <Col xs={24} lg={14}>
              <Card
                title={
                  <Space>
                    <LineChartOutlined style={{ color: token.colorPrimary }} />
                    <Text strong>任务趋势</Text>
                  </Space>
                }
                extra={
                  <Segmented
                    options={[
                      { label: '按任务数', value: 'count' },
                      { label: '按机时', value: 'hours' },
                    ]}
                    value={chartType}
                    onChange={(value) => setChartType(value as 'count' | 'hours')}
                    size="small"
                  />
                }
                style={{ borderRadius: 12, boxShadow: isDark ? '0 4px 12px rgba(0, 0, 0, 0.3)' : '0 4px 12px rgba(0, 0, 0, 0.06)', height: '100%' }}
              >
                {trendData.length > 0 ? (
                  <Column
                    data={trendData}
                    xField="date"
                    yField="value"
                    seriesField="type"
                    isGroup={true}
                    columnStyle={{
                      radius: [4, 4, 0, 0],
                    }}
                    color={[token.colorSuccess, token.colorError, token.colorPrimary]}
                    legend={{
                      position: 'top-right',
                    }}
                    yAxis={{
                      label: {
                        formatter: (v: string) => chartType === 'hours' ? `${v}h` : v,
                      },
                      title: {
                        text: chartType === 'hours' ? 'CPU 核时 (h)' : '任务数',
                        style: {
                          fontSize: 12,
                        },
                      },
                    }}
                    theme={isDark ? 'dark' : 'light'}
                    smooth={true}
                    height={300}
                  />
                ) : (
                  <Empty description="暂无数据" style={{ padding: '60px 0' }} />
                )}
              </Card>
            </Col>

            {/* 右侧：资源与配方结构 */}
            <Col xs={24} lg={10}>
              <Space direction="vertical" size={16} style={{ width: '100%' }}>
                {/* 任务状态分布 */}
                <Card
                  title={
                    <Space>
                      <PieChartOutlined style={{ color: token.colorPrimary }} />
                      <Text strong>任务状态分布</Text>
                    </Space>
                  }
                  style={{ borderRadius: 12, boxShadow: isDark ? '0 4px 12px rgba(0, 0, 0, 0.3)' : '0 4px 12px rgba(0, 0, 0, 0.06)' }}
                >
                  {statusDistribution.length > 0 ? (
                    <Pie
                      data={statusDistribution}
                      angleField="value"
                      colorField="type"
                      radius={0.8}
                      innerRadius={0.6}
                      label={{
                        type: 'inner',
                        offset: '-30%',
                        content: '{value}',
                        style: {
                          fontSize: 14,
                          textAlign: 'center',
                        },
                      }}
                      legend={{
                        position: 'bottom',
                      }}
                      theme={isDark ? 'dark' : 'light'}
                      color={[token.colorSuccess, token.colorPrimary, token.colorWarning, token.colorError, '#8c8c8c']}
                      height={200}
                    />
                  ) : (
                    <Empty description="暂无数据" style={{ padding: '40px 0' }} />
                  )}
                </Card>

                {/* 机时占用结构 */}
                <Card
                  title={
                    <Space>
                      <BarChartOutlined style={{ color: token.colorPrimary }} />
                      <Text strong>机时占用 Top 5</Text>
                    </Space>
                  }
                  style={{ borderRadius: 12, boxShadow: isDark ? '0 4px 12px rgba(0, 0, 0, 0.3)' : '0 4px 12px rgba(0, 0, 0, 0.06)' }}
                >
                  {projectDistribution.length > 0 ? (
                    <Pie
                      data={projectDistribution}
                      angleField="value"
                      colorField="type"
                      radius={0.8}
                      label={{
                        type: 'outer',
                        content: '{type} {value}h',
                      }}
                      legend={{
                        position: 'bottom',
                      }}
                      theme={isDark ? 'dark' : 'light'}
                      height={200}
                    />
                  ) : (
                    <Empty description="暂无数据" style={{ padding: '40px 0' }} />
                  )}
                </Card>
              </Space>
            </Col>
          </Row>

          {/* 第四层：最近任务列表 */}
          <Card
            title={
              <Space>
                <ClockCircleOutlined style={{ color: token.colorPrimary }} />
                <Text strong>最近任务</Text>
              </Space>
            }
            extra={
              <Space>
                <Segmented
                  options={[
                    { label: '全部', value: 'all' },
                    { label: '运行中', value: 'running' },
                    { label: '已完成', value: 'completed' },
                    { label: '失败', value: 'failed' },
                  ]}
                  value={taskFilter}
                  onChange={(value) => setTaskFilter(value as any)}
                  size="small"
                />
                <Button type="link" onClick={() => navigate('/workspace/liquid-electrolyte/md')}>
                  查看全部
                </Button>
              </Space>
            }
            style={{ borderRadius: 12, boxShadow: isDark ? '0 4px 12px rgba(0, 0, 0, 0.3)' : '0 4px 12px rgba(0, 0, 0, 0.06)' }}
          >
            {allJobs.length === 0 ? (
              <div style={{ textAlign: 'center', padding: '60px 0' }}>
                <RocketOutlined style={{ fontSize: 48, color: '#d9d9d9', marginBottom: 16 }} />
                <div>
                  <Text type="secondary">暂无任务</Text>
                </div>
              </div>
            ) : (
              <Table
                dataSource={
                  allJobs
                    .filter(job => {
                      if (taskFilter === 'all') return true;
                      if (taskFilter === 'running') return job.status === JobStatus.RUNNING || job.status === JobStatus.QUEUED;
                      if (taskFilter === 'completed') return job.status === JobStatus.COMPLETED;
                      if (taskFilter === 'failed') return job.status === JobStatus.FAILED;
                      return true;
                    })
                    .slice(0, 10)
                }
                rowKey="id"
                pagination={false}
                onRow={(record) => ({
                  onClick: () => navigate(`/workspace/liquid-electrolyte/md/${record.id}`),
                  style: { cursor: 'pointer' },
                })}
                columns={[
                  {
                    title: '任务 ID',
                    dataIndex: 'id',
                    key: 'id',
                    width: 100,
                    render: (id) => <Text strong>#{id}</Text>,
                  },
                  {
                    title: '配方 ID',
                    dataIndex: 'system_id',
                    key: 'system_id',
                    width: 100,
                  },
                  {
                    title: '状态',
                    dataIndex: 'status',
                    key: 'status',
                    width: 120,
                    render: (status) => <StatusTag status={status} />,
                  },
                  {
                    title: '进度',
                    dataIndex: 'progress',
                    key: 'progress',
                    width: 150,
                    render: (progress, record) => (
                      <Progress
                        percent={progress}
                        status={record.status === JobStatus.RUNNING ? 'active' : 'normal'}
                        size="small"
                      />
                    ),
                  },
                  {
                    title: 'Slurm Job ID',
                    dataIndex: 'slurm_job_id',
                    key: 'slurm_job_id',
                    width: 120,
                    render: (id) => id || '-',
                  },
                  {
                    title: '提交时间',
                    dataIndex: 'created_at',
                    key: 'created_at',
                    width: 180,
                    render: (time) => dayjs(time).format('YYYY-MM-DD HH:mm:ss'),
                  },
                ]}
              />
            )}
          </Card>
        </Spin>
      </div>
    </div>
  );
}

