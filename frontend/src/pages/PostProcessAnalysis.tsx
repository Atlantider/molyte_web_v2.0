/**
 * 后处理分析列表页面
 * 展示所有 AdvancedClusterJob 任务，支持基于 MD Job 筛选
 */
import { useState, useEffect, useCallback } from 'react';
import { useNavigate, useSearchParams } from 'react-router-dom';
import {
  Card,
  Table,
  Button,
  Space,
  Tag,
  Typography,
  Select,
  Empty,
  Spin,
  Row,
  Col,
  Statistic,
  Progress,
  message,
  Tooltip,
  theme,
  Popconfirm,
} from 'antd';
import {
  PlusOutlined,
  ReloadOutlined,
  ExperimentOutlined,
  RocketOutlined,
  CheckCircleOutlined,
  ClockCircleOutlined,
  SyncOutlined,
  CloseCircleOutlined,
  ThunderboltOutlined,
  LineChartOutlined,
  DeleteOutlined,
  UserOutlined,
  EyeOutlined,
  StopOutlined,
  RedoOutlined,
} from '@ant-design/icons';
import type { ColumnsType } from 'antd/es/table';
import { listClusterAnalysisJobs, deleteClusterAnalysisJob, cancelClusterAnalysisJob, resubmitClusterAnalysisJob, CALC_TYPE_INFO, type AdvancedClusterJob, type ClusterCalcType } from '../api/clusterAnalysis';
import { getMDJobs } from '../api/jobs';
import type { MDJob } from '../types';
import { JobStatus } from '../types';
import { useAuthStore } from '../stores/authStore';
import dayjs from 'dayjs';
import relativeTime from 'dayjs/plugin/relativeTime';
import 'dayjs/locale/zh-cn';

dayjs.extend(relativeTime);
dayjs.locale('zh-cn');

const { Title, Text } = Typography;

// 状态配置
const STATUS_CONFIG: Record<string, { color: string; icon: React.ReactNode; text: string }> = {
  CREATED: { color: 'default', icon: <ClockCircleOutlined />, text: '已创建' },
  SUBMITTED: { color: 'processing', icon: <SyncOutlined spin />, text: '已提交' },
  RUNNING: { color: 'processing', icon: <SyncOutlined spin />, text: '运行中' },
  WAITING_QC: { color: 'warning', icon: <ThunderboltOutlined />, text: '等待QC' },
  CALCULATING: { color: 'processing', icon: <SyncOutlined spin />, text: '计算中' },
  COMPLETED: { color: 'success', icon: <CheckCircleOutlined />, text: '已完成' },
  FAILED: { color: 'error', icon: <CloseCircleOutlined />, text: '失败' },
  CANCELLED: { color: 'default', icon: <CloseCircleOutlined />, text: '已取消' },
};

export default function PostProcessAnalysis() {
  const navigate = useNavigate();
  const [searchParams] = useSearchParams();
  const { token } = theme.useToken();
  const { user } = useAuthStore();
  const isAdmin = user?.role === 'ADMIN';

  const [loading, setLoading] = useState(false);
  const [jobs, setJobs] = useState<AdvancedClusterJob[]>([]);
  const [mdJobs, setMdJobs] = useState<MDJob[]>([]);
  const [selectedMdJobId, setSelectedMdJobId] = useState<number | undefined>(
    searchParams.get('md_job_id') ? Number(searchParams.get('md_job_id')) : undefined
  );

  // 统计数据
  const stats = {
    total: jobs.length,
    running: jobs.filter(j => ['SUBMITTED', 'RUNNING', 'WAITING_QC', 'CALCULATING'].includes(j.status)).length,
    completed: jobs.filter(j => j.status === 'COMPLETED').length,
    failed: jobs.filter(j => j.status === 'FAILED').length,
    cancelled: jobs.filter(j => j.status === 'CANCELLED').length,
  };

  // 加载 MD Jobs（只加载已完成的）
  const loadMdJobs = useCallback(async () => {
    try {
      const data = await getMDJobs();
      setMdJobs(data.filter(j => j.status === JobStatus.COMPLETED));
    } catch (err) {
      console.error('Failed to load MD jobs:', err);
    }
  }, []);

  // 加载分析任务
  const loadJobs = useCallback(async () => {
    setLoading(true);
    try {
      const data = await listClusterAnalysisJobs(selectedMdJobId);
      setJobs(data);
    } catch (err) {
      console.error('Failed to load jobs:', err);
      message.error('加载分析任务失败');
    } finally {
      setLoading(false);
    }
  }, [selectedMdJobId]);

  // 删除任务
  const handleDelete = async (jobId: number) => {
    try {
      await deleteClusterAnalysisJob(jobId);
      message.success('任务已删除');
      loadJobs();
    } catch (err: any) {
      message.error(err.response?.data?.detail || '删除任务失败');
    }
  };

  // 取消任务
  const handleCancel = async (jobId: number) => {
    try {
      await cancelClusterAnalysisJob(jobId);
      message.success('任务已取消');
      loadJobs();
    } catch (err: any) {
      message.error(err.response?.data?.detail || '取消任务失败');
    }
  };

  // 重新提交任务
  const handleResubmit = async (jobId: number) => {
    try {
      await resubmitClusterAnalysisJob(jobId);
      message.success('任务已重新提交');
      loadJobs();
    } catch (err: any) {
      message.error(err.response?.data?.detail || '重新提交任务失败');
    }
  };

  useEffect(() => {
    loadMdJobs();
  }, [loadMdJobs]);

  useEffect(() => {
    loadJobs();
    // 只在有运行中的任务时才定时刷新（每 10 秒刷新一次）
    const hasRunningJobs = stats.running > 0;
    if (hasRunningJobs) {
      const interval = setInterval(() => {
        loadJobs();
      }, 10000);
      return () => clearInterval(interval);
    }
  }, [loadJobs, stats.running]);

  // 表格列定义
  const columns: ColumnsType<AdvancedClusterJob> = [
    {
      title: 'ID',
      dataIndex: 'id',
      width: 70,
      render: (id: number) => <Text strong>#{id}</Text>,
    },
    {
      title: 'MD Job',
      dataIndex: 'md_job_id',
      width: 100,
      render: (mdJobId: number) => (
        <Button
          type="link"
          size="small"
          icon={<RocketOutlined />}
          onClick={() => navigate(`/workspace/liquid-electrolyte/md/${mdJobId}`)}
        >
          #{mdJobId}
        </Button>
      ),
    },
    ...(isAdmin ? [{
      title: '用户',
      dataIndex: 'username',
      width: 120,
      render: (username: string, record: AdvancedClusterJob) => (
        <Tooltip title={record.user_email}>
          <Space size={4}>
            <UserOutlined />
            <Text>{username || `User #${record.user_id}`}</Text>
          </Space>
        </Tooltip>
      ),
    }] : []),
    {
      title: '计算类型',
      dataIndex: 'calc_types',
      width: 350,
      render: (calcTypes: string[]) => (
        <Space size={4} wrap={false}>
          {calcTypes.map((ct) => {
            const info = CALC_TYPE_INFO[ct as ClusterCalcType];
            return (
              <Tooltip key={ct} title={info?.description}>
                <Tag color={info?.riskLevel === 'high' ? 'red' : info?.riskLevel === 'medium' ? 'orange' : 'blue'}>
                  {info?.icon} {info?.label || ct}
                </Tag>
              </Tooltip>
            );
          })}
        </Space>
      ),
    },
    {
      title: '结构数',
      dataIndex: 'selected_structures',
      width: 80,
      align: 'center',
      render: (s: { count: number }) => <Text>{s?.count || 0}</Text>,
    },
    {
      title: '状态',
      dataIndex: 'status',
      width: 120,
      render: (status: string) => {
        const cfg = STATUS_CONFIG[status] || STATUS_CONFIG.CREATED;
        return (
          <Tag icon={cfg.icon} color={cfg.color}>{cfg.text}</Tag>
        );
      },
    },
    {
      title: 'QC 进度',
      dataIndex: 'qc_task_plan',
      width: 120,
      render: (plan: AdvancedClusterJob['qc_task_plan']) => {
        if (!plan) return <Text type="secondary">-</Text>;
        const completed = plan.completed_qc_tasks || 0;
        const total = plan.total_qc_tasks || 0;
        if (total === 0) return <Text type="secondary">-</Text>;
        return (
          <Tooltip title={`${completed}/${total} QC 任务已完成`}>
            <Progress
              percent={Math.round((completed / total) * 100)}
              size="small"
              status={completed === total ? 'success' : 'active'}
              format={() => `${completed}/${total}`}
            />
          </Tooltip>
        );
      },
    },
    {
      title: '核时',
      dataIndex: 'cpu_hours_used',
      width: 100,
      align: 'center',
      render: (cpuHours: number | undefined) => {
        if (cpuHours === undefined || cpuHours === null) {
          return <Text type="secondary">-</Text>;
        }
        return (
          <Tooltip title={`${cpuHours.toFixed(2)} 小时`}>
            <Text strong>{cpuHours.toFixed(1)}h</Text>
          </Tooltip>
        );
      },
    },
    {
      title: '创建时间',
      dataIndex: 'created_at',
      width: 140,
      render: (t: string) => (
        <Tooltip title={dayjs(t).format('YYYY-MM-DD HH:mm:ss')}>
          <Text type="secondary">{dayjs(t).fromNow()}</Text>
        </Tooltip>
      ),
    },
    {
      title: '操作',
      key: 'actions',
      width: 180,
      align: 'center',
      render: (_, record) => (
        <Space size="small" wrap>
          <Tooltip title="查看详情">
            <Button
              type="text"
              size="small"
              icon={<EyeOutlined />}
              onClick={() => navigate(`/workspace/liquid-electrolyte/analysis/${record.id}`)}
            />
          </Tooltip>
          {['CREATED', 'SUBMITTED', 'RUNNING', 'WAITING_QC', 'CALCULATING'].includes(record.status) && (
            <Popconfirm
              title="确定取消此任务？"
              description="取消后可重新提交"
              onConfirm={() => handleCancel(record.id)}
              okText="确定"
              cancelText="取消"
            >
              <Tooltip title="取消任务">
                <Button
                  type="text"
                  size="small"
                  icon={<StopOutlined />}
                />
              </Tooltip>
            </Popconfirm>
          )}
          {['FAILED', 'CANCELLED'].includes(record.status) && (
            <Tooltip title="重新提交任务">
              <Button
                type="text"
                size="small"
                icon={<RedoOutlined />}
                onClick={() => handleResubmit(record.id)}
              />
            </Tooltip>
          )}
          <Popconfirm
            title="确定删除此任务？"
            description="删除后无法恢复"
            onConfirm={() => handleDelete(record.id)}
            okText="确定"
            cancelText="取消"
          >
            <Tooltip title="删除任务">
              <Button
                type="text"
                danger
                size="small"
                icon={<DeleteOutlined />}
              />
            </Tooltip>
          </Popconfirm>
        </Space>
      ),
    },
  ];

  return (
    <div style={{ padding: 24 }}>
      {/* 页面标题 */}
      <div style={{ marginBottom: 24 }}>
        <Title level={4} style={{ margin: 0 }}>
          <LineChartOutlined style={{ marginRight: 8, color: token.colorPrimary }} />
          后处理分析
        </Title>
        <Text type="secondary">基于 MD 模拟结果进行 Binding/Desolvation/Redox/Reorg 计算</Text>
      </div>

      {/* 统计卡片 */}
      <Row gutter={8} style={{ marginBottom: 24 }}>
        <Col flex="1 1 20%">
          <Card size="small" style={{ padding: '8px' }}>
            <Statistic title="总任务数" value={stats.total} prefix={<ExperimentOutlined />} />
          </Card>
        </Col>
        <Col flex="1 1 20%">
          <Card size="small" style={{ padding: '8px' }}>
            <Statistic title="运行中" value={stats.running} valueStyle={{ color: '#1890ff' }} prefix={<SyncOutlined spin={stats.running > 0} />} />
          </Card>
        </Col>
        <Col flex="1 1 20%">
          <Card size="small" style={{ padding: '8px' }}>
            <Statistic title="已完成" value={stats.completed} valueStyle={{ color: '#52c41a' }} prefix={<CheckCircleOutlined />} />
          </Card>
        </Col>
        <Col flex="1 1 20%">
          <Card size="small" style={{ padding: '8px' }}>
            <Statistic title="失败" value={stats.failed} valueStyle={{ color: stats.failed > 0 ? '#ff4d4f' : undefined }} prefix={<CloseCircleOutlined />} />
          </Card>
        </Col>
        <Col flex="1 1 20%">
          <Card size="small" style={{ padding: '8px' }}>
            <Statistic title="已取消" value={stats.cancelled} valueStyle={{ color: stats.cancelled > 0 ? '#8c8c8c' : undefined }} prefix={<CloseCircleOutlined />} />
          </Card>
        </Col>
      </Row>

      {/* 工具栏 */}
      <Card
        style={{ marginBottom: 16 }}
        bodyStyle={{ padding: '12px 16px' }}
      >
        <Row justify="space-between" align="middle">
          <Col>
            <Space>
              <Text>数据来源：</Text>
              <Select
                style={{ width: 300 }}
                placeholder="选择 MD Job（可选）"
                allowClear
                value={selectedMdJobId}
                onChange={(v) => setSelectedMdJobId(v)}
                options={mdJobs.map(j => ({
                  value: j.id,
                  label: `#${j.id} - ${j.config?.job_name || 'MD Job'}`,
                }))}
              />
            </Space>
          </Col>
          <Col>
            <Space>
              <Button icon={<ReloadOutlined />} onClick={loadJobs}>刷新</Button>
              <Button
                type="primary"
                icon={<PlusOutlined />}
                onClick={() => navigate('/workspace/liquid-electrolyte/analysis/create')}
              >
                新建分析
              </Button>
            </Space>
          </Col>
        </Row>
      </Card>

      {/* 任务列表 */}
      <Card>
        <Spin spinning={loading}>
          {jobs.length === 0 && !loading ? (
            <Empty
              image={Empty.PRESENTED_IMAGE_SIMPLE}
              description={
                <Space direction="vertical">
                  <Text>暂无分析任务</Text>
                  <Button
                    type="primary"
                    icon={<PlusOutlined />}
                    onClick={() => navigate('/workspace/liquid-electrolyte/analysis/create')}
                  >
                    创建第一个分析任务
                  </Button>
                </Space>
              }
            />
          ) : (
            <Table
              columns={columns}
              dataSource={jobs}
              rowKey="id"
              pagination={{
                pageSize: 20,
                showSizeChanger: true,
                showTotal: (total) => `共 ${total} 条`,
              }}
              size="middle"
            />
          )}
        </Spin>
      </Card>
    </div>
  );
}

