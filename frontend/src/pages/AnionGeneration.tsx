/**
 * 溶盐力场构建页面
 */
import { useState, useEffect } from 'react';
import {
  Button,
  Space,
  message,
  Modal,
  Form,
  Input,
  Select,
  InputNumber,
  Spin,
  Empty,
  Card,
  Table,
  Tag,
  Tooltip,
  Row,
  Col,
  Divider,
  Alert,
  Progress,
  Descriptions,
  Typography,
  Statistic,
} from 'antd';
import {
  PlusOutlined,
  ReloadOutlined,
  DeleteOutlined,
  EyeOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  ClockCircleOutlined,
  ExperimentOutlined,
  ThunderboltOutlined,
} from '@ant-design/icons';
import type { ColumnsType } from 'antd/es/table';
import { useThemeStore } from '../stores/themeStore';
import {
  getAnionGenerationJobs,
  submitAnionGenerationJob,
  deleteAnionGenerationJob,
  type AnionGenerationJob as AnionGenerationJobType,
} from '../api/anion';

const { Title, Text } = Typography;
const { TextArea } = Input;

// 任务状态映射
const statusMap: Record<string, { color: string; text: string; icon: React.ReactNode }> = {
  pending: { color: 'default', text: '待处理', icon: <ClockCircleOutlined /> },
  running: { color: 'processing', text: '运行中', icon: <ExperimentOutlined /> },
  success: { color: 'success', text: '已完成', icon: <CheckCircleOutlined /> },
  failed: { color: 'error', text: '失败', icon: <CloseCircleOutlined /> },
};

interface AnionGenerationJob {
  id: number;
  job_id: string;
  anion_name: string;
  display_name: string;
  charge: number;
  identifier_type: string;
  identifier_value: string;
  status: string;
  message: string;
  lt_path?: string;
  pdb_path?: string;
  created_at: string;
  started_at?: string;
  finished_at?: string;
  cpu_hours_used?: number;
}

export default function AnionGeneration() {
  const [jobs, setJobs] = useState<AnionGenerationJob[]>([]);
  const [loading, setLoading] = useState(false);
  const [isModalVisible, setIsModalVisible] = useState(false);
  const [form] = Form.useForm();
  const [selectedJob, setSelectedJob] = useState<AnionGenerationJob | null>(null);
  const [detailModalVisible, setDetailModalVisible] = useState(false);
  const { isDark } = useThemeStore();
  const [autoRefreshEnabled, setAutoRefreshEnabled] = useState(false);

  // 获取任务列表
  const fetchJobs = async () => {
    setLoading(true);
    try {
      const data = await getAnionGenerationJobs();
      setJobs(data.jobs || []);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '获取任务列表失败');
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    fetchJobs();
  }, []);

  // 自动刷新逻辑（仅当启用时）
  useEffect(() => {
    if (!autoRefreshEnabled) return;

    const interval = setInterval(fetchJobs, 5000); // 每5秒刷新一次
    return () => clearInterval(interval);
  }, [autoRefreshEnabled]);

  // 提交新任务
  const handleSubmit = async (values: any) => {
    try {
      await submitAnionGenerationJob(values);
      message.success('任务提交成功');
      setIsModalVisible(false);
      form.resetFields();
      fetchJobs();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '任务提交失败');
    }
  };

  // 删除任务
  const handleDelete = async (jobId: number) => {
    try {
      await deleteAnionGenerationJob(jobId);
      message.success('任务已删除');
      fetchJobs();
    } catch (error: any) {
      message.error(error.response?.data?.detail || '删除任务失败');
    }
  };

  // 计算进度百分比
  const getProgress = (job: AnionGenerationJob): number => {
    const statusProgress: Record<string, number> = {
      pending: 0,
      running: 50,
      success: 100,
      failed: 0,
    };
    return statusProgress[job.status] || 0;
  };

  const columns: ColumnsType<AnionGenerationJob> = [
    {
      title: '阴离子名称',
      dataIndex: 'anion_name',
      key: 'anion_name',
      width: 100,
    },
    {
      title: '标识符',
      dataIndex: 'identifier_value',
      key: 'identifier_value',
      width: 150,
      render: (text) => <Text ellipsis>{text}</Text>,
    },
    {
      title: '状态',
      dataIndex: 'status',
      key: 'status',
      width: 100,
      render: (status) => {
        const config = statusMap[status] || statusMap.pending;
        return <Tag icon={config.icon} color={config.color}>{config.text}</Tag>;
      },
    },
    {
      title: '核时 (h)',
      dataIndex: 'cpu_hours_used',
      key: 'cpu_hours_used',
      width: 100,
      align: 'right' as const,
      render: (hours) => (
        <span style={{ display: 'flex', alignItems: 'center', justifyContent: 'flex-end', gap: '4px' }}>
          <ThunderboltOutlined style={{ color: '#faad14' }} />
          {hours !== undefined && hours !== null ? hours.toFixed(2) : '-'}
        </span>
      ),
    },
    {
      title: '进度',
      dataIndex: 'status',
      key: 'progress',
      width: 120,
      render: (status, record) => (
        <Progress percent={getProgress(record)} size="small" />
      ),
    },
    {
      title: '消息',
      dataIndex: 'message',
      key: 'message',
      width: 200,
      render: (text) => <Text ellipsis>{text}</Text>,
    },
    {
      title: '操作',
      key: 'action',
      width: 120,
      fixed: 'right' as const,
      render: (_, record) => (
        <Space>
          <Tooltip title="查看详情">
            <Button
              type="text"
              icon={<EyeOutlined />}
              onClick={() => {
                setSelectedJob(record);
                setDetailModalVisible(true);
              }}
            />
          </Tooltip>
          <Tooltip title="删除">
            <Button
              type="text"
              danger
              icon={<DeleteOutlined />}
              onClick={() => {
                Modal.confirm({
                  title: '确认删除',
                  content: `确定要删除任务 ${record.anion_name} 吗？`,
                  okText: '确定',
                  cancelText: '取消',
                  onOk: () => handleDelete(record.id),
                });
              }}
            />
          </Tooltip>
        </Space>
      ),
    },
  ];

  return (
    <div style={{ padding: '24px' }}>
      <Row gutter={[16, 16]}>
        {/* 统计卡片 */}
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic
              title="任务总数"
              value={jobs.length}
              suffix="个"
              prefix={<ExperimentOutlined />}
              valueStyle={{ color: '#1890ff' }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic
              title="已完成"
              value={jobs.filter(j => j.status === 'success').length}
              suffix="个"
              prefix={<CheckCircleOutlined />}
              valueStyle={{ color: '#52c41a' }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic
              title="运行中"
              value={jobs.filter(j => j.status === 'running').length}
              suffix="个"
              prefix={<ClockCircleOutlined />}
              valueStyle={{ color: '#faad14' }}
            />
          </Card>
        </Col>
        <Col xs={24} sm={12} md={6}>
          <Card>
            <Statistic
              title="总核时"
              value={jobs.reduce((sum, j) => sum + (j.cpu_hours_used || 0), 0).toFixed(2)}
              suffix="h"
              prefix={<ThunderboltOutlined />}
              valueStyle={{ color: '#ff7a45' }}
            />
          </Card>
        </Col>

        {/* 操作栏 */}
        <Col span={24}>
          <Card>
            <Row justify="space-between" align="middle">
              <Col>
                <div>
                  <Title level={2} style={{ margin: 0, marginBottom: 4 }}>
                    溶盐力场构建
                  </Title>
                  <Text type="secondary" style={{ fontSize: 14 }}>
                    溶盐力场参数自动生成任务管理
                  </Text>
                </div>
              </Col>
              <Col>
                <Space>
                  <Button
                    type={autoRefreshEnabled ? 'primary' : 'default'}
                    onClick={() => setAutoRefreshEnabled(!autoRefreshEnabled)}
                  >
                    {autoRefreshEnabled ? '⏸ 停止自动刷新' : '▶ 自动刷新'}
                  </Button>
                  <Button
                    type="primary"
                    icon={<PlusOutlined />}
                    onClick={() => setIsModalVisible(true)}
                  >
                    新建任务
                  </Button>
                  <Button
                    icon={<ReloadOutlined />}
                    onClick={fetchJobs}
                    loading={loading}
                  >
                    刷新
                  </Button>
                </Space>
              </Col>
            </Row>
          </Card>
        </Col>

        {/* 任务列表 */}
        <Col span={24}>
          <Card loading={loading}>
            {jobs.length === 0 ? (
              <Empty description="暂无任务" />
            ) : (
              <Table
                columns={columns}
                dataSource={jobs}
                rowKey="id"
                pagination={{ pageSize: 10 }}
                scroll={{ x: 1200 }}
              />
            )}
          </Card>
        </Col>
      </Row>

      {/* 新建任务模态框 */}
      <Modal
        title="新建溶盐力场开发任务"
        open={isModalVisible}
        onCancel={() => {
          setIsModalVisible(false);
          form.resetFields();
        }}
        footer={null}
      >
        <Form
          form={form}
          layout="vertical"
          onFinish={handleSubmit}
        >
          <Form.Item
            label="阴离子名称"
            name="anion_name"
            rules={[{ required: true, message: '请输入阴离子名称' }]}
          >
            <Input placeholder="例如: Cl, FSI, TFSI" />
          </Form.Item>

          <Form.Item
            label="显示名称"
            name="display_name"
            rules={[{ required: true, message: '请输入显示名称' }]}
          >
            <Input placeholder="例如: Chloride, Bis(fluorosulfonyl)imide" />
          </Form.Item>

          <Form.Item
            label="电荷"
            name="charge"
            initialValue={-1}
            rules={[{ required: true, message: '请输入电荷' }]}
          >
            <InputNumber />
          </Form.Item>

          <Form.Item
            label="标识符类型"
            name="identifier_type"
            initialValue="smiles"
            rules={[{ required: true, message: '请选择标识符类型' }]}
          >
            <Select>
              <Select.Option value="smiles">SMILES</Select.Option>
              <Select.Option value="cas">CAS号</Select.Option>
            </Select>
          </Form.Item>

          <Form.Item
            label="标识符值"
            name="identifier_value"
            rules={[{ required: true, message: '请输入标识符值' }]}
          >
            <TextArea
              placeholder="例如: [Cl-] (SMILES) 或 16887-00-6 (CAS号)"
              rows={3}
            />
          </Form.Item>

          <Form.Item>
            <Button type="primary" htmlType="submit" block>
              提交任务
            </Button>
          </Form.Item>
        </Form>
      </Modal>

      {/* 任务详情模态框 */}
      <Modal
        title={`任务详情: ${selectedJob?.anion_name}`}
        open={detailModalVisible}
        onCancel={() => setDetailModalVisible(false)}
        footer={null}
        width={700}
      >
        {selectedJob && (
          <Descriptions column={1} bordered>
            <Descriptions.Item label="任务ID">{selectedJob.job_id}</Descriptions.Item>
            <Descriptions.Item label="阴离子名称">{selectedJob.display_name}</Descriptions.Item>
            <Descriptions.Item label="标识符">{selectedJob.identifier_value}</Descriptions.Item>
            <Descriptions.Item label="电荷">{selectedJob.charge}</Descriptions.Item>
            <Descriptions.Item label="状态">
              <Tag color={statusMap[selectedJob.status]?.color}>
                {statusMap[selectedJob.status]?.text}
              </Tag>
            </Descriptions.Item>
            <Descriptions.Item label="进度">
              <Progress percent={getProgress(selectedJob)} />
            </Descriptions.Item>
            <Descriptions.Item label="核时消耗 (h)">
              <span style={{ display: 'flex', alignItems: 'center', gap: '6px' }}>
                <ThunderboltOutlined style={{ color: '#faad14' }} />
                {selectedJob.cpu_hours_used !== undefined && selectedJob.cpu_hours_used !== null
                  ? selectedJob.cpu_hours_used.toFixed(2)
                  : '-'}
              </span>
            </Descriptions.Item>
            <Descriptions.Item label="消息">{selectedJob.message}</Descriptions.Item>
            <Descriptions.Item label="创建时间">
              {new Date(selectedJob.created_at).toLocaleString()}
            </Descriptions.Item>
            {selectedJob.started_at && (
              <Descriptions.Item label="开始时间">
                {new Date(selectedJob.started_at).toLocaleString()}
              </Descriptions.Item>
            )}
            {selectedJob.finished_at && (
              <Descriptions.Item label="完成时间">
                {new Date(selectedJob.finished_at).toLocaleString()}
              </Descriptions.Item>
            )}
            {selectedJob.lt_path && (
              <Descriptions.Item label=".lt文件">{selectedJob.lt_path}</Descriptions.Item>
            )}
            {selectedJob.pdb_path && (
              <Descriptions.Item label=".pdb文件">{selectedJob.pdb_path}</Descriptions.Item>
            )}
          </Descriptions>
        )}
      </Modal>
    </div>
  );
}

