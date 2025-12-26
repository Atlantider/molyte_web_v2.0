/**
 * 电解质配方卡片组件
 */
import { Card, Space, Button, Popconfirm, Typography, Tag, Descriptions, Tooltip } from 'antd';
import { EditOutlined, DeleteOutlined, ExperimentOutlined, ThunderboltOutlined, CheckCircleOutlined, SyncOutlined, CloseCircleOutlined, CopyOutlined, CalendarOutlined } from '@ant-design/icons';
import type { ElectrolyteSystem, MDJob, JobStatus } from '../types';
import { UserRole } from '../types';
import dayjs from 'dayjs';
import { useAuthStore } from '../stores/authStore';
import { useThemeStore } from '../stores/themeStore';

const { Text } = Typography;

interface ElectrolyteCardProps {
  electrolyte: ElectrolyteSystem;
  jobs?: MDJob[];
  onEdit: (electrolyte: ElectrolyteSystem) => void;
  onCopy: (electrolyte: ElectrolyteSystem) => void;
  onDelete: (id: number) => void;
  onCreateJob: (electrolyte: ElectrolyteSystem) => void;
}

export default function ElectrolyteCard({
  electrolyte,
  jobs = [],
  onEdit,
  onCopy,
  onDelete,
  onCreateJob,
}: ElectrolyteCardProps) {
  const { user } = useAuthStore();
  const { mode } = useThemeStore();
  const isDark = mode === 'dark';

  // 解析组成信息（兼容新旧格式）
  const cations = electrolyte.cations || electrolyte.composition?.cations || [];
  const anions = electrolyte.anions || electrolyte.composition?.anions || [];
  const solvents = electrolyte.solvents || electrolyte.composition?.solvents || [];

  // 查找该配方的任务
  const relatedJobs = jobs.filter(job => job.system_id === electrolyte.id);
  const hasJobs = relatedJobs.length > 0;
  const latestJob = relatedJobs.length > 0 ? relatedJobs[0] : null;

  // 任务状态标签
  const getJobStatusTag = (status: JobStatus) => {
    const statusConfig: Record<string, { color: string; icon: React.ReactNode; text: string }> = {
      CREATED: { color: 'default', icon: <SyncOutlined />, text: '已创建' },
      QUEUED: { color: 'processing', icon: <SyncOutlined spin />, text: '排队中' },
      RUNNING: { color: 'processing', icon: <SyncOutlined spin />, text: '运行中' },
      POSTPROCESSING: { color: 'processing', icon: <SyncOutlined spin />, text: '后处理' },
      COMPLETED: { color: 'success', icon: <CheckCircleOutlined />, text: '已完成' },
      FAILED: { color: 'error', icon: <CloseCircleOutlined />, text: '失败' },
      CANCELLED: { color: 'warning', icon: <CloseCircleOutlined />, text: '已取消' },
    };
    const config = statusConfig[status] || statusConfig.CREATED;
    return (
      <Tag
        color={config.color}
        icon={config.icon}
        style={{ margin: 0, fontSize: '12px', fontWeight: 600 }}
      >
        {config.text}
      </Tag>
    );
  };

  // 处理卡片点击 - 整个卡片可点击
  const handleCardClick = (e: React.MouseEvent) => {
    // 如果点击的是按钮或其子元素，不触发卡片跳转
    const target = e.target as HTMLElement;
    if (target.closest('button') || target.closest('.ant-popconfirm') || target.closest('.ant-card-actions')) {
      return;
    }
    // 点击卡片时创建任务
    onCreateJob(electrolyte);
  };

  return (
    <Card
      hoverable
      onClick={handleCardClick}
      style={{
        height: '100%',
        cursor: 'pointer',
        transition: 'all 0.3s cubic-bezier(0.4, 0, 0.2, 1)',
        border: isDark ? '1px solid rgba(255, 255, 255, 0.12)' : '1px solid rgba(0, 0, 0, 0.06)',
        borderRadius: 16,
        background: isDark
          ? 'rgba(30, 30, 30, 0.8)'
          : 'rgba(255, 255, 255, 0.9)',
        backdropFilter: 'blur(12px)',
        WebkitBackdropFilter: 'blur(12px)',
        boxShadow: isDark
          ? '0 4px 12px rgba(0, 0, 0, 0.3)'
          : '0 4px 12px rgba(0, 0, 0, 0.08)',
      }}
      styles={{
        body: { padding: '24px' },
      }}
      onMouseEnter={(e) => {
        e.currentTarget.style.transform = 'translateY(-4px)';
        e.currentTarget.style.boxShadow = isDark
          ? '0 12px 24px rgba(0, 0, 0, 0.4), 0 0 40px rgba(240, 147, 251, 0.15)'
          : '0 12px 24px rgba(0, 0, 0, 0.12), 0 0 40px rgba(240, 147, 251, 0.2)';
      }}
      onMouseLeave={(e) => {
        e.currentTarget.style.transform = 'translateY(0)';
        e.currentTarget.style.boxShadow = isDark
          ? '0 4px 12px rgba(0, 0, 0, 0.3)'
          : '0 4px 12px rgba(0, 0, 0, 0.08)';
      }}
      actions={[
        <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'center', gap: '8px' }}>
          <Button
            type="link"
            icon={<ThunderboltOutlined />}
            onClick={(e) => { e.stopPropagation(); onCreateJob(electrolyte); }}
            style={{ padding: 0 }}
          >
            创建任务
          </Button>
          {hasJobs && (
            <div style={{ display: 'flex', alignItems: 'center', gap: '4px' }}>
              {getJobStatusTag(latestJob!.status)}
              {relatedJobs.length > 1 && (
                <span style={{ fontSize: '12px', color: '#999' }}>
                  +{relatedJobs.length - 1}
                </span>
              )}
            </div>
          )}
        </div>,
        <Button type="link" icon={<EditOutlined />} onClick={(e) => { e.stopPropagation(); onEdit(electrolyte); }}>
          编辑
        </Button>,
        <Button type="link" icon={<CopyOutlined />} onClick={(e) => { e.stopPropagation(); onCopy(electrolyte); }}>
          复制
        </Button>,
        <Popconfirm
          title="确定要删除这个电解质配方吗？"
          description="删除后将无法恢复，相关的计算任务也会被删除。"
          onConfirm={() => onDelete(electrolyte.id)}
          okText="确定"
          cancelText="取消"
        >
          <Button type="link" danger icon={<DeleteOutlined />} onClick={(e) => e.stopPropagation()}>
            删除
          </Button>
        </Popconfirm>,
      ]}
    >
      <div style={{ display: 'flex', gap: 20 }}>
        {/* 左侧图标 - 渐变 + 阴影 */}
        <div style={{
          width: 56,
          height: 56,
          borderRadius: 14,
          background: 'linear-gradient(135deg, #f093fb 0%, #f5576c 100%)',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          flexShrink: 0,
          boxShadow: '0 8px 20px rgba(240, 147, 251, 0.4), 0 0 30px rgba(240, 147, 251, 0.2)',
          transition: 'transform 0.3s cubic-bezier(0.4, 0, 0.2, 1)',
        }}
          onMouseEnter={(e) => e.currentTarget.style.transform = 'scale(1.1) rotate(5deg)'}
          onMouseLeave={(e) => e.currentTarget.style.transform = 'scale(1) rotate(0deg)'}
        >
          <ExperimentOutlined style={{ fontSize: 28, color: '#fff' }} />
        </div>

        {/* 右侧内容 */}
        <div style={{ flex: 1, minWidth: 0 }}>
          <Text strong style={{ fontSize: 16, display: 'block', marginBottom: 4 }}>
            {electrolyte.name}
          </Text>

          {/* 显示用户备注 */}
          {electrolyte.user_note && (
            <Text type="secondary" style={{ fontSize: 12, display: 'block', marginBottom: 8 }}>
              📝 {electrolyte.user_note}
            </Text>
          )}

          <div style={{ marginBottom: 12, display: 'flex', alignItems: 'center', gap: 4, flexWrap: 'wrap' }}>
            <Tag color="cyan" style={{ margin: 0 }}>{electrolyte.temperature} K</Tag>
            {/* 管理员可见：提交用户 */}
            {user?.role === UserRole.ADMIN && electrolyte.username && (
              <Tooltip title={`提交用户: ${electrolyte.user_email || '未知邮箱'}`}>
                <Tag color="purple" style={{ margin: 0 }}>
                  👤 {electrolyte.username}
                </Tag>
              </Tooltip>
            )}
          </div>

          <div style={{ marginBottom: 8 }}>
            <Text type="secondary" style={{ fontSize: 12, display: 'block', marginBottom: 4 }}>组成:</Text>
            <div style={{ display: 'flex', flexWrap: 'wrap', gap: 4 }}>
              {cations.slice(0, 2).map((c: any, i: number) => (
                <Tag key={`c-${i}`} color="blue" style={{ margin: 0, fontSize: 11 }}>
                  {c.name?.replace('+', '').replace('-', '')}×{c.number}
                </Tag>
              ))}
              {anions.slice(0, 2).map((a: any, i: number) => (
                <Tag key={`a-${i}`} color="red" style={{ margin: 0, fontSize: 11 }}>
                  {a.name?.replace('+', '').replace('-', '')}×{a.number}
                </Tag>
              ))}
              {solvents.slice(0, 3).map((s: any, i: number) => (
                <Tag key={`s-${i}`} color="green" style={{ margin: 0, fontSize: 11 }}>
                  {s.name}×{s.number}
                </Tag>
              ))}
              {(cations.length + anions.length + solvents.length > 7) && (
                <Tag style={{ margin: 0, fontSize: 11 }}>...</Tag>
              )}
            </div>
          </div>

          <Space size={4}>
            <CalendarOutlined style={{ color: '#bfbfbf', fontSize: 12 }} />
            <Text type="secondary" style={{ fontSize: 12 }}>
              {dayjs(electrolyte.created_at).format('YYYY-MM-DD HH:mm')}
            </Text>
          </Space>
        </div>
      </div>
    </Card>
  );
}

