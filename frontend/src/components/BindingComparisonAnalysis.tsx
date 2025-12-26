/**
 * Binding 对比分析组件
 * 
 * 展示：
 * 1. 结构间对比（柱状图）
 * 2. 配体类型对比
 * 3. 趋势分析
 */
import React, { useMemo } from 'react';
import {
  Card,
  Row,
  Col,
  Table,
  Space,
  Tag,
  Empty,
  Tooltip,
  Divider,
  Typography,
  theme,
  Progress,
  Statistic,
  Select,
  Button,
} from 'antd';
import {
  DownloadOutlined,
  LineChartOutlined,
  SortAscendingOutlined,
} from '@ant-design/icons';
import type { ColumnsType } from 'antd/es/table';

const { Text, Title } = Typography;

interface ClusterResult {
  structure_id: number;
  composition: string;
  e_cluster_au: number;
  e_ion_au: number;
  total_ligand_energy_au: number;
  ligand_details: Record<string, { count: number; energy_au: number }>;
  e_bind_au: number;
  e_bind_ev: number;
  e_bind_kcal_mol: number;
}

interface Props {
  results?: ClusterResult[];
}

export default function BindingComparisonAnalysis({ results }: Props) {
  const { token } = theme.useToken();
  const [sortBy, setSortBy] = React.useState<'binding' | 'composition'>('binding');

  if (!results || results.length === 0) {
    return <Empty description="暂无对比数据" />;
  }

  const sortedResults = useMemo(() => {
    const sorted = [...results];
    if (sortBy === 'binding') {
      sorted.sort((a, b) => a.e_bind_kcal_mol - b.e_bind_kcal_mol);
    } else {
      sorted.sort((a, b) => a.composition.localeCompare(b.composition));
    }
    return sorted;
  }, [results, sortBy]);

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: 24 }}>
      {/* 对比分析工具栏 */}
      <Card style={{ boxShadow: '0 2px 8px rgba(0,0,0,0.08)' }}>
        <Row gutter={16} align="middle">
          <Col flex="auto">
            <Space>
              <LineChartOutlined style={{ fontSize: 18, color: token.colorPrimary }} />
              <Title level={5} style={{ margin: 0 }}>结构对比分析</Title>
            </Space>
          </Col>
          <Col>
            <Space>
              <Text type="secondary">排序方式：</Text>
              <Select
                value={sortBy}
                onChange={setSortBy}
                style={{ width: 150 }}
                options={[
                  { label: 'Binding 能量', value: 'binding' },
                  { label: '结构名称', value: 'composition' },
                ]}
              />
              <Button
                type="primary"
                ghost
                icon={<DownloadOutlined />}
                onClick={() => exportData(sortedResults)}
              >
                导出
              </Button>
            </Space>
          </Col>
        </Row>
      </Card>

      {/* 对比表格 */}
      <Card
        title={
          <Space>
            <SortAscendingOutlined style={{ fontSize: 20, color: token.colorPrimary }} />
            <Title level={5} style={{ margin: 0 }}>详细对比</Title>
          </Space>
        }
        style={{ boxShadow: '0 2px 8px rgba(0,0,0,0.08)' }}
      >
        <ComparisonTable results={sortedResults} />
      </Card>

      {/* 统计摘要 */}
      <Card
        title={
          <Space>
            <Title level={5} style={{ margin: 0 }}>统计摘要</Title>
          </Space>
        }
        style={{ boxShadow: '0 2px 8px rgba(0,0,0,0.08)' }}
      >
        <StatisticsSummary results={sortedResults} />
      </Card>
    </div>
  );
}

// 对比表格
function ComparisonTable({ results }: { results: ClusterResult[] }) {
  const { token } = theme.useToken();

  const columns: ColumnsType<ClusterResult> = [
    {
      title: '排名',
      key: 'rank',
      width: 60,
      align: 'center',
      render: (_, __, index) => (
        <Text strong style={{ color: token.colorPrimary }}>
          #{index + 1}
        </Text>
      ),
    },
    {
      title: '结构',
      dataIndex: 'composition',
      key: 'composition',
      render: (text) => <Text code>{text}</Text>,
    },
    {
      title: 'E_bind (kcal/mol)',
      dataIndex: 'e_bind_kcal_mol',
      key: 'e_bind_kcal_mol',
      align: 'right',
      render: (v: number) => {
        const isNegative = v < 0;
        return (
          <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'flex-end', gap: 8 }}>
            <Text strong style={{ color: isNegative ? token.colorSuccess : token.colorError, fontSize: 14 }}>
              {v.toFixed(2)}
            </Text>
            <Progress
              type="circle"
              percent={Math.abs(v) > 100 ? 100 : Math.abs(v)}
              width={40}
              strokeColor={isNegative ? token.colorSuccess : token.colorError}
            />
          </div>
        );
      },
    },
    {
      title: 'E_bind (eV)',
      dataIndex: 'e_bind_ev',
      key: 'e_bind_ev',
      align: 'right',
      render: (v: number) => <Text>{v.toFixed(4)}</Text>,
    },
    {
      title: '配体数',
      key: 'ligand_count',
      align: 'center',
      render: (_, record) => {
        const count = Object.values(record.ligand_details || {}).reduce((sum, d) => sum + d.count, 0);
        return <Tag color="blue">{count}</Tag>;
      },
    },
  ];

  return <Table columns={columns} dataSource={results} rowKey="structure_id" pagination={false} />;
}

// 统计摘要
function StatisticsSummary({ results }: { results: ClusterResult[] }) {
  const { token } = theme.useToken();

  const stats = useMemo(() => {
    const bindings = results.map(r => r.e_bind_kcal_mol);
    const mean = bindings.reduce((a, b) => a + b, 0) / bindings.length;
    const variance = bindings.reduce((a, b) => a + Math.pow(b - mean, 2), 0) / bindings.length;
    const std = Math.sqrt(variance);

    return {
      mean,
      std,
      min: Math.min(...bindings),
      max: Math.max(...bindings),
      median: bindings.sort((a, b) => a - b)[Math.floor(bindings.length / 2)],
      range: Math.max(...bindings) - Math.min(...bindings),
    };
  }, [results]);

  return (
    <Row gutter={[16, 16]}>
      <Col xs={24} sm={12} md={6}>
        <div style={{
          background: `linear-gradient(135deg, ${token.colorPrimary}15 0%, ${token.colorPrimary}05 100%)`,
          padding: 16,
          borderRadius: 8,
          border: `1px solid ${token.colorPrimary}20`,
        }}>
          <Statistic
            title="平均值"
            value={stats.mean}
            precision={2}
            suffix="kcal/mol"
            valueStyle={{ color: token.colorPrimary }}
          />
        </div>
      </Col>

      <Col xs={24} sm={12} md={6}>
        <div style={{
          background: `linear-gradient(135deg, ${token.colorWarning}15 0%, ${token.colorWarning}05 100%)`,
          padding: 16,
          borderRadius: 8,
          border: `1px solid ${token.colorWarning}20`,
        }}>
          <Statistic
            title="标准差"
            value={stats.std}
            precision={2}
            suffix="kcal/mol"
            valueStyle={{ color: token.colorWarning }}
          />
        </div>
      </Col>

      <Col xs={24} sm={12} md={6}>
        <div style={{
          background: `linear-gradient(135deg, ${token.colorSuccess}15 0%, ${token.colorSuccess}05 100%)`,
          padding: 16,
          borderRadius: 8,
          border: `1px solid ${token.colorSuccess}20`,
        }}>
          <Statistic
            title="最小值"
            value={stats.min}
            precision={2}
            suffix="kcal/mol"
            valueStyle={{ color: token.colorSuccess }}
          />
        </div>
      </Col>

      <Col xs={24} sm={12} md={6}>
        <div style={{
          background: `linear-gradient(135deg, ${token.colorError}15 0%, ${token.colorError}05 100%)`,
          padding: 16,
          borderRadius: 8,
          border: `1px solid ${token.colorError}20`,
        }}>
          <Statistic
            title="最大值"
            value={stats.max}
            precision={2}
            suffix="kcal/mol"
            valueStyle={{ color: token.colorError }}
          />
        </div>
      </Col>

      <Col xs={24} sm={12} md={6}>
        <div style={{
          background: `linear-gradient(135deg, ${token.colorInfo}15 0%, ${token.colorInfo}05 100%)`,
          padding: 16,
          borderRadius: 8,
          border: `1px solid ${token.colorInfo}20`,
        }}>
          <Statistic
            title="中位数"
            value={stats.median}
            precision={2}
            suffix="kcal/mol"
            valueStyle={{ color: token.colorInfo }}
          />
        </div>
      </Col>

      <Col xs={24} sm={12} md={6}>
        <div style={{
          background: `linear-gradient(135deg, #722ed1 15 0%, #722ed1 05 100%)`,
          padding: 16,
          borderRadius: 8,
          border: `1px solid #722ed120`,
        }}>
          <Statistic
            title="范围"
            value={stats.range}
            precision={2}
            suffix="kcal/mol"
            valueStyle={{ color: '#722ed1' }}
          />
        </div>
      </Col>
    </Row>
  );
}

// 导出数据
function exportData(results: ClusterResult[]) {
  const csv = [
    ['结构', 'E_bind (kcal/mol)', 'E_bind (eV)', 'E_bind (a.u.)', 'E_cluster (a.u.)', 'E_ion (a.u.)', 'E_ligand (a.u.)'],
    ...results.map(r => [
      r.composition,
      r.e_bind_kcal_mol.toFixed(2),
      r.e_bind_ev.toFixed(4),
      r.e_bind_au.toFixed(6),
      r.e_cluster_au.toFixed(6),
      r.e_ion_au.toFixed(6),
      r.total_ligand_energy_au.toFixed(6),
    ]),
  ].map(row => row.join(',')).join('\n');

  const blob = new Blob([csv], { type: 'text/csv' });
  const url = window.URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = `binding-results-${new Date().toISOString().split('T')[0]}.csv`;
  a.click();
}

