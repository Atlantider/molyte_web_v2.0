/**
 * Binding 结果可视化组件
 * 
 * 展示：
 * 1. 总 Binding Energy 概览（多个结构的平均值）
 * 2. 分子-Li Pairwise Binding（各配体与 Li+ 的结合能）
 * 3. 能量分解详情（Cluster、离子、配体能量）
 * 4. 统计分析（均值、标准差、最值）
 */
import React, { useMemo } from 'react';
import {
  Card,
  Row,
  Col,
  Statistic,
  Table,
  Space,
  Tag,
  Progress,
  Empty,
  Tooltip,
  Divider,
  Typography,
  theme,
  Segmented,
} from 'antd';
import {
  ArrowUpOutlined,
  ArrowDownOutlined,
  ThunderboltOutlined,
  BgColorsOutlined,
  BarChartOutlined,
} from '@ant-design/icons';
import type { ColumnsType } from 'antd/es/table';

const { Text, Title } = Typography;

interface BindingTotalResult {
  e_ion_au?: number;
  ligand_energies_au?: Record<string, number>;
  cluster_binding_results?: Array<{
    structure_id: number;
    composition: string;
    e_cluster_au: number;
    e_ion_au: number;
    total_ligand_energy_au: number;
    ligand_details: Record<string, { count: number; energy_au: number }>;
    e_bind_au: number;
    e_bind_ev: number;
    e_bind_kcal_mol: number;
  }>;
  average_e_bind_au?: number;
  average_e_bind_ev?: number;
  average_e_bind_kcal_mol?: number;
}

interface PairwiseBinding {
  ligand: string;
  e_dimer_au: number;
  e_ligand_au: number;
  e_ion_au: number;
  e_bind_au: number;
  e_bind_ev: number;
  e_bind_kcal_mol: number;
}

interface Props {
  bindingTotalResult?: BindingTotalResult;
  pairwiseResult?: { pairwise_bindings?: PairwiseBinding[] };
}

export default function BindingResultsVisualization({ bindingTotalResult, pairwiseResult }: Props) {
  const { token } = theme.useToken();
  const [viewMode, setViewMode] = React.useState<'overview' | 'detailed'>('overview');

  // 计算统计数据
  const stats = useMemo(() => {
    if (!bindingTotalResult?.cluster_binding_results) return null;
    
    const results = bindingTotalResult.cluster_binding_results;
    const bindings = results.map(r => r.e_bind_kcal_mol);
    
    const mean = bindings.reduce((a, b) => a + b, 0) / bindings.length;
    const variance = bindings.reduce((a, b) => a + Math.pow(b - mean, 2), 0) / bindings.length;
    const std = Math.sqrt(variance);
    
    return {
      mean,
      std,
      min: Math.min(...bindings),
      max: Math.max(...bindings),
      count: bindings.length,
    };
  }, [bindingTotalResult]);

  if (!bindingTotalResult && !pairwiseResult) {
    return <Empty description="暂无 Binding 结果" />;
  }

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: 24 }}>
      {/* 总 Binding Energy 概览 */}
      {bindingTotalResult && (
        <Card
          title={
            <Space>
              <ThunderboltOutlined style={{ fontSize: 20, color: token.colorPrimary }} />
              <Title level={5} style={{ margin: 0 }}>总 Binding Energy</Title>
            </Space>
          }
          style={{ boxShadow: '0 2px 8px rgba(0,0,0,0.08)' }}
        >
          <Row gutter={[24, 24]}>
            {/* 主要指标 */}
            <Col xs={24} sm={12} md={6}>
              <div style={{
                background: `linear-gradient(135deg, ${token.colorPrimary}15 0%, ${token.colorPrimary}05 100%)`,
                padding: 16,
                borderRadius: 8,
                border: `1px solid ${token.colorPrimary}20`,
              }}>
                <Statistic
                  title="平均 Binding Energy"
                  value={bindingTotalResult.average_e_bind_kcal_mol || 0}
                  precision={2}
                  suffix="kcal/mol"
                  valueStyle={{ color: token.colorPrimary, fontSize: 24, fontWeight: 600 }}
                  prefix={bindingTotalResult.average_e_bind_kcal_mol! < 0 ? <ArrowDownOutlined /> : <ArrowUpOutlined />}
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
                  title="eV"
                  value={bindingTotalResult.average_e_bind_ev || 0}
                  precision={4}
                  valueStyle={{ color: token.colorSuccess, fontSize: 20 }}
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
                  title="Hartree"
                  value={bindingTotalResult.average_e_bind_au || 0}
                  precision={6}
                  valueStyle={{ color: token.colorWarning, fontSize: 20 }}
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
                  title="结构数"
                  value={stats?.count || 0}
                  valueStyle={{ color: token.colorInfo, fontSize: 20 }}
                />
              </div>
            </Col>
          </Row>

          {/* 统计分析 */}
          {stats && (
            <>
              <Divider />
              <Row gutter={16}>
                <Col xs={12} sm={6}>
                  <Text type="secondary">标准差</Text>
                  <div style={{ fontSize: 18, fontWeight: 600, marginTop: 4 }}>
                    {stats.std.toFixed(2)} kcal/mol
                  </div>
                </Col>
                <Col xs={12} sm={6}>
                  <Text type="secondary">最小值</Text>
                  <div style={{ fontSize: 18, fontWeight: 600, marginTop: 4, color: token.colorError }}>
                    {stats.min.toFixed(2)} kcal/mol
                  </div>
                </Col>
                <Col xs={12} sm={6}>
                  <Text type="secondary">最大值</Text>
                  <div style={{ fontSize: 18, fontWeight: 600, marginTop: 4, color: token.colorSuccess }}>
                    {stats.max.toFixed(2)} kcal/mol
                  </div>
                </Col>
                <Col xs={12} sm={6}>
                  <Text type="secondary">范围</Text>
                  <div style={{ fontSize: 18, fontWeight: 600, marginTop: 4 }}>
                    {(stats.max - stats.min).toFixed(2)} kcal/mol
                  </div>
                </Col>
              </Row>
            </>
          )}
        </Card>
      )}

      {/* 分子-Li Pairwise Binding */}
      {pairwiseResult?.pairwise_bindings && pairwiseResult.pairwise_bindings.length > 0 && (
        <Card
          title={
            <Space>
              <BgColorsOutlined style={{ fontSize: 20, color: token.colorPrimary }} />
              <Title level={5} style={{ margin: 0 }}>分子-Li Binding Energy</Title>
            </Space>
          }
          style={{ boxShadow: '0 2px 8px rgba(0,0,0,0.08)' }}
        >
          <PairwiseBindingTable bindings={pairwiseResult.pairwise_bindings} />
        </Card>
      )}

      {/* 详细结果表格 */}
      {bindingTotalResult?.cluster_binding_results && bindingTotalResult.cluster_binding_results.length > 0 && (
        <Card
          title={
            <Space>
              <BarChartOutlined style={{ fontSize: 20, color: token.colorPrimary }} />
              <Title level={5} style={{ margin: 0 }}>各结构详细结果</Title>
            </Space>
          }
          style={{ boxShadow: '0 2px 8px rgba(0,0,0,0.08)' }}
        >
          <ClusterBindingTable results={bindingTotalResult.cluster_binding_results} />
        </Card>
      )}
    </div>
  );
}

// 分子-Li Binding 表格
function PairwiseBindingTable({ bindings }: { bindings: PairwiseBinding[] }) {
  const columns: ColumnsType<PairwiseBinding> = [
    {
      title: '配体',
      dataIndex: 'ligand',
      key: 'ligand',
      render: (text) => <Tag color="blue">{text}</Tag>,
    },
    {
      title: 'E_bind (kcal/mol)',
      dataIndex: 'e_bind_kcal_mol',
      key: 'e_bind_kcal_mol',
      align: 'right',
      render: (v: number) => (
        <Text strong style={{ color: v < 0 ? '#52c41a' : '#ff4d4f' }}>
          {v.toFixed(2)}
        </Text>
      ),
    },
    {
      title: 'E_bind (eV)',
      dataIndex: 'e_bind_ev',
      key: 'e_bind_ev',
      align: 'right',
      render: (v: number) => v.toFixed(4),
    },
  ];

  return <Table columns={columns} dataSource={bindings} rowKey="ligand" pagination={false} />;
}

// Cluster Binding 表格
function ClusterBindingTable({ results }: { results: BindingTotalResult['cluster_binding_results'] }) {
  const columns: ColumnsType<any> = [
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
      render: (v: number) => (
        <Text strong style={{ color: v < 0 ? '#52c41a' : '#ff4d4f' }}>
          {v.toFixed(2)}
        </Text>
      ),
    },
    {
      title: 'E_bind (eV)',
      dataIndex: 'e_bind_ev',
      key: 'e_bind_ev',
      align: 'right',
      render: (v: number) => v.toFixed(4),
    },
  ];

  return <Table columns={columns} dataSource={results} rowKey="structure_id" pagination={false} />;
}

