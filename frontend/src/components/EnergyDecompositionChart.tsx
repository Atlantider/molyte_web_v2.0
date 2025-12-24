/**
 * 能量分解可视化组件
 * 
 * 展示：
 * 1. 能量分解瀑布图（Cluster → 离子 + 配体 → Binding）
 * 2. 配体贡献度分析
 * 3. 结构对比分析
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
} from 'antd';
import {
  InfoCircleOutlined,
  PercentageOutlined,
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
  ionEnergy?: number;
  ligandEnergies?: Record<string, number>;
}

export default function EnergyDecompositionChart({ results, ionEnergy, ligandEnergies }: Props) {
  const { token } = theme.useToken();

  if (!results || results.length === 0) {
    return <Empty description="暂无能量分解数据" />;
  }

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: 24 }}>
      {/* 能量分解详情表 */}
      <Card
        title={
          <Space>
            <PercentageOutlined style={{ fontSize: 20, color: token.colorPrimary }} />
            <Title level={5} style={{ margin: 0 }}>能量分解详情</Title>
          </Space>
        }
        style={{ boxShadow: '0 2px 8px rgba(0,0,0,0.08)' }}
      >
        <EnergyDecompositionTable results={results} />
      </Card>

      {/* 配体贡献度分析 */}
      <Card
        title={
          <Space>
            <InfoCircleOutlined style={{ fontSize: 20, color: token.colorPrimary }} />
            <Title level={5} style={{ margin: 0 }}>配体贡献度分析</Title>
          </Space>
        }
        style={{ boxShadow: '0 2px 8px rgba(0,0,0,0.08)' }}
      >
        <LigandContributionAnalysis results={results} />
      </Card>
    </div>
  );
}

// 能量分解表格
function EnergyDecompositionTable({ results }: { results: ClusterResult[] }) {
  const columns: ColumnsType<ClusterResult> = [
    {
      title: '结构',
      dataIndex: 'composition',
      key: 'composition',
      render: (text) => <Text code>{text}</Text>,
      width: 150,
    },
    {
      title: 'E_cluster (a.u.)',
      dataIndex: 'e_cluster_au',
      key: 'e_cluster_au',
      align: 'right',
      render: (v: number) => (
        <Tooltip title={`${(v * 627.509).toFixed(2)} kcal/mol`}>
          <Text>{v.toFixed(6)}</Text>
        </Tooltip>
      ),
    },
    {
      title: 'E_ion (a.u.)',
      dataIndex: 'e_ion_au',
      key: 'e_ion_au',
      align: 'right',
      render: (v: number) => (
        <Tooltip title={`${(v * 627.509).toFixed(2)} kcal/mol`}>
          <Text>{v.toFixed(6)}</Text>
        </Tooltip>
      ),
    },
    {
      title: 'E_ligand (a.u.)',
      dataIndex: 'total_ligand_energy_au',
      key: 'total_ligand_energy_au',
      align: 'right',
      render: (v: number) => (
        <Tooltip title={`${(v * 627.509).toFixed(2)} kcal/mol`}>
          <Text>{v.toFixed(6)}</Text>
        </Tooltip>
      ),
    },
    {
      title: 'E_bind (kcal/mol)',
      dataIndex: 'e_bind_kcal_mol',
      key: 'e_bind_kcal_mol',
      align: 'right',
      render: (v: number) => (
        <Text strong style={{ color: v < 0 ? '#52c41a' : '#ff4d4f', fontSize: 14 }}>
          {v.toFixed(2)}
        </Text>
      ),
    },
  ];

  return <Table columns={columns} dataSource={results} rowKey="structure_id" pagination={false} />;
}

// 配体贡献度分析
function LigandContributionAnalysis({ results }: { results: ClusterResult[] }) {
  const { token } = theme.useToken();

  // 计算平均配体贡献
  const ligandStats = useMemo(() => {
    const stats: Record<string, { count: number; totalEnergy: number; occurrences: number }> = {};

    results.forEach(result => {
      Object.entries(result.ligand_details || {}).forEach(([ligand, detail]) => {
        if (!stats[ligand]) {
          stats[ligand] = { count: 0, totalEnergy: 0, occurrences: 0 };
        }
        stats[ligand].count = detail.count;
        stats[ligand].totalEnergy += detail.energy_au;
        stats[ligand].occurrences += 1;
      });
    });

    return Object.entries(stats).map(([ligand, stat]) => ({
      ligand,
      count: stat.count,
      avgEnergy: stat.totalEnergy / stat.occurrences,
      totalEnergy: stat.totalEnergy,
      occurrences: stat.occurrences,
    }));
  }, [results]);

  const columns: ColumnsType<any> = [
    {
      title: '配体',
      dataIndex: 'ligand',
      key: 'ligand',
      render: (text) => <Tag color="blue">{text}</Tag>,
    },
    {
      title: '平均数量',
      dataIndex: 'count',
      key: 'count',
      align: 'center',
      render: (v: number) => <Text strong>{v}</Text>,
    },
    {
      title: '平均能量 (a.u.)',
      dataIndex: 'avgEnergy',
      key: 'avgEnergy',
      align: 'right',
      render: (v: number) => (
        <Tooltip title={`${(v * 627.509).toFixed(2)} kcal/mol`}>
          <Text>{v.toFixed(6)}</Text>
        </Tooltip>
      ),
    },
    {
      title: '出现次数',
      dataIndex: 'occurrences',
      key: 'occurrences',
      align: 'center',
      render: (v: number) => <Text type="secondary">{v}</Text>,
    },
  ];

  return (
    <div>
      <Table columns={columns} dataSource={ligandStats} rowKey="ligand" pagination={false} />
      
      <Divider />
      
      {/* 配体贡献度可视化 */}
      <div style={{ marginTop: 16 }}>
        <Text strong>配体能量分布</Text>
        <div style={{ marginTop: 12 }}>
          {ligandStats.map(stat => {
            const maxEnergy = Math.max(...ligandStats.map(s => Math.abs(s.avgEnergy)));
            const percentage = (Math.abs(stat.avgEnergy) / maxEnergy) * 100;
            
            return (
              <div key={stat.ligand} style={{ marginBottom: 12 }}>
                <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: 4 }}>
                  <Text>{stat.ligand}</Text>
                  <Text type="secondary">{stat.avgEnergy.toFixed(6)} a.u.</Text>
                </div>
                <Progress
                  percent={percentage}
                  strokeColor={stat.avgEnergy < 0 ? token.colorSuccess : token.colorError}
                  size="small"
                />
              </div>
            );
          })}
        </div>
      </div>
    </div>
  );
}

