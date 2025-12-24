import React, { useState } from 'react';
import { Card, Row, Col, Statistic, Space, Button, Tabs, Table, Tag, Progress, Typography, Tooltip } from 'antd';
import { LineChartOutlined, BarChartOutlined, TableOutlined, StarOutlined, QuestionCircleOutlined, CopyOutlined } from '@ant-design/icons';
import ParetoFrontierChart from './ParetoFrontierChart';
import TradeOffAnalysisChart from './TradeOffAnalysisChart';
import ParetoDrawer from './ParetoDrawer';
import { useThemeStore } from '../stores/themeStore';

const { Text, Title } = Typography;
const { TabPane } = Tabs;

interface MultiObjectiveDashboardProps {
  optimizeResults: any[];
  paretoFrontier: any[];
  objectiveWeights: { [key: string]: { direction: string; weight: number } };
  statistics?: any;
  baselineMolecule?: any;
}

const MultiObjectiveDashboard: React.FC<MultiObjectiveDashboardProps> = ({
  optimizeResults,
  paretoFrontier,
  objectiveWeights,
  statistics = {},
  baselineMolecule
}) => {
  const { mode } = useThemeStore();
  const isDark = mode === 'dark';
  const [activeTab, setActiveTab] = useState('pareto');
  const [showParetoDrawer, setShowParetoDrawer] = useState(false);

  // 获取活跃的优化目标
  const activeObjectives = Object.entries(objectiveWeights)
    .filter(([_, config]) => config.direction !== 'none')
    .map(([key, _]) => key);

  // 计算优化统计信息
  const totalCandidates = optimizeResults.length;
  const paretoCount = paretoFrontier.length;
  const improvementRate = totalCandidates > 0 ? (paretoCount / totalCandidates * 100) : 0;

  // 准备表格数据
  const tableColumns = [
    {
      title: '排名',
      dataIndex: 'rank',
      key: 'rank',
      width: 60,
      render: (text: any, record: any, index: number) => (
        <Space>
          {record.isParetoOptimal && <StarOutlined style={{ color: '#ff4d4f' }} />}
          {index + 1}
        </Space>
      )
    },
    {
      title: '分子信息',
      dataIndex: 'smiles',
      key: 'molecule_info',
      width: 280,
      render: (smiles: string, record: any) => {
        const moleculeInfo = record.molecule_info || {};
        const name = moleculeInfo.name || moleculeInfo.molecular_formula || '未知分子';
        const formula = moleculeInfo.molecular_formula;

        return (
          <div style={{ lineHeight: '1.2' }}>
            <div style={{ fontWeight: 'bold', fontSize: '13px', marginBottom: '2px' }}>
              <Tooltip title={name}>
                {name.length > 25 ? `${name.substring(0, 25)}...` : name}
              </Tooltip>
            </div>
            {formula && formula !== name && (
              <div style={{ fontSize: '11px', color: '#666', marginBottom: '2px' }}>
                <Tooltip title={formula}>
                  {formula.length > 30 ? `${formula.substring(0, 30)}...` : formula}
                </Tooltip>
              </div>
            )}
            <div style={{ display: 'flex', alignItems: 'center', gap: '4px' }}>
              <Tooltip title={`完整SMILES: ${smiles}`}>
                <Text
                  code
                  style={{
                    fontSize: '10px',
                    color: '#999',
                    maxWidth: '200px',
                    overflow: 'hidden',
                    textOverflow: 'ellipsis',
                    whiteSpace: 'nowrap',
                    display: 'inline-block'
                  }}
                >
                  {smiles}
                </Text>
              </Tooltip>
              <Tooltip title="复制SMILES">
                <Button
                  type="text"
                  size="small"
                  icon={<CopyOutlined />}
                  style={{ padding: '0 2px', height: '16px', minWidth: '16px' }}
                  onClick={() => {
                    navigator.clipboard.writeText(smiles);
                    // 可以添加成功提示
                  }}
                />
              </Tooltip>
            </div>
          </div>
        );
      }
    },
    ...activeObjectives.map(obj => ({
      title: obj.toUpperCase(),
      dataIndex: ['properties', obj],
      key: obj,
      width: 100,
      render: (value: number, record: any) => {
        const direction = objectiveWeights[obj]?.direction;
        const baselineValue = baselineMolecule?.properties?.[obj] || 0;
        const improvement = direction === 'up' ? value - baselineValue : baselineValue - value;
        const isImproved = improvement > 0;
        
        return (
          <Space direction="vertical" size={0}>
            <Text strong>{value?.toFixed(3) || '-'}</Text>
            <Text 
              style={{ 
                fontSize: '10px', 
                color: isImproved ? '#52c41a' : '#ff4d4f' 
              }}
            >
              {isImproved ? '↑' : '↓'} {Math.abs(improvement).toFixed(3)}
            </Text>
          </Space>
        );
      }
    })),
    {
      title: '相似度',
      dataIndex: 'similarity',
      key: 'similarity',
      width: 100,
      render: (similarity: number) => (
        <Space direction="vertical" size={0}>
          <Text>{(similarity * 100).toFixed(1)}%</Text>
          <Progress 
            percent={similarity * 100} 
            size="small" 
            showInfo={false}
            strokeColor={similarity > 0.7 ? '#52c41a' : similarity > 0.4 ? '#faad14' : '#ff4d4f'}
          />
        </Space>
      )
    },
    {
      title: '状态',
      key: 'status',
      width: 100,
      render: (record: any) => (
        <Tag color={record.isParetoOptimal ? 'red' : 'blue'}>
          {record.isParetoOptimal ? 'Pareto最优' : '候选分子'}
        </Tag>
      )
    }
  ];

  // 准备表格数据
  const tableData = optimizeResults.slice(0, 20).map((mol, index) => ({
    key: index,
    rank: index + 1,
    smiles: mol.smiles,
    properties: mol.properties || {},
    similarity: mol.similarity || 0,
    isParetoOptimal: paretoFrontier.some(p => p.smiles === mol.smiles),
    ...mol
  }));

  return (
    <div>
      {/* 概览统计 */}
      <Row gutter={16} style={{ marginBottom: '16px' }}>
        <Col span={6}>
          <Card>
            <Statistic
              title="总候选分子"
              value={totalCandidates}
              prefix="🧪"
              valueStyle={{ color: '#1890ff' }}
            />
          </Card>
        </Col>
        <Col span={6}>
          <Card>
            <Statistic
              title="Pareto最优解"
              value={paretoCount}
              prefix="⭐"
              valueStyle={{ color: '#ff4d4f' }}
            />
          </Card>
        </Col>
        <Col span={6}>
          <Card>
            <Statistic
              title="优化效率"
              value={improvementRate}
              precision={1}
              suffix="%"
              prefix="📈"
              valueStyle={{ color: '#52c41a' }}
            />
          </Card>
        </Col>
        <Col span={6}>
          <Card>
            <Statistic
              title="优化目标"
              value={activeObjectives.length}
              prefix="🎯"
              valueStyle={{ color: '#faad14' }}
            />
          </Card>
        </Col>
      </Row>

      {/* 详细分析标签页 */}
      <Card>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '16px' }}>
          <div></div>
          <Button
            type="primary"
            size="small"
            onClick={() => setShowParetoDrawer(true)}
            icon={<QuestionCircleOutlined />}
            style={{
              fontSize: '12px',
              background: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)',
              border: 'none',
              borderRadius: '6px',
              boxShadow: '0 2px 4px rgba(0,0,0,0.1)'
            }}
          >
            Pareto说明
          </Button>
        </div>
        <Tabs activeKey={activeTab} onChange={setActiveTab}>
          <TabPane 
            tab={
              <span>
                <LineChartOutlined />
                Pareto前沿
              </span>
            } 
            key="pareto"
          >
            <ParetoFrontierChart
              data={optimizeResults}
              paretoFrontier={paretoFrontier}
              objectives={activeObjectives}
              height={450}
            />
          </TabPane>

          <TabPane 
            tab={
              <span>
                <BarChartOutlined />
                权衡分析
              </span>
            } 
            key="tradeoff"
          >
            <TradeOffAnalysisChart
              data={optimizeResults}
              objectives={activeObjectives}
              height={400}
            />
          </TabPane>

          <TabPane 
            tab={
              <span>
                <TableOutlined />
                详细结果
              </span>
            } 
            key="table"
          >
            <div style={{ marginBottom: '16px' }}>
              <Space>
                <Text strong>优化结果列表</Text>
                <Tag color="red">⭐ = Pareto最优解</Tag>
                <Tag color="blue">🔵 = 候选分子</Tag>
              </Space>
            </div>
            <Table
              columns={tableColumns}
              dataSource={tableData}
              pagination={{
                pageSize: 10,
                showSizeChanger: true,
                showQuickJumper: true,
                showTotal: (total, range) => 
                  `第 ${range[0]}-${range[1]} 条，共 ${total} 条结果`
              }}
              scroll={{ x: 800 }}
              size="small"
            />
          </TabPane>

          <TabPane 
            tab={
              <span>
                <StarOutlined />
                优化摘要
              </span>
            } 
            key="summary"
          >
            <Row gutter={16}>
              <Col span={12}>
                <Card title="目标配置" size="small">
                  <Space direction="vertical" style={{ width: '100%' }}>
                    {Object.entries(objectiveWeights)
                      .filter(([_, config]) => config.direction !== 'none')
                      .map(([obj, config]) => (
                        <div key={obj} style={{ 
                          display: 'flex', 
                          justifyContent: 'space-between',
                          padding: '8px',
                          background: isDark ? '#1f1f1f' : '#fafafa',
                          borderRadius: '4px'
                        }}>
                          <Text strong>{obj.toUpperCase()}</Text>
                          <Space>
                            <Tag color={config.direction === 'up' ? 'green' : 'orange'}>
                              {config.direction === 'up' ? '提高 ↑' : '降低 ↓'}
                            </Tag>
                            <Text>权重: {config.weight}</Text>
                          </Space>
                        </div>
                      ))}
                  </Space>
                </Card>
              </Col>
              <Col span={12}>
                <Card title="优化统计" size="small">
                  <Space direction="vertical" style={{ width: '100%' }}>
                    <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                      <Text>搜索范围:</Text>
                      <Text strong>{statistics.scope || 'global'}</Text>
                    </div>
                    <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                      <Text>优化方法:</Text>
                      <Text strong>{statistics.optimization_method || 'pareto'}</Text>
                    </div>
                    <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                      <Text>Pareto效率:</Text>
                      <Text strong style={{ color: '#52c41a' }}>
                        {improvementRate.toFixed(1)}%
                      </Text>
                    </div>
                    <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                      <Text>平均相似度:</Text>
                      <Text strong>
                        {optimizeResults.length > 0 
                          ? (optimizeResults.reduce((sum, mol) => sum + (mol.similarity || 0), 0) / optimizeResults.length * 100).toFixed(1) + '%'
                          : 'N/A'
                        }
                      </Text>
                    </div>
                  </Space>
                </Card>
              </Col>
            </Row>
          </TabPane>
        </Tabs>
      </Card>

      {/* Pareto说明抽屉 */}
      <ParetoDrawer
        open={showParetoDrawer}
        onClose={() => setShowParetoDrawer(false)}
        showDetailed={true}
      />
    </div>
  );
};

export default MultiObjectiveDashboard;
