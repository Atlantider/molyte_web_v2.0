/**
 * 分子相似性搜索组件
 * 支持基于中心分子的加权相似性搜索
 */
import React, { useState } from 'react';
import {
  Card, Input, Button, Slider, Row, Col, Table, Tag, Space, 
  Typography, Divider, Alert, Spin, Tooltip, Switch, message
} from 'antd';
import { SearchOutlined, SettingOutlined, InfoCircleOutlined } from '@ant-design/icons';
import apiClient from '../api/client';

const { Title, Text } = Typography;

interface PropertyWeights {
  melting_point: number;
  boiling_point: number;
  flash_point: number;
  dipole_moment: number;
  polarizability: number;
  homo: number;
  lumo: number;
  gap: number;
}

interface MoleculeProperties {
  melting_point?: number;
  boiling_point?: number;
  flash_point?: number;
  dipole_moment?: number;
  polarizability?: number;
  homo?: number;
  lumo?: number;
  gap?: number;
  real_values: Record<string, boolean>;
}

interface SimilarMolecule {
  smiles: string;
  name?: string;
  similarity_score: number;
  weighted_distance: number;
  properties: MoleculeProperties;
  is_center: boolean;
}

interface SimilaritySearchResponse {
  center_molecule: SimilarMolecule;
  similar_molecules: SimilarMolecule[];
  weights_used: PropertyWeights;
  total_found: number;
}

const SimilaritySearch: React.FC = () => {
  const [centerSmiles, setCenterSmiles] = useState('CC(=O)O'); // 默认醋酸
  const [loading, setLoading] = useState(false);
  const [results, setResults] = useState<SimilaritySearchResponse | null>(null);
  const [showWeights, setShowWeights] = useState(false);
  const [nSimilar, setNSimilar] = useState(20);
  const [includePredictions, setIncludePredictions] = useState(true);
  
  // 权重配置
  const [weights, setWeights] = useState<PropertyWeights>({
    melting_point: 1.0,
    boiling_point: 1.0,
    flash_point: 1.0,
    dipole_moment: 1.0,
    polarizability: 1.0,
    homo: 1.0,
    lumo: 1.0,
    gap: 1.0,
  });

  const handleSearch = async () => {
    if (!centerSmiles.trim()) {
      message.error('请输入中心分子的SMILES');
      return;
    }

    setLoading(true);
    try {
      const response = await apiClient.post('/similarity/search', {
        center_smiles: centerSmiles.trim(),
        n_similar: nSimilar,
        weights: weights,
        include_predictions: includePredictions,
      });
      setResults(response.data);
    } catch (error: any) {
      message.error(error.response?.data?.detail || '搜索失败');
    } finally {
      setLoading(false);
    }
  };

  const handleWeightChange = (property: keyof PropertyWeights, value: number) => {
    setWeights(prev => ({
      ...prev,
      [property]: value,
    }));
  };

  const loadPreset = async (presetName: string) => {
    try {
      const response = await apiClient.get('/similarity/weights/presets');
      const preset = response.data[presetName];
      if (preset) {
        setWeights(preset);
        message.success(`已加载${presetName}权重预设`);
      }
    } catch (error) {
      message.error('加载预设失败');
    }
  };

  const propertyLabels = {
    melting_point: '熔点',
    boiling_point: '沸点',
    flash_point: '闪点',
    dipole_moment: '偶极矩',
    polarizability: '极化率',
    homo: 'HOMO',
    lumo: 'LUMO',
    gap: 'HOMO-LUMO能隙',
  };

  const columns = [
    {
      title: 'SMILES',
      dataIndex: 'smiles',
      key: 'smiles',
      width: 200,
      render: (smiles: string, record: SimilarMolecule) => (
        <div>
          <Text code style={{ fontSize: 12 }}>{smiles}</Text>
          {record.is_center && <Tag color="gold" style={{ marginLeft: 8 }}>中心分子</Tag>}
        </div>
      ),
    },
    {
      title: '分子名称',
      dataIndex: 'name',
      key: 'name',
      width: 150,
      render: (name: string) => name || '-',
    },
    {
      title: '相似度',
      dataIndex: 'similarity_score',
      key: 'similarity_score',
      width: 100,
      render: (score: number) => (
        <Tag color={score > 0.8 ? 'green' : score > 0.6 ? 'orange' : 'red'}>
          {(score * 100).toFixed(1)}%
        </Tag>
      ),
      sorter: (a: SimilarMolecule, b: SimilarMolecule) => b.similarity_score - a.similarity_score,
    },
    {
      title: '熔点 (°C)',
      key: 'melting_point',
      width: 100,
      render: (record: SimilarMolecule) => {
        const value = record.properties.melting_point;
        const isReal = record.properties.real_values.melting_point;
        return value !== undefined ? (
          <Tooltip title={isReal ? '真实值' : '预测值'}>
            <Text type={isReal ? 'success' : 'secondary'}>
              {value.toFixed(1)}
            </Text>
          </Tooltip>
        ) : '-';
      },
    },
    {
      title: '沸点 (°C)',
      key: 'boiling_point',
      width: 100,
      render: (record: SimilarMolecule) => {
        const value = record.properties.boiling_point;
        const isReal = record.properties.real_values.boiling_point;
        return value !== undefined ? (
          <Tooltip title={isReal ? '真实值' : '预测值'}>
            <Text type={isReal ? 'success' : 'secondary'}>
              {value.toFixed(1)}
            </Text>
          </Tooltip>
        ) : '-';
      },
    },
    {
      title: '闪点 (°C)',
      key: 'flash_point',
      width: 100,
      render: (record: SimilarMolecule) => {
        const value = record.properties.flash_point;
        const isReal = record.properties.real_values.flash_point;
        return value !== undefined ? (
          <Tooltip title={isReal ? '真实值' : '预测值'}>
            <Text type={isReal ? 'success' : 'secondary'}>
              {value.toFixed(1)}
            </Text>
          </Tooltip>
        ) : '-';
      },
    },
    {
      title: '偶极矩',
      key: 'dipole_moment',
      width: 100,
      render: (record: SimilarMolecule) => {
        const value = record.properties.dipole_moment;
        return value !== undefined ? value.toFixed(2) : '-';
      },
    },
  ];

  return (
    <div style={{ padding: 24 }}>
      <Title level={3}>分子相似性搜索</Title>
      <Text type="secondary">
        基于中心分子搜索相似分子，支持多属性加权相似性计算
      </Text>

      <Card style={{ marginTop: 16 }}>
        <Row gutter={16} align="middle">
          <Col span={8}>
            <Text strong>中心分子 SMILES:</Text>
            <Input
              value={centerSmiles}
              onChange={(e) => setCenterSmiles(e.target.value)}
              placeholder="输入SMILES，如 CC(=O)O"
              style={{ marginTop: 8 }}
            />
          </Col>
          <Col span={4}>
            <Text strong>相似分子数量:</Text>
            <Slider
              min={5}
              max={50}
              value={nSimilar}
              onChange={setNSimilar}
              style={{ marginTop: 8 }}
            />
            <Text type="secondary">{nSimilar}</Text>
          </Col>
          <Col span={4}>
            <Text strong>包含预测值:</Text>
            <br />
            <Switch
              checked={includePredictions}
              onChange={setIncludePredictions}
              style={{ marginTop: 8 }}
            />
          </Col>
          <Col span={4}>
            <Button
              type="primary"
              icon={<SearchOutlined />}
              onClick={handleSearch}
              loading={loading}
              style={{ marginTop: 24 }}
            >
              搜索
            </Button>
          </Col>
          <Col span={4}>
            <Button
              icon={<SettingOutlined />}
              onClick={() => setShowWeights(!showWeights)}
              style={{ marginTop: 24 }}
            >
              权重设置
            </Button>
          </Col>
        </Row>

        {showWeights && (
          <>
            <Divider />
            <Title level={5}>属性权重配置</Title>
            <Row gutter={16}>
              <Col span={6}>
                <Space direction="vertical" size="small">
                  <Button size="small" onClick={() => loadPreset('thermal_properties')}>
                    热力学性质优先
                  </Button>
                  <Button size="small" onClick={() => loadPreset('electronic_properties')}>
                    电子性质优先
                  </Button>
                  <Button size="small" onClick={() => loadPreset('balanced')}>
                    平衡权重
                  </Button>
                </Space>
              </Col>
              <Col span={18}>
                <Row gutter={[16, 16]}>
                  {Object.entries(propertyLabels).map(([key, label]) => (
                    <Col span={6} key={key}>
                      <Text strong>{label}:</Text>
                      <Slider
                        min={0}
                        max={5}
                        step={0.1}
                        value={weights[key as keyof PropertyWeights]}
                        onChange={(value) => handleWeightChange(key as keyof PropertyWeights, value)}
                        style={{ marginTop: 4 }}
                      />
                      <Text type="secondary">{weights[key as keyof PropertyWeights].toFixed(1)}</Text>
                    </Col>
                  ))}
                </Row>
              </Col>
            </Row>
          </>
        )}
      </Card>

      {results && (
        <Card style={{ marginTop: 16 }}>
          <Title level={4}>搜索结果</Title>
          <Alert
            message={
              <div>
                找到 {results.total_found} 个相似分子
                <Tooltip title="绿色表示真实值，灰色表示预测值">
                  <InfoCircleOutlined style={{ marginLeft: 8 }} />
                </Tooltip>
              </div>
            }
            type="info"
            style={{ marginBottom: 16 }}
          />
          
          <Spin spinning={loading}>
            <Table
              columns={columns}
              dataSource={[results.center_molecule, ...results.similar_molecules]}
              rowKey="smiles"
              pagination={{ pageSize: 20 }}
              size="small"
              scroll={{ x: 1000 }}
            />
          </Spin>
        </Card>
      )}
    </div>
  );
};

export default SimilaritySearch;
