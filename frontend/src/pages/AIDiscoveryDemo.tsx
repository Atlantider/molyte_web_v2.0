import React, { useState } from 'react';
import { Card, Button, Input, Space, Typography, message, Spin } from 'antd';
import { ExperimentOutlined } from '@ant-design/icons';

const { Title, Text } = Typography;

const AIDiscoveryDemo: React.FC = () => {
  const [loading, setLoading] = useState(false);
  const [smiles, setSmiles] = useState('CCO');
  const [results, setResults] = useState<any>(null);

  const testPrediction = async () => {
    setLoading(true);
    try {
      // 直接调用后端API，不使用认证
      const response = await fetch('/api/v1/ai-discovery/predict-batch', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          smiles_list: [smiles],
          properties: ['bp', 'mp', 'fp', 'alpha', 'mu', 'gap', 'homo', 'lumo']
        })
      });

      if (response.ok) {
        const data = await response.json();
        setResults(data);
        message.success('预测成功！');
      } else {
        const errorData = await response.json();
        message.error(`预测失败: ${errorData.detail || '未知错误'}`);
      }
    } catch (error) {
      console.error('预测失败:', error);
      message.error('网络错误，请重试');
    } finally {
      setLoading(false);
    }
  };

  const testSimilarMolecules = async () => {
    setLoading(true);
    try {
      const response = await fetch(`/api/v1/ai-discovery/similar-molecules?smiles=${encodeURIComponent(smiles)}&n_similar=10&scope=global`);
      
      if (response.ok) {
        const data = await response.json();
        setResults(data);
        message.success('搜索成功！');
      } else {
        const errorData = await response.json();
        message.error(`搜索失败: ${errorData.detail || '未知错误'}`);
      }
    } catch (error) {
      console.error('搜索失败:', error);
      message.error('网络错误，请重试');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div style={{ padding: '24px', maxWidth: '1200px', margin: '0 auto' }}>
      <Card>
        <Title level={2}>
          <ExperimentOutlined /> AI Discovery 功能测试
        </Title>
        
        <Space direction="vertical" style={{ width: '100%' }} size="large">
          <div>
            <Text strong>输入SMILES:</Text>
            <Input
              value={smiles}
              onChange={(e) => setSmiles(e.target.value)}
              placeholder="输入SMILES字符串，如: CCO"
              style={{ marginTop: '8px' }}
            />
          </div>

          <Space>
            <Button 
              type="primary" 
              onClick={testPrediction}
              loading={loading}
            >
              测试属性预测
            </Button>
            <Button 
              onClick={testSimilarMolecules}
              loading={loading}
            >
              测试相似分子搜索
            </Button>
          </Space>

          {results && (
            <Card title="测试结果" style={{ marginTop: '16px' }}>
              <pre style={{ 
                background: '#f5f5f5', 
                padding: '16px', 
                borderRadius: '4px',
                overflow: 'auto',
                maxHeight: '400px'
              }}>
                {JSON.stringify(results, null, 2)}
              </pre>
            </Card>
          )}
        </Space>
      </Card>
    </div>
  );
};

export default AIDiscoveryDemo;
