import React, { useState } from 'react';
import { Card, Typography, Space, Button, Collapse, Row, Col, Divider, Tag } from 'antd';
import { QuestionCircleOutlined, BookOutlined, ExperimentOutlined } from '@ant-design/icons';
import { useThemeStore } from '../stores/themeStore';

const { Title, Text, Paragraph } = Typography;
const { Panel } = Collapse;

const ParetoExplanation: React.FC = () => {
  const { mode } = useThemeStore();
  const isDark = mode === 'dark';
  const [showDetailed, setShowDetailed] = useState(false);

  return (
    <Card 
      title={
        <Space>
          <QuestionCircleOutlined style={{ color: '#1890ff' }} />
          <span>什么是Pareto最优？</span>
        </Space>
      }
      size="small"
      style={{ 
        marginBottom: '16px',
        border: `1px solid ${isDark ? '#434343' : '#d9d9d9'}`,
        background: isDark ? '#1f1f1f' : '#fafafa'
      }}
    >
      <Space direction="vertical" style={{ width: '100%' }} size="middle">
        {/* 核心概念 */}
        <div>
          <Text strong style={{ fontSize: '14px', color: '#1890ff' }}>
            🎯 核心概念
          </Text>
          <Paragraph style={{ marginTop: '8px', marginBottom: '12px' }}>
            <Text>
              <strong>Pareto最优解</strong>是指在多目标优化中，无法在不恶化其他目标的情况下改善任何一个目标的解。
              这些解构成了<strong>Pareto前沿</strong>，代表了所有可能的最优权衡方案。
            </Text>
          </Paragraph>
        </div>

        {/* 直观示例 */}
        <div style={{
          background: isDark ? '#262626' : '#f0f8ff',
          padding: '12px',
          borderRadius: '6px',
          border: `1px solid ${isDark ? '#434343' : '#91d5ff'}`
        }}>
          <Text strong style={{ color: '#1890ff' }}>💡 直观理解：</Text>
          <div style={{ marginTop: '8px' }}>
            <Text>
              假设我们要优化分子的<strong>能隙(gap)</strong>和<strong>闪点(fp)</strong>：
            </Text>
            <ul style={{ marginTop: '8px', paddingLeft: '20px' }}>
              <li><Text>分子A: gap=8.0, fp=150 ⭐</Text></li>
              <li><Text>分子B: gap=7.5, fp=180 ⭐</Text></li>
              <li><Text>分子C: gap=7.0, fp=160 ❌ (被B支配)</Text></li>
            </ul>
            <Text type="secondary" style={{ fontSize: '12px' }}>
              分子A和B都是Pareto最优的，因为A的gap更高，B的fp更高，无法简单判断哪个更好。
              分子C被B支配，因为B在两个目标上都不差于C，且至少在一个目标上更好。
            </Text>
          </div>
        </div>

        {/* 快速切换按钮 */}
        <div style={{ textAlign: 'center' }}>
          <Button 
            type="link" 
            size="small"
            onClick={() => setShowDetailed(!showDetailed)}
            icon={<BookOutlined />}
          >
            {showDetailed ? '收起详细说明' : '查看详细说明'}
          </Button>
        </div>

        {/* 详细说明 */}
        {showDetailed && (
          <Collapse ghost>
            <Panel header="📚 详细理论说明" key="theory">
              <Space direction="vertical" style={{ width: '100%' }} size="middle">
                
                {/* 数学定义 */}
                <div>
                  <Title level={5} style={{ color: '#52c41a' }}>🔢 数学定义</Title>
                  <div style={{
                    background: isDark ? '#1f1f1f' : '#f6ffed',
                    padding: '12px',
                    borderRadius: '6px',
                    border: `1px solid ${isDark ? '#434343' : '#b7eb8f'}`,
                    fontFamily: 'monospace'
                  }}>
                    <Text>
                      对于解x₁和x₂，如果满足：<br/>
                      • ∀i: f_i(x₁) ≥ f_i(x₂) （所有目标都不差）<br/>
                      • ∃j: f_j(x₁) &gt; f_j(x₂) （至少一个目标更好）<br/>
                      则称x₁<strong>支配</strong>x₂，记作x₁ ≻ x₂
                    </Text>
                  </div>
                  <Text type="secondary" style={{ fontSize: '12px' }}>
                    Pareto最优解是不被任何其他解支配的解
                  </Text>
                </div>

                {/* 应用场景 */}
                <div>
                  <Title level={5} style={{ color: '#fa8c16' }}>🧪 在分子设计中的应用</Title>
                  <Row gutter={16}>
                    <Col span={12}>
                      <div style={{
                        background: isDark ? '#2b1d16' : '#fff7e6',
                        padding: '10px',
                        borderRadius: '6px',
                        border: `1px solid ${isDark ? '#594214' : '#ffd591'}`
                      }}>
                        <Text strong>高温电解质设计</Text>
                        <ul style={{ marginTop: '6px', fontSize: '12px' }}>
                          <li>能隙 ↑ (电化学稳定性)</li>
                          <li>闪点 ↑ (安全性)</li>
                          <li>沸点 ↑ (热稳定性)</li>
                        </ul>
                      </div>
                    </Col>
                    <Col span={12}>
                      <div style={{
                        background: isDark ? '#161f2d' : '#f0f5ff',
                        padding: '10px',
                        borderRadius: '6px',
                        border: `1px solid ${isDark ? '#1c3a5c' : '#adc6ff'}`
                      }}>
                        <Text strong>宽电位窗口优化</Text>
                        <ul style={{ marginTop: '6px', fontSize: '12px' }}>
                          <li>HOMO ↓ (氧化稳定性)</li>
                          <li>LUMO ↑ (还原稳定性)</li>
                          <li>介电常数 ↑ (溶解性)</li>
                        </ul>
                      </div>
                    </Col>
                  </Row>
                </div>

                {/* 优势说明 */}
                <div>
                  <Title level={5} style={{ color: '#722ed1' }}>✨ 相比传统方法的优势</Title>
                  <Row gutter={16}>
                    <Col span={8}>
                      <Tag color="red">传统单目标</Tag>
                      <ul style={{ fontSize: '12px', marginTop: '6px' }}>
                        <li>只能优化一个属性</li>
                        <li>忽略其他重要指标</li>
                        <li>结果可能不实用</li>
                      </ul>
                    </Col>
                    <Col span={8}>
                      <Tag color="orange">加权求和</Tag>
                      <ul style={{ fontSize: '12px', marginTop: '6px' }}>
                        <li>需要预设权重</li>
                        <li>权重选择主观</li>
                        <li>可能错过最优解</li>
                      </ul>
                    </Col>
                    <Col span={8}>
                      <Tag color="green">Pareto优化</Tag>
                      <ul style={{ fontSize: '12px', marginTop: '6px' }}>
                        <li>找到所有最优解</li>
                        <li>提供权衡选择</li>
                        <li>客观数学标准</li>
                      </ul>
                    </Col>
                  </Row>
                </div>

                {/* 实际应用指导 */}
                <div>
                  <Title level={5} style={{ color: '#13c2c2' }}>🎯 如何使用Pareto结果</Title>
                  <div style={{
                    background: isDark ? '#1b2a2a' : '#e6fffb',
                    padding: '12px',
                    borderRadius: '6px',
                    border: `1px solid ${isDark ? '#2e5c5c' : '#87e8de'}`
                  }}>
                    <ol style={{ paddingLeft: '20px', margin: 0 }}>
                      <li><Text strong>分析Pareto前沿</Text> - 观察目标间的权衡关系</li>
                      <li><Text strong>识别关键区域</Text> - 找到性能平衡点</li>
                      <li><Text strong>结合实际需求</Text> - 根据应用场景选择最适合的解</li>
                      <li><Text strong>进一步优化</Text> - 在选定区域内精细调优</li>
                    </ol>
                  </div>
                </div>
              </Space>
            </Panel>
          </Collapse>
        )}

        {/* 图表说明 */}
        <Divider style={{ margin: '12px 0' }} />
        <div style={{ textAlign: 'center' }}>
          <Space>
            <div style={{ display: 'flex', alignItems: 'center' }}>
              <div style={{ 
                width: '12px', 
                height: '12px', 
                borderRadius: '50%', 
                background: '#ff4d4f',
                border: '2px solid #ffffff',
                marginRight: '6px' 
              }} />
              <Text style={{ fontSize: '12px' }}>Pareto最优解</Text>
            </div>
            <div style={{ display: 'flex', alignItems: 'center' }}>
              <div style={{ 
                width: '8px', 
                height: '8px', 
                borderRadius: '50%', 
                background: '#1890ff',
                opacity: 0.6,
                marginRight: '6px' 
              }} />
              <Text style={{ fontSize: '12px' }}>其他候选分子</Text>
            </div>
          </Space>
        </div>
      </Space>
    </Card>
  );
};

export default ParetoExplanation;
