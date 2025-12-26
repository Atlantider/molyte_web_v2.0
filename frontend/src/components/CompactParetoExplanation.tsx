import React, { useState } from 'react';
import { Card, Typography, Space, Button, Collapse, Tag, Divider } from 'antd';
import { QuestionCircleOutlined, BookOutlined, CloseOutlined } from '@ant-design/icons';
import { useThemeStore } from '../stores/themeStore';

const { Text, Paragraph } = Typography;
const { Panel } = Collapse;

interface CompactParetoExplanationProps {
  onClose?: () => void;
}

const CompactParetoExplanation: React.FC<CompactParetoExplanationProps> = ({ onClose }) => {
  const { mode } = useThemeStore();
  const isDark = mode === 'dark';
  const [showDetailed, setShowDetailed] = useState(false);

  return (
    <Card 
      title={
        <Space style={{ justifyContent: 'space-between', width: '100%' }}>
          <Space>
            <QuestionCircleOutlined style={{ color: '#1890ff' }} />
            <span style={{ fontSize: '13px' }}>什么是Pareto最优？</span>
          </Space>
          {onClose && (
            <Button 
              type="text" 
              size="small" 
              icon={<CloseOutlined />} 
              onClick={onClose}
              style={{ padding: '0 4px' }}
            />
          )}
        </Space>
      }
      size="small"
      style={{
        border: `1px solid ${isDark ? '#434343' : '#d9d9d9'}`,
        background: isDark ? '#1f1f1f' : '#fafafa',
        boxShadow: onClose ? '0 4px 12px rgba(0,0,0,0.15)' : '0 2px 8px rgba(0,0,0,0.1)'
      }}
      bodyStyle={{ padding: '12px' }}
    >
      <Space direction="vertical" style={{ width: '100%' }} size="small">
        {/* 核心概念 - 紧凑版 */}
        <div style={{ display: 'flex', gap: '16px', alignItems: 'flex-start' }}>
          <div style={{ flex: 1 }}>
            <Paragraph style={{ margin: 0, fontSize: '12px' }}>
              <Text strong style={{ color: '#1890ff' }}>Pareto最优解</Text>
              是指无法在不恶化其他目标的情况下改善任何一个目标的解。
            </Paragraph>
          </div>

          {/* 快速示例 */}
          <div style={{
            background: isDark ? '#262626' : '#f0f8ff',
            padding: '8px',
            borderRadius: '4px',
            border: `1px solid ${isDark ? '#434343' : '#91d5ff'}`,
            flex: 1.5,
            minWidth: '200px'
          }}>
            <Text style={{ fontSize: '11px' }}>
              <Text strong>💡 例子：</Text> 分子A(gap=8.0, fp=150)和分子B(gap=7.5, fp=180)都是Pareto最优的，
              因为A的gap更高，B的fp更高，无法简单判断哪个更好。
            </Text>
          </div>
        </div>

        {/* 图表说明 */}
        <div style={{ textAlign: 'center' }}>
          <Space size="small">
            <div style={{ display: 'flex', alignItems: 'center' }}>
              <div style={{ 
                width: '8px', 
                height: '8px', 
                borderRadius: '50%', 
                background: '#ff4d4f',
                border: '1px solid #ffffff',
                marginRight: '4px' 
              }} />
              <Text style={{ fontSize: '10px' }}>Pareto最优</Text>
            </div>
            <div style={{ display: 'flex', alignItems: 'center' }}>
              <div style={{ 
                width: '6px', 
                height: '6px', 
                borderRadius: '50%', 
                background: '#1890ff',
                opacity: 0.6,
                marginRight: '4px' 
              }} />
              <Text style={{ fontSize: '10px' }}>其他候选</Text>
            </div>
          </Space>
        </div>

        {/* 展开/收起按钮 */}
        <div style={{ textAlign: 'center' }}>
          <Button 
            type="link" 
            size="small"
            onClick={() => setShowDetailed(!showDetailed)}
            icon={<BookOutlined />}
            style={{ fontSize: '11px', padding: '0' }}
          >
            {showDetailed ? '收起详细说明' : '查看详细说明'}
          </Button>
        </div>

        {/* 详细说明 - 可折叠 */}
        {showDetailed && (
          <Collapse ghost size="small">
            <Panel header="📚 详细说明" key="theory" style={{ fontSize: '11px' }}>
              <Space direction="vertical" style={{ width: '100%' }} size="small">
                
                {/* 数学定义 */}
                <div>
                  <Text strong style={{ fontSize: '11px', color: '#52c41a' }}>🔢 数学定义</Text>
                  <div style={{
                    background: isDark ? '#1f1f1f' : '#f6ffed',
                    padding: '6px',
                    borderRadius: '4px',
                    border: `1px solid ${isDark ? '#434343' : '#b7eb8f'}`,
                    fontSize: '10px',
                    fontFamily: 'monospace'
                  }}>
                    对于解x₁和x₂，如果满足：<br/>
                    • ∀i: f_i(x₁) ≥ f_i(x₂) （所有目标都不差）<br/>
                    • ∃j: f_j(x₁) &gt; f_j(x₂) （至少一个目标更好）<br/>
                    则称x₁支配x₂
                  </div>
                </div>

                {/* 优势对比 */}
                <div>
                  <Text strong style={{ fontSize: '11px', color: '#722ed1' }}>✨ 优势对比</Text>
                  <div style={{ display: 'flex', gap: '4px', marginTop: '4px' }}>
                    <Tag color="red" style={{ fontSize: '10px' }}>传统单目标</Tag>
                    <Tag color="orange" style={{ fontSize: '10px' }}>加权求和</Tag>
                    <Tag color="green" style={{ fontSize: '10px' }}>Pareto优化</Tag>
                  </div>
                  <Text style={{ fontSize: '10px', color: '#666' }}>
                    Pareto优化找到所有最优解，提供权衡选择，客观数学标准
                  </Text>
                </div>

                {/* 使用指导 */}
                <div>
                  <Text strong style={{ fontSize: '11px', color: '#13c2c2' }}>🎯 使用指导</Text>
                  <ol style={{ paddingLeft: '16px', margin: '4px 0', fontSize: '10px' }}>
                    <li>分析Pareto前沿 - 观察目标间权衡关系</li>
                    <li>识别关键区域 - 找到性能平衡点</li>
                    <li>结合实际需求 - 根据应用场景选择最适合的解</li>
                  </ol>
                </div>
              </Space>
            </Panel>
          </Collapse>
        )}
      </Space>
    </Card>
  );
};

export default CompactParetoExplanation;
