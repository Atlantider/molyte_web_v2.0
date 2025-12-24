import React from 'react';
import ReactECharts from 'echarts-for-react';
import { Card, Typography, Space, Tag, Row, Col, Statistic } from 'antd';
import { useThemeStore } from '../stores/themeStore';
import {
  getScientificChartConfig,
  getTradeOffConfig,
  generateScientificTooltip,
  generateTrendLineTooltip
} from '../utils/scientificChartTheme';

const { Text } = Typography;

interface TradeOffAnalysisChartProps {
  data: any[];
  objectives: string[];
  title?: string;
  height?: number;
}

const TradeOffAnalysisChart: React.FC<TradeOffAnalysisChartProps> = ({
  data,
  objectives,
  title = "目标权衡分析",
  height = 350
}) => {
  const { mode } = useThemeStore();
  const isDark = mode === 'dark';

  if (!data || data.length === 0 || objectives.length < 2) {
    return (
      <Card>
        <div style={{ textAlign: 'center', padding: '40px' }}>
          <Text type="secondary">需要至少2个目标和数据才能显示权衡分析</Text>
        </div>
      </Card>
    );
  }

  const [obj1, obj2] = objectives.slice(0, 2);

  // 计算相关性
  const values1 = data.map(mol => mol.properties?.[obj1] || 0);
  const values2 = data.map(mol => mol.properties?.[obj2] || 0);
  
  const mean1 = values1.reduce((a, b) => a + b, 0) / values1.length;
  const mean2 = values2.reduce((a, b) => a + b, 0) / values2.length;
  
  const numerator = values1.reduce((sum, val1, i) => sum + (val1 - mean1) * (values2[i] - mean2), 0);
  const denominator = Math.sqrt(
    values1.reduce((sum, val) => sum + Math.pow(val - mean1, 2), 0) *
    values2.reduce((sum, val) => sum + Math.pow(val - mean2, 2), 0)
  );
  
  const correlation = denominator > 0 ? numerator / denominator : 0;

  // 准备散点图数据
  const scatterData = data.map((mol, index) => ({
    value: [values1[index], values2[index]],
    smiles: mol.smiles,
    similarity: mol.similarity || 0,
    isParetoOptimal: mol.pareto_rank === 1 || mol.paretoRank === 1
  }));

  // 计算趋势线（简单线性回归）
  const slope = correlation * (Math.sqrt(values2.reduce((sum, val) => sum + Math.pow(val - mean2, 2), 0) / values2.length) / 
                               Math.sqrt(values1.reduce((sum, val) => sum + Math.pow(val - mean1, 2), 0) / values1.length));
  const intercept = mean2 - slope * mean1;

  const minX = Math.min(...values1);
  const maxX = Math.max(...values1);
  const trendLineData = [
    [minX, slope * minX + intercept],
    [maxX, slope * maxX + intercept]
  ];

  // 计算坐标轴范围
  const xPadding = (maxX - minX) * 0.1; // 10% 边距
  const yMin = Math.min(...values2);
  const yMax = Math.max(...values2);
  const yPadding = (yMax - yMin) * 0.1; // 10% 边距

  // 使用科学期刊主题配置
  const chartConfig = getScientificChartConfig(isDark);
  const tradeOffConfig = getTradeOffConfig(isDark);

  const option = {
    ...chartConfig,
    title: {
      ...chartConfig.title,
      text: title
    },
    tooltip: {
      ...chartConfig.tooltip,
      trigger: 'item',
      formatter: (params: any) => {
        if (params.seriesName === 'Trend Line') {
          return generateTrendLineTooltip(slope, intercept, correlation, isDark);
        }
        const data = params.data;
        return generateScientificTooltip(data, obj1, obj2, isDark, data.isParetoOptimal);
      }
    },
    legend: {
      ...chartConfig.legend,
      data: ['Data Points', 'Trend Line']
    },
    grid: chartConfig.grid,
    xAxis: {
      type: 'value',
      name: obj1.toUpperCase(),
      nameLocation: 'middle',
      nameGap: 35,
      min: minX - xPadding,
      max: maxX + xPadding,
      ...chartConfig.axisCommon,
      axisLabel: {
        ...chartConfig.axisCommon.axisLabel,
        formatter: (value: number) => value.toFixed(2)
      }
    },
    yAxis: {
      type: 'value',
      name: obj2.toUpperCase(),
      nameLocation: 'middle',
      nameGap: 55,
      min: yMin - yPadding,
      max: yMax + yPadding,
      ...chartConfig.axisCommon,
      axisLabel: {
        ...chartConfig.axisCommon.axisLabel,
        formatter: (value: number) => value.toFixed(2)
      }
    },
    series: [
      {
        ...tradeOffConfig.scatterSeries,
        data: scatterData,
        symbolSize: (data: any) => data.isParetoOptimal ? 14 : 10,
        symbol: (data: any) => data.isParetoOptimal ? 'diamond' : 'circle',
        itemStyle: {
          ...tradeOffConfig.scatterSeries.itemStyle,
          color: (params: any) => params.data.isParetoOptimal ? '#e74c3c' : '#3498db',
          borderColor: (params: any) => params.data.isParetoOptimal ? '#ffffff' : (isDark ? '#2c3e50' : '#ffffff'),
          borderWidth: (params: any) => params.data.isParetoOptimal ? 2.5 : 1.5,
          shadowColor: (params: any) => params.data.isParetoOptimal ? 'rgba(231, 76, 60, 0.3)' : 'rgba(52, 152, 219, 0.2)'
        },
        emphasis: {
          ...tradeOffConfig.scatterSeries.emphasis,
          itemStyle: {
            ...tradeOffConfig.scatterSeries.emphasis.itemStyle,
            borderWidth: (params: any) => params.data.isParetoOptimal ? 3 : 2,
            shadowColor: (params: any) => params.data.isParetoOptimal ? 'rgba(231, 76, 60, 0.5)' : 'rgba(52, 152, 219, 0.4)'
          }
        }
      },
      {
        ...tradeOffConfig.trendLineSeries,
        data: trendLineData
      }
    ]
  };

  // 相关性分析结果
  const getCorrelationDescription = (corr: number) => {
    if (corr > 0.7) return { text: '强正相关', color: '#52c41a', icon: '📈' };
    if (corr > 0.3) return { text: '中等正相关', color: '#1890ff', icon: '📊' };
    if (corr > -0.3) return { text: '弱相关', color: '#faad14', icon: '➡️' };
    if (corr > -0.7) return { text: '中等负相关', color: '#fa8c16', icon: '📉' };
    return { text: '强负相关', color: '#ff4d4f', icon: '🔄' };
  };

  const corrDesc = getCorrelationDescription(correlation);

  return (
    <Card>
      <Space direction="vertical" style={{ width: '100%' }} size="middle">
        {/* 相关性统计 */}
        <Row gutter={16}>
          <Col span={8}>
            <Statistic
              title="相关系数"
              value={correlation}
              precision={3}
              valueStyle={{ color: corrDesc.color }}
              prefix={corrDesc.icon}
            />
          </Col>
          <Col span={8}>
            <Statistic
              title="关系类型"
              value={corrDesc.text}
              valueStyle={{ color: corrDesc.color, fontSize: '16px' }}
            />
          </Col>
          <Col span={8}>
            <Statistic
              title="数据点数"
              value={data.length}
              prefix="📊"
            />
          </Col>
        </Row>

        {/* 权衡分析图表 */}
        <ReactECharts
          option={option}
          style={{ height: `${height}px`, width: '100%' }}
          theme={isDark ? 'dark' : 'light'}
        />

        {/* 分析说明 */}
        <div style={{ 
          background: isDark ? '#1f1f1f' : '#fafafa',
          padding: '12px',
          borderRadius: '6px',
          border: `1px solid ${isDark ? '#434343' : '#d9d9d9'}`
        }}>
          <Text style={{ fontSize: '12px', color: isDark ? '#ffffff' : '#666666' }}>
            💡 <strong>权衡分析说明：</strong>
            {correlation > 0.3 && '两个目标呈正相关，可以协同优化。'}
            {correlation < -0.3 && '两个目标呈负相关，存在权衡关系，提高一个目标可能会降低另一个。'}
            {Math.abs(correlation) <= 0.3 && '两个目标相对独立，可以分别优化。'}
            绿色虚线表示数据的整体趋势。
          </Text>
        </div>
      </Space>
    </Card>
  );
};

export default TradeOffAnalysisChart;
