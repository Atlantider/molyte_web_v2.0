import React from 'react';
import ReactECharts from 'echarts-for-react';
import { Card, Typography, Space, Tag, Tooltip } from 'antd';
import { useThemeStore } from '../stores/themeStore';
import {
  getScientificChartConfig,
  getParetoScatterConfig,
  generateScientificTooltip
} from '../utils/scientificChartTheme';

const { Text } = Typography;

interface ParetoFrontierChartProps {
  data: any[];
  paretoFrontier: any[];
  objectives: string[];
  title?: string;
  height?: number;
}

const ParetoFrontierChart: React.FC<ParetoFrontierChartProps> = ({
  data,
  paretoFrontier,
  objectives,
  title = "Pareto前沿分析",
  height = 400
}) => {
  const { mode } = useThemeStore();
  const isDark = mode === 'dark';

  if (!data || data.length === 0 || objectives.length < 2) {
    return (
      <Card>
        <div style={{ textAlign: 'center', padding: '40px' }}>
          <Text type="secondary">需要至少2个优化目标和数据才能显示Pareto前沿图</Text>
        </div>
      </Card>
    );
  }

  const [obj1, obj2, obj3] = objectives.slice(0, 3);

  // 准备图表数据
  const allPoints = data.map(mol => ({
    value: [
      mol.properties?.[obj1] || 0,
      mol.properties?.[obj2] || 0
    ],
    smiles: mol.smiles,
    similarity: mol.similarity || 0,
    isParetoOptimal: paretoFrontier.some(p => p.smiles === mol.smiles),
    properties: mol.properties || {},
    obj3Value: obj3 ? (mol.properties?.[obj3] || 0) : null
  }));

  const paretoPoints = allPoints.filter(p => p.isParetoOptimal);
  const regularPoints = allPoints.filter(p => !p.isParetoOptimal);

  // 计算数据范围用于坐标轴
  const obj1Values = allPoints.map(p => p.value[0]);
  const obj2Values = allPoints.map(p => p.value[1]);
  const obj1Range = [Math.min(...obj1Values), Math.max(...obj1Values)];
  const obj2Range = [Math.min(...obj2Values), Math.max(...obj2Values)];

  // 使用科学期刊主题配置
  const chartConfig = getScientificChartConfig(isDark);
  const scatterConfig = getParetoScatterConfig(isDark);

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
        const data = params.data;
        const tooltipContent = [
          `<div style="padding: 8px; font-family: Arial, sans-serif;">`,
          `<div style="font-weight: bold; margin-bottom: 8px; color: ${data.isParetoOptimal ? '#ff4d4f' : '#1890ff'};">`,
          `${data.isParetoOptimal ? '⭐ Pareto最优解' : '🔵 候选分子'}`,
          `</div>`,
          `<div style="margin-bottom: 6px;">`,
          `<strong>${obj1.toUpperCase()}:</strong> ${data.value[0].toFixed(3)}`,
          `</div>`,
          `<div style="margin-bottom: 6px;">`,
          `<strong>${obj2.toUpperCase()}:</strong> ${data.value[1].toFixed(3)}`,
          `</div>`,
          obj3 && data.obj3Value !== null ? [
            `<div style="margin-bottom: 6px;">`,
            `<strong>${obj3.toUpperCase()}:</strong> ${data.obj3Value.toFixed(3)}`,
            `</div>`
          ].join('') : '',
          `<div style="margin-bottom: 6px;">`,
          `<strong>相似度:</strong> ${(data.similarity * 100).toFixed(1)}%`,
          `</div>`,
          `<div style="font-size: 11px; color: #666; margin-top: 8px;">`,
          `SMILES: ${data.smiles.length > 30 ? data.smiles.substring(0, 30) + '...' : data.smiles}`,
          `</div>`,
          `</div>`
        ].join('');
        return tooltipContent;
      }
    },
    legend: {
      ...chartConfig.legend,
      data: ['Pareto Optimal', 'Candidates']
    },
    grid: chartConfig.grid,
    xAxis: {
      type: 'value',
      name: obj1.toUpperCase(),
      nameLocation: 'middle',
      nameGap: 35,
      ...chartConfig.axisCommon,
      axisLabel: {
        ...chartConfig.axisCommon.axisLabel,
        formatter: (value: number) => value.toFixed(2)
      },
      min: (value: any) => Math.max(0, value.min - (value.max - value.min) * 0.05),
      max: (value: any) => value.max + (value.max - value.min) * 0.05
    },
    yAxis: {
      type: 'value',
      name: obj2.toUpperCase(),
      nameLocation: 'middle',
      nameGap: 55,
      ...chartConfig.axisCommon,
      axisLabel: {
        ...chartConfig.axisCommon.axisLabel,
        formatter: (value: number) => value.toFixed(2)
      },
      min: (value: any) => Math.max(0, value.min - (value.max - value.min) * 0.05),
      max: (value: any) => value.max + (value.max - value.min) * 0.05
    },
    series: [
      {
        ...scatterConfig.paretoSeries,
        data: paretoPoints
      },
      {
        ...scatterConfig.candidateSeries,
        data: regularPoints
      }
    ]
  };

  return (
    <Card>
      <Space direction="vertical" style={{ width: '100%' }} size="middle">
        {/* 图表标题和统计信息 */}
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
          <Space>
            <Tag color="red">⭐ Pareto最优: {paretoPoints.length}</Tag>
            <Tag color="blue">🔵 其他候选: {regularPoints.length}</Tag>
          </Space>
          <Space>
            <Text type="secondary" style={{ fontSize: '12px' }}>
              {obj1.toUpperCase()} 范围: {obj1Range[0].toFixed(2)} - {obj1Range[1].toFixed(2)}
            </Text>
            <Text type="secondary" style={{ fontSize: '12px' }}>
              {obj2.toUpperCase()} 范围: {obj2Range[0].toFixed(2)} - {obj2Range[1].toFixed(2)}
            </Text>
            {obj3 && (
              <Text type="secondary" style={{ fontSize: '12px' }}>
                {obj3.toUpperCase()} 范围: {Math.min(...allPoints.map(p => p.obj3Value || 0)).toFixed(2)} - {Math.max(...allPoints.map(p => p.obj3Value || 0)).toFixed(2)}
              </Text>
            )}
          </Space>
        </div>

        {/* ECharts图表 */}
        <ReactECharts
          option={option}
          style={{ height: `${height}px`, width: '100%' }}
          theme={isDark ? 'dark' : 'light'}
        />

        {/* 图表说明 */}
        <div style={{ 
          background: isDark ? '#1f1f1f' : '#fafafa',
          padding: '12px',
          borderRadius: '6px',
          border: `1px solid ${isDark ? '#434343' : '#d9d9d9'}`
        }}>
          <Text style={{ fontSize: '12px', color: isDark ? '#ffffff' : '#666666' }}>
            💡 <strong>图表说明：</strong>
            红色星形点表示Pareto最优解（在所有目标上都不被其他解支配），
            蓝色圆点表示其他候选分子。鼠标悬停可查看详细信息。
          </Text>
        </div>
      </Space>
    </Card>
  );
};

export default ParetoFrontierChart;
