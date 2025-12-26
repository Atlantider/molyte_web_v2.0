/**
 * 科学期刊风格图表主题配置
 * 参考Nature、Science、RDF等期刊的绘图标准
 */

// Nature期刊推荐的颜色方案
export const NATURE_COLORS = {
  // 主色调 - 适合色盲友好
  primary: '#1f77b4',      // 蓝色
  secondary: '#ff7f0e',    // 橙色
  success: '#2ca02c',      // 绿色
  danger: '#d62728',       // 红色
  warning: '#9467bd',      // 紫色
  info: '#8c564b',         // 棕色
  
  // Pareto前沿专用色彩
  paretoOptimal: '#e74c3c',     // 鲜红色 - 突出最优解
  candidates: '#3498db',        // 蓝色 - 候选解
  trendLine: '#27ae60',         // 绿色 - 趋势线
  
  // 渐变色方案
  gradients: {
    pareto: ['#e74c3c', '#c0392b'],
    candidate: ['#3498db', '#2980b9'],
    trend: ['#27ae60', '#229954']
  }
};

// RDF期刊风格的字体配置
export const SCIENTIFIC_FONTS = {
  title: {
    family: 'Arial, Helvetica, sans-serif',
    size: 18,
    weight: '600'
  },
  axis: {
    family: 'Arial, Helvetica, sans-serif',
    size: 14,
    weight: '600'
  },
  label: {
    family: 'Arial, Helvetica, sans-serif',
    size: 12,
    weight: '400'
  },
  legend: {
    family: 'Arial, Helvetica, sans-serif',
    size: 13,
    weight: '500'
  },
  tooltip: {
    family: 'Arial, Helvetica, sans-serif',
    size: 13,
    weight: '400'
  }
};

// 科学期刊标准的图表配置
export const getScientificChartConfig = (isDark: boolean = false) => ({
  // 背景和边框
  backgroundColor: 'transparent',
  
  // 网格配置 - 参考Nature期刊
  grid: {
    left: '12%',
    right: '8%',
    bottom: '18%',
    top: '28%',
    containLabel: true,
    backgroundColor: 'transparent',
    borderWidth: 0
  },
  
  // 标题样式
  title: {
    textStyle: {
      color: isDark ? '#ecf0f1' : '#2c3e50',
      fontSize: SCIENTIFIC_FONTS.title.size,
      fontWeight: SCIENTIFIC_FONTS.title.weight,
      fontFamily: SCIENTIFIC_FONTS.title.family
    },
    left: 'center',
    top: '2%'
  },

  // 图例样式
  legend: {
    textStyle: {
      color: isDark ? '#ecf0f1' : '#2c3e50',
      fontSize: SCIENTIFIC_FONTS.legend.size,
      fontFamily: SCIENTIFIC_FONTS.legend.family,
      fontWeight: SCIENTIFIC_FONTS.legend.weight
    },
    top: '10%',
    left: 'center',
    itemGap: 30,
    itemWidth: 16,
    itemHeight: 16
  },
  
  // 坐标轴通用配置
  axisCommon: {
    nameTextStyle: {
      color: isDark ? '#ecf0f1' : '#2c3e50',
      fontSize: SCIENTIFIC_FONTS.axis.size,
      fontWeight: SCIENTIFIC_FONTS.axis.weight,
      fontFamily: SCIENTIFIC_FONTS.axis.family
    },
    axisLabel: {
      color: isDark ? '#bdc3c7' : '#34495e',
      fontSize: SCIENTIFIC_FONTS.label.size,
      fontFamily: SCIENTIFIC_FONTS.label.family,
      margin: 12
    },
    axisLine: {
      show: true,
      lineStyle: {
        color: isDark ? '#5a6c7d' : '#7f8c8d',
        width: 1.5
      }
    },
    axisTick: {
      show: true,
      length: 6,
      lineStyle: {
        color: isDark ? '#5a6c7d' : '#7f8c8d',
        width: 1
      }
    },
    splitLine: {
      show: true,
      lineStyle: {
        color: isDark ? '#3a4a5c' : '#ecf0f1',
        width: 1,
        type: 'solid',
        opacity: 0.6
      }
    }
  },
  
  // 提示框样式
  tooltip: {
    backgroundColor: isDark ? 'rgba(31, 31, 31, 0.95)' : 'rgba(255, 255, 255, 0.95)',
    borderColor: isDark ? '#434343' : '#d9d9d9',
    borderWidth: 1,
    borderRadius: 8,
    textStyle: {
      color: isDark ? '#ffffff' : '#2c3e50',
      fontSize: SCIENTIFIC_FONTS.tooltip.size,
      fontFamily: SCIENTIFIC_FONTS.tooltip.family
    },
    extraCssText: 'box-shadow: 0 4px 12px rgba(0,0,0,0.15);'
  }
});

// Pareto前沿散点图专用配置
export const getParetoScatterConfig = (isDark: boolean = false) => ({
  paretoSeries: {
    name: 'Pareto Optimal',
    type: 'scatter',
    symbolSize: 14,
    symbol: 'diamond',
    itemStyle: {
      color: NATURE_COLORS.paretoOptimal,
      borderColor: '#ffffff',
      borderWidth: 2.5,
      shadowBlur: 8,
      shadowColor: 'rgba(231, 76, 60, 0.3)',
      shadowOffsetX: 0,
      shadowOffsetY: 2
    },
    emphasis: {
      itemStyle: {
        color: '#c0392b',
        borderWidth: 3,
        shadowBlur: 12,
        shadowColor: 'rgba(231, 76, 60, 0.5)',
        scale: 1.2
      }
    },
    zlevel: 3
  },
  
  candidateSeries: {
    name: 'Candidates',
    type: 'scatter',
    symbolSize: 10,
    symbol: 'circle',
    itemStyle: {
      color: NATURE_COLORS.candidates,
      borderColor: isDark ? '#2c3e50' : '#ffffff',
      borderWidth: 1.5,
      opacity: 0.75,
      shadowBlur: 4,
      shadowColor: 'rgba(52, 152, 219, 0.2)',
      shadowOffsetX: 0,
      shadowOffsetY: 1
    },
    emphasis: {
      itemStyle: {
        color: '#2980b9',
        opacity: 0.9,
        borderWidth: 2,
        shadowBlur: 8,
        shadowColor: 'rgba(52, 152, 219, 0.4)',
        scale: 1.15
      }
    },
    zlevel: 2
  }
});

// 权衡分析图专用配置
export const getTradeOffConfig = (isDark: boolean = false) => ({
  scatterSeries: {
    name: 'Data Points',
    type: 'scatter',
    itemStyle: {
      shadowBlur: 4,
      shadowOffsetX: 0,
      shadowOffsetY: 1
    },
    emphasis: {
      itemStyle: {
        shadowBlur: 8,
        scale: 1.15
      }
    }
  },
  
  trendLineSeries: {
    name: 'Trend Line',
    type: 'line',
    lineStyle: {
      color: NATURE_COLORS.trendLine,
      width: 2.5,
      type: 'dashed',
      dashOffset: 5,
      cap: 'round'
    },
    symbol: 'none',
    emphasis: {
      lineStyle: {
        width: 3.5,
        shadowBlur: 8,
        shadowColor: 'rgba(39, 174, 96, 0.3)'
      }
    },
    zlevel: 1
  }
});

// 生成科学期刊风格的提示框HTML
export const generateScientificTooltip = (
  data: any, 
  obj1: string, 
  obj2: string, 
  isDark: boolean = false,
  isPareto: boolean = false
) => {
  const statusColor = isPareto ? NATURE_COLORS.paretoOptimal : NATURE_COLORS.candidates;
  const statusIcon = isPareto ? '★' : '●';
  const statusText = isPareto ? 'Pareto Optimal' : 'Candidate';
  
  return `
    <div style="padding: 12px; font-family: ${SCIENTIFIC_FONTS.tooltip.family};">
      <div style="font-weight: 600; margin-bottom: 8px; color: ${statusColor}; font-size: 14px;">
        ${statusIcon} ${statusText}
      </div>
      <div style="font-family: 'Courier New', monospace; font-size: 11px; margin-bottom: 10px; 
                  background: ${isDark ? '#2c2c2c' : '#f8f9fa'}; padding: 6px; border-radius: 4px;">
        ${data.smiles.substring(0, 35)}${data.smiles.length > 35 ? '...' : ''}
      </div>
      <div style="margin-bottom: 6px; font-size: 12px;">
        <span style="font-weight: 600; color: #e67e22;">${obj1.toUpperCase()}:</span> 
        <span style="font-weight: 500;">${data.value[0].toFixed(3)}</span>
      </div>
      <div style="margin-bottom: 6px; font-size: 12px;">
        <span style="font-weight: 600; color: #9b59b6;">${obj2.toUpperCase()}:</span> 
        <span style="font-weight: 500;">${data.value[1].toFixed(3)}</span>
      </div>
      <div style="font-size: 12px; color: #7f8c8d;">
        <span style="font-weight: 600;">Similarity:</span> ${(data.similarity * 100).toFixed(1)}%
      </div>
    </div>
  `;
};

// 生成趋势线提示框
export const generateTrendLineTooltip = (
  slope: number, 
  intercept: number, 
  correlation: number,
  isDark: boolean = false
) => {
  return `
    <div style="padding: 10px; font-family: ${SCIENTIFIC_FONTS.tooltip.family};">
      <div style="font-weight: 600; margin-bottom: 6px; color: ${NATURE_COLORS.trendLine};">
        📈 Linear Regression
      </div>
      <div style="font-family: 'Courier New', monospace; font-size: 12px;">
        y = ${slope.toFixed(3)}x + ${intercept.toFixed(3)}
      </div>
      <div style="font-size: 11px; color: #7f8c8d; margin-top: 4px;">
        R = ${correlation.toFixed(3)}
      </div>
    </div>
  `;
};
