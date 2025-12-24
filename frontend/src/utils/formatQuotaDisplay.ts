/**
 * 配额显示格式化工具
 * 统一所有核时、配额等数值的显示格式
 * 规则：最多两位小数
 */

/**
 * 格式化核时显示 - 最多两位小数
 * @param value 核时数值
 * @param defaultValue 默认值（当 value 为 null/undefined 时）
 * @returns 格式化后的字符串
 */
export const formatCpuHours = (value: number | null | undefined, defaultValue: string = '0.00'): string => {
  if (value === null || value === undefined) {
    return defaultValue;
  }
  return value.toFixed(2);
};

/**
 * 格式化配额显示 - 最多两位小数
 * @param value 配额数值
 * @param defaultValue 默认值
 * @returns 格式化后的字符串
 */
export const formatQuota = (value: number | null | undefined, defaultValue: string = '0.00'): string => {
  if (value === null || value === undefined) {
    return defaultValue;
  }
  return value.toFixed(2);
};

/**
 * 格式化余额显示 - 最多两位小数
 * @param value 余额数值
 * @param defaultValue 默认值
 * @returns 格式化后的字符串
 */
export const formatBalance = (value: number | null | undefined, defaultValue: string = '0.00'): string => {
  if (value === null || value === undefined) {
    return defaultValue;
  }
  return value.toFixed(2);
};

/**
 * 格式化价格显示 - 最多两位小数
 * @param value 价格数值
 * @param defaultValue 默认值
 * @returns 格式化后的字符串
 */
export const formatPrice = (value: number | null | undefined, defaultValue: string = '0.00'): string => {
  if (value === null || value === undefined) {
    return defaultValue;
  }
  return value.toFixed(2);
};

/**
 * 格式化百分比显示 - 最多两位小数
 * @param value 百分比数值（0-100）
 * @param defaultValue 默认值
 * @returns 格式化后的字符串
 */
export const formatPercentage = (value: number | null | undefined, defaultValue: string = '0.00'): string => {
  if (value === null || value === undefined) {
    return defaultValue;
  }
  return value.toFixed(2);
};

/**
 * 获取 Ant Design Statistic 组件的 precision 值
 * 统一使用 2 位小数
 */
export const QUOTA_PRECISION = 2;

/**
 * 格式化数值为指定小数位数
 * @param value 数值
 * @param decimals 小数位数（默认 2）
 * @returns 格式化后的字符串
 */
export const formatNumber = (value: number | null | undefined, decimals: number = 2): string => {
  if (value === null || value === undefined) {
    return '0.' + '0'.repeat(decimals);
  }
  return value.toFixed(decimals);
};

/**
 * 格式化核时显示（带单位）
 * @param value 核时数值
 * @param unit 单位（默认 '核时'）
 * @returns 格式化后的字符串
 */
export const formatCpuHoursWithUnit = (value: number | null | undefined, unit: string = '核时'): string => {
  return `${formatCpuHours(value)} ${unit}`;
};

/**
 * 安全的 toFixed 包装函数 - 防止 undefined/null 导致的错误
 * @param value 数值，可以是 number | null | undefined
 * @param digits 小数位数（默认 2）
 * @returns 格式化后的字符串
 */
export const safeToFixed = (value: number | null | undefined, digits: number = 2): string => {
  if (value === null || value === undefined) {
    return '0.' + '0'.repeat(digits);
  }
  return value.toFixed(digits);
};

/**
 * 格式化配额显示（带单位）
 * @param value 配额数值
 * @param unit 单位（默认 'h'）
 * @returns 格式化后的字符串
 */
export const formatQuotaWithUnit = (value: number | null | undefined, unit: string = 'h'): string => {
  return `${formatQuota(value)} ${unit}`;
};

/**
 * 格式化价格显示（带单位）
 * @param value 价格数值
 * @param unit 单位（默认 '¥'）
 * @returns 格式化后的字符串
 */
export const formatPriceWithUnit = (value: number | null | undefined, unit: string = '¥'): string => {
  return `${unit}${formatPrice(value)}/核时`;
};

