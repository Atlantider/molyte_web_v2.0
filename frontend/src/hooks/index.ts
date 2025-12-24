/**
 * Hooks 导出文件
 */

export { useAccountType, getAccountTypeLabel, getAccountTypeColor, hasPermission } from './useAccountType';
export type { AccountInfo, QuotaInfo } from './useAccountType';

export { useQuota, formatQuota, getQuotaStatusColor, getQuotaStatusText, getQuotaSourceLabel } from './useQuota';
export type { QuotaData, QuotaSource } from './useQuota';

