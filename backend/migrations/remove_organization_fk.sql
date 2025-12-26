-- 移除 master_accounts 表中的 organization_id 外键
-- 原因：organizations 表不存在，导致删除主账号时出现外键约束错误

-- 1. 删除依赖的视图
DROP VIEW IF EXISTS master_account_stats CASCADE;

-- 2. 删除外键约束
ALTER TABLE master_accounts DROP CONSTRAINT IF EXISTS master_accounts_organization_id_fkey;

-- 3. 删除 organization_id 列
ALTER TABLE master_accounts DROP COLUMN IF EXISTS organization_id;

-- 完成
SELECT 'Foreign key and dependent objects removed successfully!' as message;

