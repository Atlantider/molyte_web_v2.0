-- 修复配额系统迁移脚本
-- 目标：
-- 1. 删除 MasterAccount 表中重复的配额字段
-- 2. 添加 SubAccount.allocated_quota 字段
-- 3. 删除 SubAccount 中不需要的字段
-- 日期：2025-12-18

-- ============ 第一步：SubAccount 表修改 ============

-- 1. 添加 allocated_quota 字段（如果不存在）
ALTER TABLE sub_accounts ADD COLUMN IF NOT EXISTS allocated_quota FLOAT DEFAULT 0.0 NOT NULL;

-- 2. 删除不需要的字段（如果存在）
-- 注意：personal_used 和 allocated_used 不再使用，因为配额消费直接从 User 表的 balance_cpu_hours 扣除
ALTER TABLE sub_accounts DROP COLUMN IF EXISTS personal_used;
ALTER TABLE sub_accounts DROP COLUMN IF EXISTS allocated_used;

-- 3. 删除 personal_quota 字段（不再使用）
ALTER TABLE sub_accounts DROP COLUMN IF EXISTS personal_quota;

-- 添加注释
COMMENT ON COLUMN sub_accounts.allocated_quota IS '主账号分配给子账号的配额（子账号实际可用 = min(User.balance_cpu_hours, allocated_quota)）';

-- ============ 第二步：MasterAccount 表修改 ============

-- 1. 删除重复的配额字段（如果存在）
-- 主账号的配额应该来自 User 表，不在 MasterAccount 表中重复存储
ALTER TABLE master_accounts DROP COLUMN IF EXISTS total_cpu_hours;
ALTER TABLE master_accounts DROP COLUMN IF EXISTS balance_cpu_hours;
ALTER TABLE master_accounts DROP COLUMN IF EXISTS frozen_cpu_hours;
ALTER TABLE master_accounts DROP COLUMN IF EXISTS used_cpu_hours;

-- ============ 第三步：创建索引 ============

-- 为 SubAccount 添加索引以加快查询
CREATE INDEX IF NOT EXISTS idx_sub_accounts_allocated_quota ON sub_accounts(allocated_quota);

-- ============ 第四步：数据迁移 ============

-- 如果 SubAccount 中有 personal_quota 数据，迁移到 allocated_quota
-- （这个操作在删除 personal_quota 之前执行）
-- UPDATE sub_accounts SET allocated_quota = COALESCE(personal_quota, 0) WHERE allocated_quota = 0;

-- ============ 第五步：验证 ============

-- 验证 SubAccount 表结构
-- SELECT column_name, data_type, is_nullable FROM information_schema.columns 
-- WHERE table_name = 'sub_accounts' ORDER BY ordinal_position;

-- 验证 MasterAccount 表结构
-- SELECT column_name, data_type, is_nullable FROM information_schema.columns 
-- WHERE table_name = 'master_accounts' ORDER BY ordinal_position;

