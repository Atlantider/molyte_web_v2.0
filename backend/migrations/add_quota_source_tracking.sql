-- 添加配额消费来源追踪字段
-- 目的：记录子账号配额消费来自哪个来源（个人池或分配池）
-- 日期：2025-12-18

-- ============ 修改 quota_transactions 表 ============

-- 1. 添加消费来源字段
ALTER TABLE quota_transactions 
ADD COLUMN IF NOT EXISTS source VARCHAR(50) DEFAULT 'personal' NOT NULL;

-- 2. 添加个人池消费记录
ALTER TABLE quota_transactions 
ADD COLUMN IF NOT EXISTS personal_consumed FLOAT DEFAULT 0.0 NOT NULL;

-- 3. 添加分配池消费记录
ALTER TABLE quota_transactions 
ADD COLUMN IF NOT EXISTS allocated_consumed FLOAT DEFAULT 0.0 NOT NULL;

-- ============ 创建索引 ============

-- 为 source 字段创建索引以加快查询
CREATE INDEX IF NOT EXISTS idx_quota_transactions_source ON quota_transactions(source);

-- 为 user_id 和 type 创建复合索引
CREATE INDEX IF NOT EXISTS idx_quota_transactions_user_type ON quota_transactions(user_id, type);

-- ============ 添加注释 ============

COMMENT ON COLUMN quota_transactions.source IS '配额消费来源：personal（个人池）| allocated（分配池）| mixed（混合）';
COMMENT ON COLUMN quota_transactions.personal_consumed IS '从个人充值池消费的核时数';
COMMENT ON COLUMN quota_transactions.allocated_consumed IS '从主账号分配池消费的核时数';

-- ============ 验证迁移 ============

-- 检查新字段是否添加成功
-- SELECT column_name, data_type, is_nullable, column_default 
-- FROM information_schema.columns 
-- WHERE table_name = 'quota_transactions' 
-- AND column_name IN ('source', 'personal_consumed', 'allocated_consumed');

