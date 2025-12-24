-- 统一核时系统迁移
-- 目标：将 balance_cpu_hours 作为唯一的真实来源，删除 total_cpu_hours 和 debt_cpu_hours
-- 时间：2025-12-18

-- 1. 添加缺失的字段（如果不存在）
ALTER TABLE users ADD COLUMN IF NOT EXISTS recharge_cpu_hours FLOAT DEFAULT 0.0 NOT NULL;

-- 2. 更新 recharge_cpu_hours：从 QuotaTransaction 表统计充值记录
UPDATE users u SET recharge_cpu_hours = (
    SELECT COALESCE(SUM(amount), 0) FROM quota_transactions 
    WHERE user_id = u.id AND type = 'recharge'
) WHERE recharge_cpu_hours = 0;

-- 3. 删除 total_cpu_hours 和 debt_cpu_hours 字段（可选，保留兼容性）
-- ALTER TABLE users DROP COLUMN IF EXISTS total_cpu_hours;
-- ALTER TABLE users DROP COLUMN IF EXISTS debt_cpu_hours;

-- 4. 添加注释说明新的核时系统
COMMENT ON COLUMN users.balance_cpu_hours IS '可用余额（核时）- 正数表示可用，负数表示欠费，0表示无可用核时。这是唯一的真实来源。';
COMMENT ON COLUMN users.frozen_cpu_hours IS '冻结核时 - 运行中的任务占用的核时';
COMMENT ON COLUMN users.free_cpu_hours_granted IS '初始赠送的免费核时 - 用于统计核时来源';
COMMENT ON COLUMN users.recharge_cpu_hours IS '充值获得的核时 - 用于统计核时来源';
COMMENT ON COLUMN users.admin_granted_cpu_hours IS '管理员赠送的核时 - 用于统计核时来源';

-- 5. 创建视图用于查询核时统计
CREATE OR REPLACE VIEW user_cpu_hours_summary AS
SELECT 
    u.id,
    u.username,
    u.balance_cpu_hours,
    u.frozen_cpu_hours,
    u.free_cpu_hours_granted,
    u.recharge_cpu_hours,
    u.admin_granted_cpu_hours,
    (u.free_cpu_hours_granted + u.recharge_cpu_hours + u.admin_granted_cpu_hours) AS total_granted,
    COALESCE((SELECT SUM(amount) FROM quota_transactions WHERE user_id = u.id AND type = 'consume'), 0) AS total_consumed,
    CASE 
        WHEN u.balance_cpu_hours < 0 THEN '欠费'
        WHEN u.balance_cpu_hours = 0 THEN '无可用'
        ELSE '可用'
    END AS status
FROM users u;

-- 6. 创建索引以提高查询性能
CREATE INDEX IF NOT EXISTS idx_users_balance_cpu_hours ON users(balance_cpu_hours);
CREATE INDEX IF NOT EXISTS idx_users_frozen_cpu_hours ON users(frozen_cpu_hours);

