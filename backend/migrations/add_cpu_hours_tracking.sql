-- 添加核时来源追踪字段
-- 用于区分免费核时、充值核时、管理员赠送核时

-- 1. 添加新字段到 users 表
ALTER TABLE users
ADD COLUMN IF NOT EXISTS recharge_cpu_hours FLOAT DEFAULT 0.0,
ADD COLUMN IF NOT EXISTS admin_granted_cpu_hours FLOAT DEFAULT 0.0;

-- 2. 初始化现有用户的数据
-- 将 free_cpu_hours_granted 作为初始值
UPDATE users 
SET recharge_cpu_hours = 0.0, 
    admin_granted_cpu_hours = 0.0
WHERE recharge_cpu_hours IS NULL OR admin_granted_cpu_hours IS NULL;

-- 3. 添加注释
COMMENT ON COLUMN users.free_cpu_hours_granted IS '初始赠送的免费核时';
COMMENT ON COLUMN users.recharge_cpu_hours IS '用户充值的核时';
COMMENT ON COLUMN users.admin_granted_cpu_hours IS '管理员赠送的核时';
COMMENT ON COLUMN users.balance_cpu_hours IS '可用余额 = free_cpu_hours_granted + recharge_cpu_hours + admin_granted_cpu_hours - 已消费';

-- 4. 创建视图用于查询核时来源分解
CREATE OR REPLACE VIEW user_cpu_hours_breakdown AS
SELECT 
    u.id,
    u.username,
    u.free_cpu_hours_granted,
    u.recharge_cpu_hours,
    u.admin_granted_cpu_hours,
    (u.free_cpu_hours_granted + u.recharge_cpu_hours + u.admin_granted_cpu_hours) AS total_granted,
    u.balance_cpu_hours,
    u.frozen_cpu_hours,
    u.debt_cpu_hours,
    u.total_cpu_hours AS quota_limit
FROM users u;

-- 5. 添加索引以提高查询性能
CREATE INDEX IF NOT EXISTS idx_users_balance_cpu_hours ON users(balance_cpu_hours);
CREATE INDEX IF NOT EXISTS idx_users_admin_granted_cpu_hours ON users(admin_granted_cpu_hours);

