-- 简化账户体系迁移脚本
-- 目标：去掉组织管理，只保留主账号和子账号的简单关系
-- 日期：2025-12-18

-- 1. 为 sub_accounts 表添加分配配额字段
ALTER TABLE sub_accounts ADD COLUMN IF NOT EXISTS allocated_quota FLOAT DEFAULT 0.0 NOT NULL;
ALTER TABLE sub_accounts ADD COLUMN IF NOT EXISTS allocated_used FLOAT DEFAULT 0.0 NOT NULL;

-- 2. 添加注释说明
COMMENT ON COLUMN sub_accounts.personal_quota IS '子账号的个人配额，由子账号自己管理';
COMMENT ON COLUMN sub_accounts.personal_used IS '子账号个人配额已使用部分';
COMMENT ON COLUMN sub_accounts.allocated_quota IS '主账号分配给子账号的配额';
COMMENT ON COLUMN sub_accounts.allocated_used IS '子账号已使用的分配配额';

-- 3. 创建索引以加快查询
CREATE INDEX IF NOT EXISTS idx_sub_accounts_master_account_id ON sub_accounts(master_account_id);
CREATE INDEX IF NOT EXISTS idx_sub_accounts_user_id ON sub_accounts(user_id);
CREATE INDEX IF NOT EXISTS idx_master_accounts_user_id ON master_accounts(user_id);

-- 4. 为 master_accounts 表添加个人配额字段（用于追踪主账号自己的配额使用）
ALTER TABLE master_accounts ADD COLUMN IF NOT EXISTS personal_quota FLOAT DEFAULT 0.0 NOT NULL;
ALTER TABLE master_accounts ADD COLUMN IF NOT EXISTS personal_used FLOAT DEFAULT 0.0 NOT NULL;

COMMENT ON COLUMN master_accounts.personal_quota IS '主账号自己的个人配额（从充值中扣除）';
COMMENT ON COLUMN master_accounts.personal_used IS '主账号个人配额已使用部分';

-- 5. 创建视图：用户配额汇总
-- 这个视图可以快速获取任何用户的总配额信息
CREATE OR REPLACE VIEW user_quota_summary AS
SELECT 
    u.id,
    u.username,
    u.account_type,
    u.balance_cpu_hours as personal_balance,
    CASE 
        WHEN u.account_type = 'master_account' THEN (
            SELECT COALESCE(SUM(allocated_quota - allocated_used), 0)
            FROM sub_accounts
            WHERE master_account_id = (
                SELECT id FROM master_accounts WHERE user_id = u.id
            )
        )
        ELSE 0
    END as allocated_balance,
    CASE 
        WHEN u.account_type = 'sub_account' THEN (
            SELECT COALESCE(allocated_quota - allocated_used, 0)
            FROM sub_accounts
            WHERE user_id = u.id
        )
        ELSE 0
    END as sub_allocated_balance,
    u.balance_cpu_hours + 
    CASE 
        WHEN u.account_type = 'master_account' THEN (
            SELECT COALESCE(SUM(allocated_quota - allocated_used), 0)
            FROM sub_accounts
            WHERE master_account_id = (
                SELECT id FROM master_accounts WHERE user_id = u.id
            )
        )
        WHEN u.account_type = 'sub_account' THEN (
            SELECT COALESCE(allocated_quota - allocated_used, 0)
            FROM sub_accounts
            WHERE user_id = u.id
        )
        ELSE 0
    END as total_available
FROM users u;

-- 6. 创建函数：获取用户的总可用配额
CREATE OR REPLACE FUNCTION get_user_total_quota(p_user_id INTEGER)
RETURNS FLOAT AS $$
DECLARE
    v_account_type VARCHAR;
    v_personal_balance FLOAT;
    v_allocated_balance FLOAT;
BEGIN
    SELECT account_type, balance_cpu_hours INTO v_account_type, v_personal_balance
    FROM users
    WHERE id = p_user_id;
    
    IF v_account_type = 'master_account' THEN
        -- 主账号：个人配额 + 分配给子账号的可用配额
        SELECT COALESCE(SUM(allocated_quota - allocated_used), 0)
        INTO v_allocated_balance
        FROM sub_accounts
        WHERE master_account_id = (
            SELECT id FROM master_accounts WHERE user_id = p_user_id
        );
        RETURN v_personal_balance + v_allocated_balance;
    ELSIF v_account_type = 'sub_account' THEN
        -- 子账号：个人配额 + 分配配额
        SELECT COALESCE(allocated_quota - allocated_used, 0)
        INTO v_allocated_balance
        FROM sub_accounts
        WHERE user_id = p_user_id;
        RETURN v_personal_balance + v_allocated_balance;
    ELSE
        -- 个人账户：只有个人配额
        RETURN v_personal_balance;
    END IF;
END;
$$ LANGUAGE plpgsql;

-- 7. 创建函数：为子账号分配配额
CREATE OR REPLACE FUNCTION allocate_quota_to_sub_account(
    p_sub_account_id INTEGER,
    p_allocated_quota FLOAT
)
RETURNS BOOLEAN AS $$
BEGIN
    UPDATE sub_accounts
    SET allocated_quota = p_allocated_quota,
        updated_at = CURRENT_TIMESTAMP
    WHERE id = p_sub_account_id;
    
    RETURN TRUE;
END;
$$ LANGUAGE plpgsql;

-- 8. 创建函数：消费子账号的配额
CREATE OR REPLACE FUNCTION consume_sub_account_quota(
    p_sub_account_id INTEGER,
    p_amount FLOAT
)
RETURNS BOOLEAN AS $$
DECLARE
    v_allocated_available FLOAT;
    v_personal_available FLOAT;
    v_total_needed FLOAT;
BEGIN
    SELECT allocated_quota - allocated_used, personal_quota - personal_used
    INTO v_allocated_available, v_personal_available
    FROM sub_accounts
    WHERE id = p_sub_account_id;
    
    v_total_needed := v_allocated_available + v_personal_available;
    
    IF v_total_needed < p_amount THEN
        RETURN FALSE;
    END IF;
    
    -- 先消费分配配额，再消费个人配额
    IF v_allocated_available >= p_amount THEN
        UPDATE sub_accounts
        SET allocated_used = allocated_used + p_amount,
            updated_at = CURRENT_TIMESTAMP
        WHERE id = p_sub_account_id;
    ELSE
        UPDATE sub_accounts
        SET allocated_used = allocated_quota,
            personal_used = personal_used + (p_amount - v_allocated_available),
            updated_at = CURRENT_TIMESTAMP
        WHERE id = p_sub_account_id;
    END IF;
    
    RETURN TRUE;
END;
$$ LANGUAGE plpgsql;

