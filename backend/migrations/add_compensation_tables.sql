-- ============================================================================
-- 创建自动化补偿规则表
-- Migration: Add compensation rules and records tables
-- Date: 2025-12-17
-- ============================================================================

-- 1. 创建补偿规则类型枚举
DO $$ BEGIN
    CREATE TYPE compensation_rule_type AS ENUM (
        'JOB_FAILURE',
        'EXPIRATION',
        'MANUAL'
    );
EXCEPTION
    WHEN duplicate_object THEN null;
END $$;

-- 2. 创建补偿状态枚举
DO $$ BEGIN
    CREATE TYPE compensation_status AS ENUM (
        'PENDING',
        'APPROVED',
        'REJECTED',
        'COMPLETED'
    );
EXCEPTION
    WHEN duplicate_object THEN null;
END $$;

-- 3. 创建补偿规则表
CREATE TABLE IF NOT EXISTS compensation_rules (
    id SERIAL PRIMARY KEY,
    name VARCHAR(200) NOT NULL,
    description TEXT,
    
    -- 规则类型和配置
    rule_type compensation_rule_type NOT NULL,
    is_active BOOLEAN DEFAULT TRUE NOT NULL,
    
    -- 规则配置（JSON格式）
    config JSONB DEFAULT '{}' NOT NULL,
    
    -- 时间戳
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP NOT NULL,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_comp_rules_type ON compensation_rules(rule_type);
CREATE INDEX IF NOT EXISTS idx_comp_rules_active ON compensation_rules(is_active);

-- 4. 创建补偿记录表
CREATE TABLE IF NOT EXISTS compensation_records (
    id SERIAL PRIMARY KEY,
    rule_id INTEGER REFERENCES compensation_rules(id) ON DELETE SET NULL,
    user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    
    -- 补偿信息
    status compensation_status DEFAULT 'PENDING' NOT NULL,
    reason VARCHAR(500) NOT NULL,
    amount FLOAT NOT NULL,
    
    -- 关联信息
    reference_id INTEGER,
    reference_type VARCHAR(50),
    
    -- 审批信息
    approved_by INTEGER REFERENCES users(id) ON DELETE SET NULL,
    approval_reason VARCHAR(500),
    
    -- 时间戳
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP NOT NULL,
    approved_at TIMESTAMP WITH TIME ZONE,
    completed_at TIMESTAMP WITH TIME ZONE
);

CREATE INDEX IF NOT EXISTS idx_comp_records_user_id ON compensation_records(user_id);
CREATE INDEX IF NOT EXISTS idx_comp_records_status ON compensation_records(status);
CREATE INDEX IF NOT EXISTS idx_comp_records_created_at ON compensation_records(created_at);

-- 5. 创建核时过期记录表
CREATE TABLE IF NOT EXISTS cpu_hours_expirations (
    id SERIAL PRIMARY KEY,
    user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    
    -- 过期信息
    amount FLOAT NOT NULL,
    reason VARCHAR(500) NOT NULL,
    
    -- 时间戳
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP NOT NULL,
    expired_at TIMESTAMP WITH TIME ZONE NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_cpu_exp_user_id ON cpu_hours_expirations(user_id);
CREATE INDEX IF NOT EXISTS idx_cpu_exp_expired_at ON cpu_hours_expirations(expired_at);

-- 6. 添加注释
COMMENT ON TABLE compensation_rules IS '补偿规则配置表 - 定义自动补偿的规则';
COMMENT ON TABLE compensation_records IS '补偿记录表 - 记录所有补偿操作';
COMMENT ON TABLE cpu_hours_expirations IS '核时过期记录表 - 记录过期的核时';

COMMENT ON COLUMN compensation_rules.config IS '规则配置，JSON格式，根据rule_type不同而不同';
COMMENT ON COLUMN compensation_records.status IS '补偿状态：待处理/已批准/已拒绝/已完成';
COMMENT ON COLUMN compensation_records.reference_type IS '关联类型：job/order/manual等';

-- 7. 插入默认补偿规则
INSERT INTO compensation_rules (name, description, rule_type, is_active, config)
VALUES (
    '任务失败自动退款',
    '当任务失败时，自动退款给用户',
    'JOB_FAILURE',
    TRUE,
    '{"failure_types": ["FAILED", "TIMEOUT"], "refund_percentage": 100, "min_cpu_hours": 0.1}'::jsonb
)
ON CONFLICT DO NOTHING;

INSERT INTO compensation_rules (name, description, rule_type, is_active, config)
VALUES (
    '核时过期规则',
    '365天未使用的核时自动过期',
    'EXPIRATION',
    TRUE,
    '{"expiration_days": 365, "warning_days": 30}'::jsonb
)
ON CONFLICT DO NOTHING;

