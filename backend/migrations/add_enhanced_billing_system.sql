-- 增强版核时计费系统数据库迁移
-- 添加任务类型差异化定价和用户类型折扣支持

-- 1. 创建任务类型定价表
CREATE TABLE IF NOT EXISTS task_type_pricing (
    id SERIAL PRIMARY KEY,
    task_type VARCHAR(50) UNIQUE NOT NULL,
    price_per_hour FLOAT NOT NULL,
    is_active BOOLEAN DEFAULT TRUE NOT NULL,
    updated_by INTEGER REFERENCES users(id) ON DELETE SET NULL,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
);

COMMENT ON TABLE task_type_pricing IS '任务类型核时单价表';
COMMENT ON COLUMN task_type_pricing.task_type IS '任务类型: MD/QC/POSTPROCESS/REACTION_NETWORK/FORCEFIELD';
COMMENT ON COLUMN task_type_pricing.price_per_hour IS '核时单价(元/核时)';

-- 2. 创建用户类型折扣表
CREATE TABLE IF NOT EXISTS user_type_discounts (
    id SERIAL PRIMARY KEY,
    user_type VARCHAR(50) UNIQUE NOT NULL,
    discount_rate FLOAT NOT NULL,
    is_active BOOLEAN DEFAULT TRUE NOT NULL,
    updated_by INTEGER REFERENCES users(id) ON DELETE SET NULL,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
);

COMMENT ON TABLE user_type_discounts IS '用户类型折扣配置表';
COMMENT ON COLUMN user_type_discounts.user_type IS '用户类型: STUDENT/RESEARCHER/COMPANY/PREMIUM';
COMMENT ON COLUMN user_type_discounts.discount_rate IS '折扣率 (0.9表示9折, 1.0表示无折扣)';

-- 3. 为users表添加自定义折扣字段
ALTER TABLE users ADD COLUMN IF NOT EXISTS custom_discount_rate FLOAT;
COMMENT ON COLUMN users.custom_discount_rate IS '特定用户自定义折扣率(优先于用户类型折扣)';

-- 4. 插入默认任务类型单价
INSERT INTO task_type_pricing (task_type, price_per_hour, is_active) VALUES
('MD', 1.0, TRUE),
('QC', 2.0, TRUE),
('POSTPROCESS', 5.0, TRUE),
('REACTION_NETWORK', 10.0, TRUE),
('FORCEFIELD', 100.0, TRUE)
ON CONFLICT (task_type) DO NOTHING;

-- 5. 插入默认用户类型折扣
INSERT INTO user_type_discounts (user_type, discount_rate, is_active) VALUES
('STUDENT', 1.0, TRUE),
('RESEARCHER', 0.9, TRUE),
('COMPANY', 0.8, TRUE),
('PREMIUM', 0.7, TRUE)
ON CONFLICT (user_type) DO NOTHING;

-- 6. 创建索引
CREATE INDEX IF NOT EXISTS idx_task_type_pricing_task_type ON task_type_pricing(task_type);
CREATE INDEX IF NOT EXISTS idx_task_type_pricing_is_active ON task_type_pricing(is_active);
CREATE INDEX IF NOT EXISTS idx_user_type_discounts_user_type ON user_type_discounts(user_type);
CREATE INDEX IF NOT EXISTS idx_user_type_discounts_is_active ON user_type_discounts(is_active);

-- 7. 创建更新时间触发器
CREATE OR REPLACE FUNCTION update_updated_at_column()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = NOW();
    RETURN NEW;
END;
$$ language 'plpgsql';

CREATE TRIGGER update_task_type_pricing_updated_at BEFORE UPDATE ON task_type_pricing
    FOR EACH ROW EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_user_type_discounts_updated_at BEFORE UPDATE ON user_type_discounts
    FOR EACH ROW EXECUTE FUNCTION update_updated_at_column();
