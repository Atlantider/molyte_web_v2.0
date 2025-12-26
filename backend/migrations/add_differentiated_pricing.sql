-- 差异化定价系统数据库迁移
-- 日期: 2025-12-26
-- 说明: 添加任务类型定价和用户折扣功能

-- 1. 创建任务类型定价表
CREATE TABLE IF NOT EXISTS task_type_pricing (
    id SERIAL PRIMARY KEY,
    task_type VARCHAR(50) UNIQUE NOT NULL,  -- MD/QC/POSTPROCESS/REACTION_NETWORK/FORCEFIELD
    price_per_hour FLOAT NOT NULL,          -- 核时单价(元/核时)
    description TEXT,                        -- 描述
    is_active BOOLEAN DEFAULT TRUE,
    updated_by INTEGER REFERENCES users(id),
    updated_at TIMESTAMP DEFAULT NOW(),
    created_at TIMESTAMP DEFAULT NOW()
);

-- 添加索引
CREATE INDEX IF NOT EXISTS idx_task_type_pricing_type ON task_type_pricing(task_type);
CREATE INDEX IF NOT EXISTS idx_task_type_pricing_active ON task_type_pricing(is_active);

-- 2. 创建用户类型折扣表
CREATE TABLE IF NOT EXISTS user_type_discounts (
    id SERIAL PRIMARY KEY,
    user_type VARCHAR(50) UNIQUE NOT NULL,  -- STUDENT/RESEARCHER/COMPANY/PREMIUM
    discount_rate FLOAT NOT NULL,            -- 折扣率 (0.9表示9折, 1.0表示无折扣)
    description TEXT,                        -- 描述
    is_active BOOLEAN DEFAULT TRUE,
    updated_by INTEGER REFERENCES users(id),
    updated_at TIMESTAMP DEFAULT NOW(),
    created_at TIMESTAMP DEFAULT NOW(),
    CONSTRAINT check_discount_rate CHECK (discount_rate > 0 AND discount_rate <= 1.0)
);

-- 添加索引
CREATE INDEX IF NOT EXISTS idx_user_type_discounts_type ON user_type_discounts(user_type);
CREATE INDEX IF NOT EXISTS idx_user_type_discounts_active ON user_type_discounts(is_active);

-- 3. 为users表添加自定义折扣率字段
DO $$
BEGIN
    IF NOT EXISTS (
        SELECT 1 FROM information_schema.columns 
        WHERE table_name = 'users' AND column_name = 'custom_discount_rate'
    ) THEN
        ALTER TABLE users ADD COLUMN custom_discount_rate FLOAT;
        ALTER TABLE users ADD CONSTRAINT check_custom_discount_rate 
            CHECK (custom_discount_rate IS NULL OR (custom_discount_rate > 0 AND custom_discount_rate <= 1.0));
    END IF;
END $$;

-- 4. 插入默认任务类型单价
INSERT INTO task_type_pricing (task_type, price_per_hour, description) VALUES
('MD', 1.0, 'MD分子动力学计算'),
('QC', 2.0, 'QC量子化学计算'),
('POSTPROCESS', 5.0, '后处理分析(RDF, MSD等)'),
('REACTION_NETWORK', 10.0, '反应网络生成'),
('FORCEFIELD', 100.0, '力场参数生成')
ON CONFLICT (task_type) DO UPDATE SET
    price_per_hour = EXCLUDED.price_per_hour,
    description = EXCLUDED.description,
    updated_at = NOW();

-- 5. 插入默认用户类型折扣
INSERT INTO user_type_discounts (user_type, discount_rate, description) VALUES
('STUDENT', 1.0, '学生用户 - 标准价格'),
('RESEARCHER', 0.9, '研究人员 - 9折优惠'),
('COMPANY', 0.8, '企业用户 - 8折优惠'),
('PREMIUM', 0.7, '高级用户 - 7折优惠')
ON CONFLICT (user_type) DO UPDATE SET
    discount_rate = EXCLUDED.discount_rate,
    description = EXCLUDED.description,
    updated_at = NOW();

-- 6. 验证数据
SELECT '任务类型定价:' as info;
SELECT task_type, price_per_hour, description FROM task_type_pricing ORDER BY price_per_hour;

SELECT '用户类型折扣:' as info;
SELECT user_type, discount_rate, description FROM user_type_discounts ORDER BY discount_rate DESC;
