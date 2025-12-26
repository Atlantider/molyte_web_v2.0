-- 简化计费配置数据库重构
-- 移除折扣概念，直接存储价格

-- 1. 重命名user_type_discounts为user_type_prices
ALTER TABLE user_type_discounts RENAME TO user_type_prices;
ALTER TABLE user_type_prices RENAME COLUMN discount_rate TO core_hour_price;
COMMENT ON TABLE user_type_prices IS '用户类型核时单价配置';
COMMENT ON COLUMN user_type_prices.core_hour_price IS '核时单价(元/核时)';

-- 2. 更新task_type_pricing注释
COMMENT ON TABLE task_type_pricing IS '任务类型价格配置';
COMMENT ON COLUMN task_type_pricing.price_per_hour IS '任务价格(元/任务，在TASK_TYPE模式下不考虑核时)';

-- 3. 更新billing_config
ALTER TABLE billing_config DROP COLUMN IF EXISTS global_core_hour_price;
COMMENT ON COLUMN billing_config.pricing_mode IS '计费模式: CORE_HOUR=按核时计费(用户类型差异化), TASK_TYPE=按任务计费(固定价格)';

-- 4. 初始化用户类型价格（基于之前的折扣率转换为绝对价格）
-- 假设之前的全局单价是0.1元/核时
UPDATE user_type_prices SET core_hour_price = 0 WHERE user_type = 'ADMIN';
UPDATE user_type_prices SET core_hour_price = 1.0 WHERE user_type = 'GUEST';
UPDATE user_type_prices SET core_hour_price = 0.9 WHERE user_type = 'USER';
UPDATE user_type_prices SET core_hour_price = 0.7 WHERE user_type = 'PREMIUM';

-- 5. 确保所有用户类型都有价格记录
INSERT INTO user_type_prices (user_type, core_hour_price, is_active, updated_by)
VALUES 
    ('ADMIN', 0, TRUE, 1),
    ('GUEST', 1.0, TRUE, 1),
    ('USER', 0.9, TRUE, 1),
    ('PREMIUM', 0.7, TRUE, 1)
ON CONFLICT (user_type) DO UPDATE
SET is_active = TRUE;
