-- 添加用户级别定价字段
-- 允许管理员为不同用户设置不同的核时单价

-- 1. 添加用户级别定价字段到 users 表
ALTER TABLE users ADD COLUMN IF NOT EXISTS custom_cpu_hour_price FLOAT;
ALTER TABLE users ADD COLUMN IF NOT EXISTS price_updated_at TIMESTAMP WITH TIME ZONE;
ALTER TABLE users ADD COLUMN IF NOT EXISTS price_updated_by INTEGER;

-- 2. 添加外键约束（price_updated_by 指向管理员用户）
ALTER TABLE users ADD CONSTRAINT fk_price_updated_by 
  FOREIGN KEY (price_updated_by) REFERENCES users(id) ON DELETE SET NULL;

-- 3. 创建索引以加快查询
CREATE INDEX IF NOT EXISTS idx_users_custom_price ON users(custom_cpu_hour_price) 
  WHERE custom_cpu_hour_price IS NOT NULL;

-- 4. 添加系统配置说明（如果不存在）
INSERT INTO system_configs (key, value, description) 
VALUES ('user_pricing_enabled', 'true', '是否启用用户级别定价功能')
ON CONFLICT (key) DO NOTHING;

-- 5. 创建审计日志表（可选，用于记录定价变更）
CREATE TABLE IF NOT EXISTS pricing_audit_log (
    id SERIAL PRIMARY KEY,
    admin_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    old_price FLOAT,
    new_price FLOAT,
    reason VARCHAR(255),
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    
    CONSTRAINT check_prices CHECK (old_price IS NULL OR old_price >= 0),
    CONSTRAINT check_new_price CHECK (new_price IS NULL OR new_price >= 0)
);

-- 6. 创建索引
CREATE INDEX IF NOT EXISTS idx_pricing_audit_user ON pricing_audit_log(user_id);
CREATE INDEX IF NOT EXISTS idx_pricing_audit_admin ON pricing_audit_log(admin_id);
CREATE INDEX IF NOT EXISTS idx_pricing_audit_created ON pricing_audit_log(created_at);

