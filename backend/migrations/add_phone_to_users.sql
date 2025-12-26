-- 添加手机号相关字段到 users 表
-- 2025-01-18: 支持手机号注册和验证

-- 添加手机号字段
ALTER TABLE users ADD COLUMN IF NOT EXISTS phone VARCHAR(20) UNIQUE;

-- 添加手机号验证状态字段
ALTER TABLE users ADD COLUMN IF NOT EXISTS phone_verified BOOLEAN DEFAULT FALSE NOT NULL;

-- 添加真实姓名字段（如果不存在）
ALTER TABLE users ADD COLUMN IF NOT EXISTS real_name VARCHAR(50);

-- 创建手机号索引
CREATE INDEX IF NOT EXISTS idx_users_phone ON users(phone);

-- 为已有用户设置默认值
UPDATE users SET phone_verified = FALSE WHERE phone_verified IS NULL;

COMMENT ON COLUMN users.phone IS '用户手机号';
COMMENT ON COLUMN users.phone_verified IS '手机号是否已验证';
COMMENT ON COLUMN users.real_name IS '用户真实姓名';

