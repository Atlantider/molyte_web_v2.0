-- 添加定价追踪字段到 users 表
-- 用于记录用户自定义定价的更新时间和更新者

ALTER TABLE users ADD COLUMN IF NOT EXISTS price_updated_at TIMESTAMP WITH TIME ZONE;
ALTER TABLE users ADD COLUMN IF NOT EXISTS price_updated_by INTEGER;

-- 完成
SELECT 'Price tracking fields added successfully!' as message;

