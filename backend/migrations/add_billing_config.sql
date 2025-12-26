-- 计费配置表
-- 用于存储全局计费模式和配置

CREATE TABLE IF NOT EXISTS billing_config (
    id SERIAL PRIMARY KEY,
    pricing_mode VARCHAR(20) NOT NULL DEFAULT 'CORE_HOUR',  -- 'CORE_HOUR' 或 'TASK_TYPE'
    global_core_hour_price FLOAT DEFAULT 0.1,                -- 按核时计费时的全局单价
    updated_by INTEGER REFERENCES users(id),
    updated_at TIMESTAMP DEFAULT NOW(),
    created_at TIMESTAMP DEFAULT NOW(),
    CONSTRAINT check_pricing_mode CHECK (pricing_mode IN ('CORE_HOUR', 'TASK_TYPE'))
);

-- 插入默认配置
INSERT INTO billing_config (pricing_mode, global_core_hour_price)
VALUES ('CORE_HOUR', 0.1)
ON CONFLICT DO NOTHING;

-- 为防止重复计费,在md_jobs表添加唯一约束索引
-- 确保同一个任务只能有一条已结算的记录
CREATE UNIQUE INDEX IF NOT EXISTS idx_md_jobs_billed_unique 
ON md_jobs (id) 
WHERE billed = true;

-- 为qc_jobs添加同样的约束
CREATE UNIQUE INDEX IF NOT EXISTS idx_qc_jobs_billed_unique 
ON qc_jobs (id) 
WHERE billed = true;

-- 查看配置
SELECT * FROM billing_config;
