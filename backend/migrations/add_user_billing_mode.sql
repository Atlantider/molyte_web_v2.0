-- Migration: Add billing_mode and custom_task_prices to users table
-- Date: 2025-12-26
-- Description: Add per-user billing mode configuration (CORE_HOUR or TASK_TYPE)
--              and custom task prices (JSON) for TASK_TYPE mode

-- Add billing_mode column (default to CORE_HOUR for existing users)
ALTER TABLE users 
ADD COLUMN IF NOT EXISTS billing_mode VARCHAR(20) NOT NULL DEFAULT 'CORE_HOUR';

-- Add custom_task_prices column (JSON for storing task-specific custom prices)
ALTER TABLE users 
ADD COLUMN IF NOT EXISTS custom_task_prices JSONB;

-- Add comment for documentation
COMMENT ON COLUMN users.billing_mode IS '计费模式: CORE_HOUR=按核时, TASK_TYPE=按任务类型';
COMMENT ON COLUMN users.custom_task_prices IS '自定义任务价格: {MD: 10.0, QC: 50.0, ...}'；
