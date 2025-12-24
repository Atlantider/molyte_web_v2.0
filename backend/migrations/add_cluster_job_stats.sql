-- Add cpu_hours_used and task_count columns to advanced_cluster_jobs table
-- This migration adds fields to track CPU hours and task count for billing and statistics

ALTER TABLE advanced_cluster_jobs
ADD COLUMN IF NOT EXISTS cpu_hours_used FLOAT DEFAULT 0.0 NOT NULL,
ADD COLUMN IF NOT EXISTS task_count INTEGER DEFAULT 0 NOT NULL;

-- Create index for task_count for efficient querying
CREATE INDEX IF NOT EXISTS idx_advanced_cluster_task_count ON advanced_cluster_jobs(task_count);

-- Create index for cpu_hours_used for efficient querying
CREATE INDEX IF NOT EXISTS idx_advanced_cluster_cpu_hours ON advanced_cluster_jobs(cpu_hours_used);

-- Create index for user_id and created_at for daily statistics
CREATE INDEX IF NOT EXISTS idx_advanced_cluster_user_created ON advanced_cluster_jobs(user_id, created_at);

