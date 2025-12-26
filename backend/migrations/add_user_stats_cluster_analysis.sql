-- Add cluster analysis task count and CPU hours tracking to user_usage_stats table
-- This migration adds fields to track cluster analysis job statistics

ALTER TABLE user_usage_stats
ADD COLUMN IF NOT EXISTS cluster_analysis_task_count INTEGER DEFAULT 0 NOT NULL,
ADD COLUMN IF NOT EXISTS cluster_analysis_cpu_hours FLOAT DEFAULT 0.0 NOT NULL,
ADD COLUMN IF NOT EXISTS max_concurrent_jobs INTEGER DEFAULT 0 NOT NULL;

-- Create indexes for efficient querying
CREATE INDEX IF NOT EXISTS idx_user_stats_cluster_task_count ON user_usage_stats(cluster_analysis_task_count);
CREATE INDEX IF NOT EXISTS idx_user_stats_cluster_cpu_hours ON user_usage_stats(cluster_analysis_cpu_hours);
CREATE INDEX IF NOT EXISTS idx_user_stats_max_concurrent ON user_usage_stats(max_concurrent_jobs);

