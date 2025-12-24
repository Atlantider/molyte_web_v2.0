-- Add CPU hours tracking to QC jobs
-- This allows accurate tracking of actual Slurm CPU hours used by each QC task

-- 添加 actual_cpu_hours 字段（真实的 Slurm CPUTimeRAW）
ALTER TABLE qc_jobs
ADD COLUMN IF NOT EXISTS actual_cpu_hours FLOAT DEFAULT 0.0 NOT NULL;

-- 添加 resp_cpu_hours 字段（RESP 电荷计算消耗的 CPU 核时数）
ALTER TABLE qc_jobs ADD COLUMN IF NOT EXISTS resp_cpu_hours FLOAT DEFAULT 0.0;

-- Add index for efficient querying
CREATE INDEX IF NOT EXISTS idx_qc_jobs_cpu_hours ON qc_jobs(actual_cpu_hours);
CREATE INDEX IF NOT EXISTS idx_qc_jobs_cluster_analysis ON qc_jobs(cluster_analysis_job_id, actual_cpu_hours);

