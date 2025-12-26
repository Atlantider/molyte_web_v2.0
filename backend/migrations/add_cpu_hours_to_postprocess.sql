-- 为后处理任务表添加CPU核时字段
-- 这样可以追踪后处理任务的实际消耗

-- 1. 为 postprocess_jobs 表添加 CPU 核时字段
ALTER TABLE postprocess_jobs ADD COLUMN IF NOT EXISTS actual_cpu_hours FLOAT DEFAULT 0.0;
ALTER TABLE postprocess_jobs ADD COLUMN IF NOT EXISTS estimated_cpu_hours FLOAT DEFAULT 0.0;

-- 2. 创建索引以加快查询
CREATE INDEX IF NOT EXISTS idx_postprocess_jobs_cpu_hours ON postprocess_jobs(actual_cpu_hours);
CREATE INDEX IF NOT EXISTS idx_resp_jobs_cpu_hours ON resp_jobs(cpu_hours);

-- 3. 添加注释
COMMENT ON COLUMN postprocess_jobs.actual_cpu_hours IS '实际消耗的CPU核时（从Slurm获取）';
COMMENT ON COLUMN postprocess_jobs.estimated_cpu_hours IS '预估的CPU核时';

