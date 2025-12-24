-- 清理QC任务表中的重复核时字段
-- 删除未使用的cpu_hours和resp_cpu_hours字段
-- 保留actual_cpu_hours作为唯一的核时字段

-- 删除未使用的cpu_hours字段
ALTER TABLE qc_jobs DROP COLUMN IF EXISTS cpu_hours;

-- 删除未使用的resp_cpu_hours字段（QC任务不需要RESP核时）
ALTER TABLE qc_jobs DROP COLUMN IF EXISTS resp_cpu_hours;

-- 确保actual_cpu_hours字段存在且有正确的定义
ALTER TABLE qc_jobs ADD COLUMN IF NOT EXISTS actual_cpu_hours FLOAT DEFAULT 0.0 NOT NULL;

-- 添加注释说明
COMMENT ON COLUMN qc_jobs.actual_cpu_hours IS 'QC任务在Slurm上实际运行的CPU核时总和（从sacct获取的CPUTimeRAW，单位：小时）';

