-- QC 自动修复规则表迁移脚本
-- 用于存储 QC 错误自动修复规则

-- 创建错误规则表
CREATE TABLE IF NOT EXISTS qc_error_rules (
    id SERIAL PRIMARY KEY,
    priority INT NOT NULL,
    error_category VARCHAR(50) NOT NULL,
    error_patterns TEXT NOT NULL,  -- JSON 数组，例如 ["pattern1", "pattern2"]
    description TEXT NOT NULL,
    max_retries INT NOT NULL DEFAULT 2,
    requires_manual_check BOOLEAN NOT NULL DEFAULT FALSE,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE(priority, error_category)
);

-- 创建修复策略表
CREATE TABLE IF NOT EXISTS qc_repair_strategies (
    id SERIAL PRIMARY KEY,
    rule_id INT NOT NULL REFERENCES qc_error_rules(id) ON DELETE CASCADE,
    strategy_index INT NOT NULL,
    name VARCHAR(100) NOT NULL,
    keywords_to_add TEXT NOT NULL,
    keywords_to_remove TEXT,  -- JSON 数组
    restart_from_chk BOOLEAN NOT NULL DEFAULT FALSE,
    use_cartesian BOOLEAN NOT NULL DEFAULT FALSE,
    description TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE(rule_id, strategy_index)
);

-- 创建任务重试历史表
CREATE TABLE IF NOT EXISTS qc_job_retry_history (
    id SERIAL PRIMARY KEY,
    job_id INT NOT NULL,
    error_type VARCHAR(50) NOT NULL,
    error_count INT NOT NULL DEFAULT 0,
    total_retries INT NOT NULL DEFAULT 0,
    last_strategy_index INT,
    last_error_message TEXT,
    last_repair_action VARCHAR(50),  -- RETRY, MANUAL_CHECK, etc.
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (job_id) REFERENCES qc_jobs(id) ON DELETE CASCADE
);

-- 创建索引以加快查询
CREATE INDEX IF NOT EXISTS idx_qc_error_rules_priority ON qc_error_rules(priority);
CREATE INDEX IF NOT EXISTS idx_qc_error_rules_category ON qc_error_rules(error_category);
CREATE INDEX IF NOT EXISTS idx_qc_repair_strategies_rule_id ON qc_repair_strategies(rule_id);
CREATE INDEX IF NOT EXISTS idx_qc_job_retry_history_job_id ON qc_job_retry_history(job_id);
CREATE INDEX IF NOT EXISTS idx_qc_job_retry_history_error_type ON qc_job_retry_history(error_type);

-- 插入默认规则数据
-- 优先级 1: 内部坐标系统崩溃（Tors failed + FormBX + 段错误）
INSERT INTO qc_error_rules (priority, error_category, error_patterns, description, max_retries, requires_manual_check)
VALUES (
    1,
    'geometry_optimization',
    '["Tors failed|Bend failed", "FormBX had a problem", "segmentation violation|segmentation fault"]',
    '内部坐标系统崩溃（Tors failed + FormBX + 段错误）',
    2,
    TRUE
) ON CONFLICT (priority, error_category) DO NOTHING;

-- 插入该规则的修复策略
INSERT INTO qc_repair_strategies (rule_id, strategy_index, name, keywords_to_add, restart_from_chk, use_cartesian, description)
SELECT id, 0, '笛卡尔坐标重启（第1级）',
    'Opt=(Cartesian,MaxStep=5) SCF=XQC Integral=UltraFine NoSymm Geom=AllCheck Guess=Read',
    TRUE, TRUE, '从 chk 文件重启，强制使用笛卡尔坐标，减小步长'
FROM qc_error_rules WHERE priority = 1 AND error_category = 'geometry_optimization'
ON CONFLICT (rule_id, strategy_index) DO NOTHING;

INSERT INTO qc_repair_strategies (rule_id, strategy_index, name, keywords_to_add, restart_from_chk, use_cartesian, description)
SELECT id, 1, '笛卡尔坐标重启（第2级）',
    'Opt=(Cartesian,CalcFC,MaxStep=5) SCF=XQC Integral=UltraFine NoSymm Geom=AllCheck Guess=Read',
    TRUE, TRUE, '计算初始 Hessian，进一步稳定优化'
FROM qc_error_rules WHERE priority = 1 AND error_category = 'geometry_optimization'
ON CONFLICT (rule_id, strategy_index) DO NOTHING;

-- 优先级 4: SCF 不收敛
INSERT INTO qc_error_rules (priority, error_category, error_patterns, description, max_retries, requires_manual_check)
VALUES (
    4,
    'scf_convergence',
    '["No convergence in SCF|SCF has not converged|Convergence failure in SCF"]',
    'SCF 迭代未收敛',
    2,
    TRUE
) ON CONFLICT (priority, error_category) DO NOTHING;

INSERT INTO qc_repair_strategies (rule_id, strategy_index, name, keywords_to_add, restart_from_chk, use_cartesian, description)
SELECT id, 0, 'SCF 收敛策略（第1级）',
    'SCF=(XQC,MaxCycle=512)',
    FALSE, FALSE, '使用二次收敛算法，增加迭代次数'
FROM qc_error_rules WHERE priority = 4 AND error_category = 'scf_convergence'
ON CONFLICT (rule_id, strategy_index) DO NOTHING;

INSERT INTO qc_repair_strategies (rule_id, strategy_index, name, keywords_to_add, restart_from_chk, use_cartesian, description)
SELECT id, 1, 'SCF 收敛策略（第2级）',
    'SCF=(XQC,MaxCycle=512) Guess=Mix Integral=UltraFine',
    FALSE, FALSE, '添加混合初始猜测，超精细积分'
FROM qc_error_rules WHERE priority = 4 AND error_category = 'scf_convergence'
ON CONFLICT (rule_id, strategy_index) DO NOTHING;

-- 优先级 8: 积分精度不足
INSERT INTO qc_error_rules (priority, error_category, error_patterns, description, max_retries, requires_manual_check)
VALUES (
    8,
    'integral_precision',
    '["Integral accuracy insufficient|Tried to raise integral accuracy"]',
    '积分精度不足',
    2,
    FALSE
) ON CONFLICT (priority, error_category) DO NOTHING;

INSERT INTO qc_repair_strategies (rule_id, strategy_index, name, keywords_to_add, restart_from_chk, use_cartesian, description)
SELECT id, 0, '提高积分精度（第1级）',
    'Integral=UltraFine',
    FALSE, FALSE, '使用超精细积分网格'
FROM qc_error_rules WHERE priority = 8 AND error_category = 'integral_precision'
ON CONFLICT (rule_id, strategy_index) DO NOTHING;

INSERT INTO qc_repair_strategies (rule_id, strategy_index, name, keywords_to_add, restart_from_chk, use_cartesian, description)
SELECT id, 1, '提高积分精度（第2级）',
    'Integral=SuperFineGrid',
    FALSE, FALSE, '使用超级精细积分网格'
FROM qc_error_rules WHERE priority = 8 AND error_category = 'integral_precision'
ON CONFLICT (rule_id, strategy_index) DO NOTHING;

-- 创建视图：获取规则及其策略
CREATE OR REPLACE VIEW v_qc_error_rules_with_strategies AS
SELECT
    r.id as rule_id,
    r.priority,
    r.error_category,
    r.error_patterns,
    r.description as rule_description,
    r.max_retries,
    r.requires_manual_check,
    s.id as strategy_id,
    s.strategy_index,
    s.name as strategy_name,
    s.keywords_to_add,
    s.keywords_to_remove,
    s.restart_from_chk,
    s.use_cartesian,
    s.description as strategy_description
FROM qc_error_rules r
LEFT JOIN qc_repair_strategies s ON r.id = s.rule_id
ORDER BY r.priority, s.strategy_index;

-- 创建视图：获取任务的重试统计
CREATE OR REPLACE VIEW v_qc_job_retry_stats AS
SELECT
    job_id,
    COUNT(*) as total_retry_attempts,
    MAX(total_retries) as max_total_retries,
    COUNT(DISTINCT error_type) as unique_error_types,
    MAX(updated_at) as last_retry_time
FROM qc_job_retry_history
GROUP BY job_id;

