-- 修复用户余额计算
-- 问题：balance_cpu_hours应该 = free_cpu_hours_granted + total_recharged - total_consumed
-- 但由于迁移脚本的问题，某些用户的balance_cpu_hours被设置为了错误的值

-- 对于管理员用户，由于不扣费，balance_cpu_hours应该 = free_cpu_hours_granted + total_recharged

-- 创建临时表来计算正确的余额
CREATE TEMP TABLE correct_balances AS
SELECT 
    u.id,
    u.balance_cpu_hours as old_balance,
    u.free_cpu_hours_granted +
    COALESCE((
        SELECT SUM(amount) FROM quota_transactions 
        WHERE user_id = u.id 
        AND type IN ('recharge', 'points_exchange', 'admin_adjust', 'bonus')
    ), 0.0) -
    COALESCE((
        SELECT ABS(SUM(amount)) FROM quota_transactions 
        WHERE user_id = u.id 
        AND type = 'consume'
    ), 0.0) as calculated_balance
FROM users u;

-- 显示需要修复的用户
SELECT id, old_balance, calculated_balance 
FROM correct_balances 
WHERE old_balance != calculated_balance
ORDER BY id;

-- 修复余额（可选，取消注释以执行）
-- UPDATE users u
-- SET balance_cpu_hours = cb.calculated_balance
-- FROM correct_balances cb
-- WHERE u.id = cb.id AND u.balance_cpu_hours != cb.calculated_balance;

-- 验证修复
-- SELECT id, balance_cpu_hours, free_cpu_hours_granted FROM users WHERE id IN (SELECT id FROM correct_balances WHERE old_balance != calculated_balance);

