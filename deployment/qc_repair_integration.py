"""
QC 自动修复规则引擎集成模块

这个模块展示如何在 Polling Worker 中集成自动修复规则引擎。
"""

import re
import logging
from typing import Dict, Optional, Tuple
from pathlib import Path
from qc_auto_repair_engine import QCAutoRepairEngine, ErrorRule, RepairStrategy

logger = logging.getLogger(__name__)


class QCRepairIntegration:
    """QC 修复规则集成"""

    def __init__(self):
        """初始化"""
        self.engine = QCAutoRepairEngine()
        self.retry_history: Dict[int, Dict] = {}  # job_id -> {error_type, count, total_retries, last_strategy}

    def analyze_and_repair(self, job_id: int, log_file_path: str, original_route: str) -> Optional[Dict]:
        """
        分析错误并生成修复方案

        Args:
            job_id: 任务 ID
            log_file_path: Gaussian 日志文件路径
            original_route: 原始 route 行

        Returns:
            修复方案字典或 None
        """
        try:
            # 1. 读取日志文件
            with open(log_file_path, 'r', encoding='utf-8', errors='ignore') as f:
                log_content = f.read()

            # 2. 分析错误
            result = self.engine.analyze_error(log_content)
            if not result:
                logger.warning(f"Job {job_id}: 无法识别错误类型")
                return None

            rule, _ = result
            logger.info(f"Job {job_id}: 匹配规则 - {rule.description}")

            # 3. 获取重试历史
            history = self.retry_history.get(job_id, {
                'error_type': rule.error_category.value,
                'error_count': 0,
                'total_retries': 0,
                'last_strategy_index': -1
            })

            # 4. 判断是否应该继续重试
            if not self.engine.should_continue_retry(
                history['error_count'],
                history['total_retries']
            ):
                logger.warning(f"Job {job_id}: 已达到最大重试次数，需要人工检查")
                return {
                    'action': 'MANUAL_CHECK',
                    'reason': f"已达到最大重试次数（{history['total_retries']} 次）",
                    'error_type': rule.error_category.value,
                    'requires_manual_check': rule.requires_manual_check
                }

            # 5. 获取修复策略
            strategy_index = history['error_count']
            strategy = self.engine.get_repair_strategy(rule, strategy_index)

            if not strategy:
                logger.warning(f"Job {job_id}: 没有更多的修复策略")
                return {
                    'action': 'MANUAL_CHECK',
                    'reason': '没有更多的修复策略',
                    'error_type': rule.error_category.value,
                    'requires_manual_check': True
                }

            logger.info(f"Job {job_id}: 应用修复策略 - {strategy.name}")

            # 6. 生成新的 route 行
            new_route = self._apply_strategy(original_route, strategy)

            # 7. 更新重试历史
            history['error_count'] += 1
            history['total_retries'] += 1
            history['last_strategy_index'] = strategy_index
            self.retry_history[job_id] = history

            # 8. 返回修复方案
            return {
                'action': 'RETRY',
                'new_route': new_route,
                'strategy_name': strategy.name,
                'strategy_description': strategy.description,
                'restart_from_chk': strategy.restart_from_chk,
                'use_cartesian': strategy.use_cartesian,
                'error_type': rule.error_category.value,
                'error_count': history['error_count'],
                'total_retries': history['total_retries'],
                'requires_manual_check': rule.requires_manual_check
            }

        except Exception as e:
            logger.error(f"Job {job_id}: 分析和修复失败 - {e}", exc_info=True)
            return None

    def _apply_strategy(self, original_route: str, strategy: RepairStrategy) -> str:
        """
        应用修复策略到 route 行

        Args:
            original_route: 原始 route 行
            strategy: 修复策略

        Returns:
            修改后的 route 行
        """
        # 移除原有的关键词
        new_route = original_route
        if strategy.keywords_to_remove:
            for keyword in strategy.keywords_to_remove:
                # 使用正则表达式移除关键词（支持带参数的关键词）
                pattern = rf'\b{re.escape(keyword)}(?:=[^\s/]*)?\b'
                new_route = re.sub(pattern, '', new_route, flags=re.IGNORECASE)

        # 添加新的关键词
        if strategy.keywords_to_add:
            # 在 route 行末尾添加新关键词
            new_route = new_route.rstrip() + ' ' + strategy.keywords_to_add

        # 清理多余的空格
        new_route = ' '.join(new_route.split())

        return new_route

    def get_retry_history(self, job_id: int) -> Dict:
        """获取任务的重试历史"""
        return self.retry_history.get(job_id, {})

    def clear_retry_history(self, job_id: int):
        """清除任务的重试历史"""
        if job_id in self.retry_history:
            del self.retry_history[job_id]


# ==================== 集成到 Polling Worker 的示例 ====================

def integrate_with_polling_worker():
    """
    这是一个示例，展示如何在 Polling Worker 中使用修复规则引擎

    在 polling_worker.py 的 _handle_job_completion 方法中，
    当检测到任务失败时，调用这个函数：

    ```python
    # 在 polling_worker.py 中
    from qc_repair_integration import QCRepairIntegration

    class PollingWorker:
        def __init__(self, ...):
            ...
            self.repair_integration = QCRepairIntegration()

        def _handle_job_completion(self, job_id, job_info):
            ...
            # 检查任务是否失败
            if job_status == 'FAILED':
                # 分析错误并生成修复方案
                repair_plan = self.repair_integration.analyze_and_repair(
                    job_id,
                    log_file_path,
                    original_route
                )

                if repair_plan:
                    if repair_plan['action'] == 'RETRY':
                        # 应用修复策略，提交新任务
                        self._submit_repaired_job(
                            job_id,
                            repair_plan['new_route'],
                            repair_plan['restart_from_chk'],
                            repair_plan['use_cartesian']
                        )
                    elif repair_plan['action'] == 'MANUAL_CHECK':
                        # 标记为需要人工检查
                        self._update_job_status(
                            job_id,
                            'FAILED',
                            error_message=repair_plan['reason']
                        )
    ```
    """
    pass


# ==================== 测试函数 ====================

def test_repair_integration():
    """测试修复规则集成"""
    integration = QCRepairIntegration()

    # 测试 1: 内部坐标崩溃
    print("=" * 60)
    print("测试 1: 内部坐标崩溃")
    print("=" * 60)

    log_content_1 = """
    Tors failed for angle 1-2-3-4
    FormBX had a problem
    segmentation violation
    Error termination via Lnk1e
    """

    original_route_1 = "#p B3LYP/6-31G(d) Opt Freq"

    result_1 = integration.analyze_and_repair(
        job_id=1001,
        log_file_path="/tmp/test_1.log",
        original_route=original_route_1
    )

    # 创建临时日志文件用于测试
    Path("/tmp/test_1.log").write_text(log_content_1)

    result_1 = integration.analyze_and_repair(
        job_id=1001,
        log_file_path="/tmp/test_1.log",
        original_route=original_route_1
    )

    if result_1:
        print(f"Action: {result_1['action']}")
        print(f"Strategy: {result_1.get('strategy_name', 'N/A')}")
        print(f"New Route: {result_1.get('new_route', 'N/A')}")
        print(f"Restart from CHK: {result_1.get('restart_from_chk', False)}")
    else:
        print("Failed to analyze error")

    # 测试 2: SCF 不收敛
    print("\n" + "=" * 60)
    print("测试 2: SCF 不收敛")
    print("=" * 60)

    log_content_2 = """
    No convergence in SCF after 128 cycles
    Error termination via Lnk1e
    """

    original_route_2 = "#p B3LYP/6-31G(d) Opt"

    Path("/tmp/test_2.log").write_text(log_content_2)

    result_2 = integration.analyze_and_repair(
        job_id=1002,
        log_file_path="/tmp/test_2.log",
        original_route=original_route_2
    )

    if result_2:
        print(f"Action: {result_2['action']}")
        print(f"Strategy: {result_2.get('strategy_name', 'N/A')}")
        print(f"New Route: {result_2.get('new_route', 'N/A')}")
    else:
        print("Failed to analyze error")


if __name__ == '__main__':
    # 配置日志
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s | %(levelname)s | %(name)s | %(message)s'
    )

    # 运行测试
    test_repair_integration()

