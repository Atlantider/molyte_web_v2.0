#!/usr/bin/env python3
"""
混合云轮询 Worker (v2.0 模块化重构版)

功能：
1. 定期轮询腾讯云 API，获取待处理任务
2. 下载任务输入数据
3. 生成 LAMMPS/Gaussian 输入文件
4. 提交到 Slurm 集群
5. 监控任务状态
6. 上传结果到 COS（支持部分结果）

重构特性：
- 模块化设计，代码清晰
- 支持部分结果上传
- 云端重试协调
- 可测试的组件
"""
import sys
import logging
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent / "backend"))
sys.path.insert(0, str(Path(__file__).parent))

# 导入模块
from worker.core import WorkerConfig, APIClient, TaskScheduler
from worker.handlers import (
    MDHandler,
    QCHandler,
    PostprocessHandler,
    BindingHandler,
    RedoxHandler,
    ClusterHandler,
    AnionHandler,
)
from worker.coordinators import JobTracker, RetryCoordinator
from worker.uploaders import COSClient

logger = logging.getLogger(__name__)


class PollingWorkerV2:
    """轮询 Worker 主类 (v2.0 模块化版本)"""
    
    def __init__(self, config_path: str = "deployment/polling_worker_config.yaml"):
        """
        初始化 Worker
        
        Args:
            config_path: 配置文件路径
        """
        # 1. 加载配置
        self.config = WorkerConfig(config_path)
        
        # 2. 初始化 API 客户端
        self.client = APIClient(self.config)
        
        # 3. 初始化任务追踪器
        self.job_tracker = JobTracker(self.config)
        
        # 4. 初始化 COS 客户端
        self.cos_client = COSClient(self.config)
        
        # 5. 初始化重试协调器
        self.retry_coordinator = RetryCoordinator(self.config, self.client)
        
        # 6. 注册任务处理器
        self.handlers = {
            'md': MDHandler(self.config, self.client, self.job_tracker),
            'qc': QCHandler(self.config, self.client, self.job_tracker),
            'postprocess': PostprocessHandler(self.config, self.client, self.job_tracker),
            'binding': BindingHandler(self.config, self.client, self.job_tracker),
            'redox': RedoxHandler(self.config, self.client, self.job_tracker),
            'cluster': ClusterHandler(self.config, self.client, self.job_tracker),
            'anion': AnionHandler(self.config, self.client, self.job_tracker),
        }
        
        # 7. 创建任务调度器
        self.scheduler = TaskScheduler(
            config=self.config,
            client=self.client,
            handlers=self.handlers,
            job_tracker=self.job_tracker
        )
        
        logger.info(f"Worker '{self.config.worker_name}' 已初始化 (v2.0 模块化版本)")
    
    def run(self):
        """启动 Worker 主循环"""
        logger.info("=== Polling Worker v2.0 启动 ===")
        logger.info(f"Worker Name: {self.config.worker_name}")
        logger.info(f"Max Concurrent Jobs: {self.config.max_concurrent_jobs}")
        logger.info(f"Poll Interval: {self.config.poll_interval}s")
        logger.info(f"API Base URL: {self.config.api_base_url}")
        
        # 开始轮询
        self.scheduler.start_polling()
    
    def stop(self):
        """停止 Worker"""
        logger.info("正在停止 Worker...")
        self.scheduler.stop()


def main():
    """主入口"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Molyte Polling Worker v2.0 (模块化重构版)'
    )
    parser.add_argument(
        '--config', '-c',
        default='deployment/polling_worker_config.yaml',
        help='配置文件路径'
    )
    parser.add_argument(
        '--version', '-v',
        action='version',
        version='Polling Worker v2.0'
    )
    
    args = parser.parse_args()
    
    try:
        worker = PollingWorkerV2(config_path=args.config)
        worker.run()
    except KeyboardInterrupt:
        logger.info("\n收到中断信号，正在退出...")
    except Exception as e:
        logger.error(f"Worker 启动失败: {e}", exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
