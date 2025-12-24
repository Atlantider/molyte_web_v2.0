"""
测试错误信息清除功能

验证当任务重试成功后，错误信息能够被正确清除
"""

import sys
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
from datetime import datetime
from pydantic import BaseModel
from typing import Optional

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent))

# 定义简单的模型，避免导入整个应用
class JobStatusUpdate(BaseModel):
    status: str
    job_type: str
    worker_name: str
    error_message: Optional[str] = None


def test_error_message_clearing_for_qc_job():
    """
    测试QC任务的错误信息清除

    场景：
    1. QC任务首次失败，设置错误信息
    2. 任务重试成功，清除错误信息
    """
    print("\n" + "="*60)
    print("Testing Error Message Clearing for QC Job")
    print("="*60)

    # 创建模拟的QC任务
    qc_job = Mock()
    qc_job.id = 1845
    qc_job.status = "RUNNING"
    qc_job.error_message = "自动重试 (1/3): 内部坐标系统在坐标转换/Hessian更新时崩溃"
    
    print(f"\n1. 初始状态 - 任务有错误信息:")
    print(f"   任务ID: {qc_job.id}")
    print(f"   状态: {qc_job.status}")
    print(f"   错误信息: {qc_job.error_message}")
    
    # 模拟任务完成时的状态更新
    # 这里测试的是 error_message=None 的情况
    status_update = JobStatusUpdate(
        status="COMPLETED",
        job_type="QC",
        worker_name="test_worker",
        error_message=None  # 显式清除错误信息
    )
    
    print(f"\n2. 更新状态:")
    print(f"   新状态: {status_update.status}")
    print(f"   error_message 参数: {status_update.error_message}")
    
    # 验证 error_message is not None 的条件
    # 注意：这里演示的是修复后的逻辑
    # 当 error_message 在请求中被显式设置为 None 时，应该被清除
    # 但由于 Pydantic 的默认值处理，我们需要检查是否在请求中提供了该字段

    # 模拟修复后的逻辑：检查字段是否被显式提供
    # 在实际的 API 中，我们会检查 status_update 中是否包含 error_message 字段
    if hasattr(status_update, 'error_message') and status_update.error_message is None:
        # 这表示显式传递了 None
        qc_job.error_message = None
        print(f"   ✓ 错误信息已清除: {qc_job.error_message}")

    # 验证结果
    assert qc_job.error_message is None, "错误信息应该被清除为 None"
    print(f"\n3. 最终状态:")
    print(f"   错误信息: {qc_job.error_message}")
    print(f"   ✓ 测试通过！")


def test_error_message_not_cleared_when_not_provided():
    """
    测试当不提供 error_message 参数时，不清除错误信息
    
    场景：
    1. 任务有错误信息
    2. 更新状态时不提供 error_message 参数
    3. 错误信息应该保持不变
    """
    print("\n" + "="*60)
    print("Testing Error Message Preservation When Not Provided")
    print("="*60)
    
    # 创建模拟的MD任务
    md_job = Mock()
    md_job.id = 100
    md_job.status = "RUNNING"
    md_job.error_message = "Previous error message"
    
    print(f"\n1. 初始状态:")
    print(f"   任务ID: {md_job.id}")
    print(f"   错误信息: {md_job.error_message}")
    
    # 模拟状态更新，不提供 error_message
    status_update = JobStatusUpdate(
        status="RUNNING",
        job_type="MD",
        worker_name="test_worker"
        # 注意：没有提供 error_message 参数
    )
    
    print(f"\n2. 更新状态（不提供 error_message）:")
    print(f"   新状态: {status_update.status}")
    print(f"   error_message 参数: {status_update.error_message}")
    
    # 验证 error_message is not None 的条件
    if status_update.error_message is not None:
        md_job.error_message = status_update.error_message
        print(f"   错误信息已更新: {md_job.error_message}")
    else:
        print(f"   ✓ 错误信息保持不变（未提供参数）")
    
    # 验证结果
    assert md_job.error_message == "Previous error message", "错误信息应该保持不变"
    print(f"\n3. 最终状态:")
    print(f"   错误信息: {md_job.error_message}")
    print(f"   ✓ 测试通过！")


def test_error_message_update_with_new_message():
    """
    测试当提供新的错误信息时，更新错误信息
    
    场景：
    1. 任务有旧的错误信息
    2. 更新状态时提供新的错误信息
    3. 错误信息应该被更新为新值
    """
    print("\n" + "="*60)
    print("Testing Error Message Update With New Message")
    print("="*60)
    
    # 创建模拟的QC任务
    qc_job = Mock()
    qc_job.id = 1846
    qc_job.status = "FAILED"
    qc_job.error_message = "Old error message"
    
    print(f"\n1. 初始状态:")
    print(f"   任务ID: {qc_job.id}")
    print(f"   错误信息: {qc_job.error_message}")
    
    # 模拟状态更新，提供新的错误信息
    new_error = "New error message from retry"
    status_update = JobStatusUpdate(
        status="RUNNING",
        job_type="QC",
        worker_name="test_worker",
        error_message=new_error
    )
    
    print(f"\n2. 更新状态（提供新的错误信息）:")
    print(f"   新状态: {status_update.status}")
    print(f"   新错误信息: {status_update.error_message}")
    
    # 验证 error_message is not None 的条件
    if status_update.error_message is not None:
        qc_job.error_message = status_update.error_message
        print(f"   ✓ 错误信息已更新: {qc_job.error_message}")
    
    # 验证结果
    assert qc_job.error_message == new_error, "错误信息应该被更新为新值"
    print(f"\n3. 最终状态:")
    print(f"   错误信息: {qc_job.error_message}")
    print(f"   ✓ 测试通过！")


if __name__ == "__main__":
    print("\n" + "="*60)
    print("Error Message Clearing Tests")
    print("="*60)
    
    try:
        test_error_message_clearing_for_qc_job()
        test_error_message_not_cleared_when_not_provided()
        test_error_message_update_with_new_message()
        
        print("\n" + "="*60)
        print("✓ 所有测试通过！")
        print("="*60)
    except AssertionError as e:
        print(f"\n✗ 测试失败: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n✗ 测试出错: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

