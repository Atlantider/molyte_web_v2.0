from pydantic import BaseModel, Field
from typing import Optional

# Reproduction of the bug
class UpdateBillingConfigRequest(BaseModel):
    """更新计费配置请求"""
    pricing_mode: Optional[str] = Field(None, description="计费模式: CORE_HOUR/TASK_TYPE")

# Mocking the request object
config = UpdateBillingConfigRequest(pricing_mode="TASK_TYPE")

print("Attempting to access config.global_core_hour_price...")
try:
    if config.global_core_hour_price is not None:
        print("Success!")
except AttributeError as e:
    print(f"Caught expected error: {e}")
except Exception as e:
    print(f"Caught unexpected error: {e}")
