import sys
import os

# Add backend directory to sys.path
sys.path.append('/opt/molyte_web_v2.0/backend')

try:
    from app.api.v1.pricing import UpdateBillingConfigRequest, BillingConfigSchema
    
    # Check if global_core_hour_price field exists
    if 'global_core_hour_price' not in UpdateBillingConfigRequest.model_fields:
        print("FAIL: global_core_hour_price missing from UpdateBillingConfigRequest")
        sys.exit(1)
        
    if 'global_core_hour_price' not in BillingConfigSchema.model_fields:
        print("FAIL: global_core_hour_price missing from BillingConfigSchema")
        sys.exit(1)
        
    print("SUCCESS: Fields present in Pydantic models.")
    
except Exception as e:
    print(f"FAIL: Error importing or checking models: {e}")
    sys.exit(1)
