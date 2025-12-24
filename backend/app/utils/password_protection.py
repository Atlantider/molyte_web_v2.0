"""
密码保护机制 - 防止未授权的批量密码重置
"""
import os
import json
from datetime import datetime
from typing import List, Dict
from app.core.logger import logger

# 密码重置保护文件
PROTECTION_FILE = "/opt/molyte_web_v1.0/backend/.password_protection.json"

class PasswordProtection:
    """密码保护类"""
    
    @staticmethod
    def log_password_reset(username: str, admin_username: str, reason: str = ""):
        """记录密码重置操作"""
        try:
            # 读取现有记录
            if os.path.exists(PROTECTION_FILE):
                with open(PROTECTION_FILE, 'r') as f:
                    data = json.load(f)
            else:
                data = {"password_resets": []}
            
            # 添加新记录
            reset_record = {
                "timestamp": datetime.now().isoformat(),
                "target_username": username,
                "admin_username": admin_username,
                "reason": reason,
                "ip_address": "unknown"  # 可以从请求中获取
            }
            
            data["password_resets"].append(reset_record)
            
            # 保存记录
            with open(PROTECTION_FILE, 'w') as f:
                json.dump(data, f, indent=2)
                
            logger.warning(f"PASSWORD RESET: Admin {admin_username} reset password for {username}. Reason: {reason}")
            
        except Exception as e:
            logger.error(f"Failed to log password reset: {e}")
    
    @staticmethod
    def get_reset_history() -> List[Dict]:
        """获取密码重置历史"""
        try:
            if os.path.exists(PROTECTION_FILE):
                with open(PROTECTION_FILE, 'r') as f:
                    data = json.load(f)
                return data.get("password_resets", [])
            return []
        except Exception as e:
            logger.error(f"Failed to read reset history: {e}")
            return []
    
    @staticmethod
    def check_bulk_reset_protection() -> bool:
        """检查是否允许批量重置（需要特殊授权）"""
        # 检查是否存在紧急授权文件
        emergency_file = "/opt/molyte_web_v1.0/backend/.emergency_reset_authorized"
        if os.path.exists(emergency_file):
            # 检查授权文件的时间戳（只允许1小时内有效）
            try:
                stat = os.stat(emergency_file)
                age = datetime.now().timestamp() - stat.st_mtime
                if age < 3600:  # 1小时
                    return True
                else:
                    os.remove(emergency_file)  # 删除过期的授权
            except:
                pass
        
        return False
    
    @staticmethod
    def create_emergency_authorization():
        """创建紧急授权（只能由系统管理员手动创建）"""
        emergency_file = "/opt/molyte_web_v1.0/backend/.emergency_reset_authorized"
        with open(emergency_file, 'w') as f:
            f.write(f"Emergency password reset authorized at {datetime.now().isoformat()}\n")
        logger.critical("EMERGENCY PASSWORD RESET AUTHORIZATION CREATED")

def require_password_reset_auth(func):
    """装饰器：要求密码重置授权"""
    def wrapper(*args, **kwargs):
        if not PasswordProtection.check_bulk_reset_protection():
            raise Exception("Bulk password reset requires emergency authorization")
        return func(*args, **kwargs)
    return wrapper
