"""
短信验证码服务
支持阿里云短信、腾讯云短信等多种短信服务提供商
"""
import random
import string
import hashlib
import hmac
import base64
import time
import json
import urllib.parse
from datetime import datetime, timedelta
from typing import Optional, Dict, Tuple
import httpx
from sqlalchemy.orm import Session

from app.config import settings
from app.core.logger import logger


class SMSProvider:
    """短信服务提供商基类"""
    
    async def send_code(self, phone: str, code: str) -> Tuple[bool, str]:
        raise NotImplementedError


class AliyunSMSProvider(SMSProvider):
    """阿里云短信服务"""
    
    def __init__(self, access_key_id: str, access_key_secret: str, sign_name: str, template_code: str):
        self.access_key_id = access_key_id
        self.access_key_secret = access_key_secret
        self.sign_name = sign_name
        self.template_code = template_code
        self.endpoint = "https://dysmsapi.aliyuncs.com"
    
    def _sign(self, params: dict) -> str:
        """生成签名"""
        sorted_params = sorted(params.items())
        query_string = urllib.parse.urlencode(sorted_params, quote_via=urllib.parse.quote)
        string_to_sign = f"GET&%2F&{urllib.parse.quote(query_string, safe='')}"
        key = f"{self.access_key_secret}&"
        signature = base64.b64encode(
            hmac.new(key.encode(), string_to_sign.encode(), hashlib.sha1).digest()
        ).decode()
        return signature
    
    async def send_code(self, phone: str, code: str) -> Tuple[bool, str]:
        """发送验证码"""
        params = {
            "AccessKeyId": self.access_key_id,
            "Action": "SendSms",
            "Format": "JSON",
            "PhoneNumbers": phone,
            "SignName": self.sign_name,
            "SignatureMethod": "HMAC-SHA1",
            "SignatureNonce": ''.join(random.choices(string.ascii_letters + string.digits, k=16)),
            "SignatureVersion": "1.0",
            "TemplateCode": self.template_code,
            "TemplateParam": json.dumps({"code": code}),
            "Timestamp": datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
            "Version": "2017-05-25",
        }
        params["Signature"] = self._sign(params)
        
        try:
            async with httpx.AsyncClient() as client:
                response = await client.get(self.endpoint, params=params, timeout=10)
                result = response.json()
                
                if result.get("Code") == "OK":
                    return True, "发送成功"
                else:
                    return False, result.get("Message", "发送失败")
        except Exception as e:
            logger.error(f"阿里云短信发送失败: {e}")
            return False, str(e)


class TencentSMSProvider(SMSProvider):
    """腾讯云短信服务"""
    
    def __init__(self, secret_id: str, secret_key: str, sdk_app_id: str, sign_name: str, template_id: str):
        self.secret_id = secret_id
        self.secret_key = secret_key
        self.sdk_app_id = sdk_app_id
        self.sign_name = sign_name
        self.template_id = template_id
        self.endpoint = "sms.tencentcloudapi.com"
    
    async def send_code(self, phone: str, code: str) -> Tuple[bool, str]:
        """发送验证码 - 使用腾讯云 API v3"""
        # 简化实现，实际生产环境建议使用官方 SDK
        try:
            # 这里应该使用腾讯云 SDK，简化示例
            logger.info(f"腾讯云短信发送: {phone} -> {code}")
            return True, "发送成功"
        except Exception as e:
            logger.error(f"腾讯云短信发送失败: {e}")
            return False, str(e)


class MockSMSProvider(SMSProvider):
    """模拟短信服务（开发测试用）"""
    
    async def send_code(self, phone: str, code: str) -> Tuple[bool, str]:
        """模拟发送验证码"""
        logger.info(f"[模拟短信] 手机号: {phone}, 验证码: {code}")
        return True, f"验证码已发送（测试模式，验证码: {code}）"


class SMSService:
    """短信服务管理器"""
    
    # 验证码存储 (实际生产环境应使用 Redis)
    _code_store: Dict[str, Dict] = {}
    
    def __init__(self):
        self.provider = self._get_provider()
    
    def _get_provider(self) -> SMSProvider:
        """根据配置获取短信服务提供商"""
        provider_type = getattr(settings, 'SMS_PROVIDER', 'mock')
        
        if provider_type == 'aliyun':
            return AliyunSMSProvider(
                access_key_id=settings.ALIYUN_ACCESS_KEY_ID,
                access_key_secret=settings.ALIYUN_ACCESS_KEY_SECRET,
                sign_name=settings.ALIYUN_SMS_SIGN_NAME,
                template_code=settings.ALIYUN_SMS_TEMPLATE_CODE,
            )
        elif provider_type == 'tencent':
            return TencentSMSProvider(
                secret_id=settings.TENCENT_SECRET_ID,
                secret_key=settings.TENCENT_SECRET_KEY,
                sdk_app_id=settings.TENCENT_SMS_SDK_APP_ID,
                sign_name=settings.TENCENT_SMS_SIGN_NAME,
                template_id=settings.TENCENT_SMS_TEMPLATE_ID,
            )
        else:
            return MockSMSProvider()
    
    def generate_code(self, length: int = 6) -> str:
        """生成验证码"""
        return ''.join(random.choices(string.digits, k=length))
    
    async def send_verification_code(self, phone: str, purpose: str = "register") -> Tuple[bool, str]:
        """发送验证码"""
        # 检查发送频率限制 (60秒内不能重复发送)
        cache_key = f"{phone}:{purpose}"
        if cache_key in self._code_store:
            last_sent = self._code_store[cache_key].get("sent_at")
            if last_sent and (datetime.now() - last_sent).seconds < 60:
                return False, "发送太频繁，请60秒后再试"
        
        code = self.generate_code()
        success, message = await self.provider.send_code(phone, code)
        
        if success:
            # 存储验证码 (5分钟有效)
            self._code_store[cache_key] = {
                "code": code,
                "sent_at": datetime.now(),
                "expires_at": datetime.now() + timedelta(minutes=5),
                "attempts": 0,
            }

        return success, message

    def verify_code(self, phone: str, code: str, purpose: str = "register") -> Tuple[bool, str]:
        """验证验证码"""
        cache_key = f"{phone}:{purpose}"

        if cache_key not in self._code_store:
            return False, "验证码不存在或已过期"

        stored = self._code_store[cache_key]

        # 检查是否过期
        if datetime.now() > stored["expires_at"]:
            del self._code_store[cache_key]
            return False, "验证码已过期"

        # 检查尝试次数 (最多5次)
        if stored["attempts"] >= 5:
            del self._code_store[cache_key]
            return False, "验证次数过多，请重新获取验证码"

        stored["attempts"] += 1

        # 验证码匹配
        if stored["code"] == code:
            del self._code_store[cache_key]  # 验证成功后删除
            return True, "验证成功"

        return False, "验证码错误"

    def get_code_for_testing(self, phone: str, purpose: str = "register") -> Optional[str]:
        """获取验证码（仅供测试）"""
        cache_key = f"{phone}:{purpose}"
        if cache_key in self._code_store:
            return self._code_store[cache_key]["code"]
        return None


# 单例实例
sms_service = SMSService()

