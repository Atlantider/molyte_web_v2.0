"""
支付服务
支持微信支付、支付宝等支付方式
"""
import hashlib
import hmac
import json
import time
import uuid
import base64
from datetime import datetime, timedelta
from typing import Optional, Dict, Tuple, Any
from decimal import Decimal
import httpx
from sqlalchemy.orm import Session

from app.config import settings
from app.core.logger import logger
from app.models.billing import RechargeOrder, PaymentStatus, PaymentMethod


class PaymentProvider:
    """支付提供商基类"""
    
    async def create_payment(self, order: RechargeOrder) -> Tuple[bool, Dict[str, Any]]:
        """创建支付订单，返回支付参数"""
        raise NotImplementedError
    
    async def query_payment(self, order_no: str) -> Tuple[bool, Dict[str, Any]]:
        """查询支付状态"""
        raise NotImplementedError
    
    def verify_callback(self, data: Dict) -> Tuple[bool, str]:
        """验证回调签名"""
        raise NotImplementedError


class WechatPayProvider(PaymentProvider):
    """微信支付 - Native 扫码支付"""
    
    def __init__(self):
        self.app_id = getattr(settings, 'WECHAT_APP_ID', '')
        self.mch_id = getattr(settings, 'WECHAT_MCH_ID', '')  # 商户号
        self.api_key = getattr(settings, 'WECHAT_API_KEY', '')  # API密钥
        self.api_v3_key = getattr(settings, 'WECHAT_API_V3_KEY', '')  # APIv3密钥
        self.serial_no = getattr(settings, 'WECHAT_SERIAL_NO', '')  # 证书序列号
        self.private_key = getattr(settings, 'WECHAT_PRIVATE_KEY', '')  # 商户私钥
        self.notify_url = getattr(settings, 'WECHAT_NOTIFY_URL', '')
    
    async def create_payment(self, order: RechargeOrder) -> Tuple[bool, Dict[str, Any]]:
        """创建微信支付订单（Native支付-扫码）"""
        try:
            url = "https://api.mch.weixin.qq.com/v3/pay/transactions/native"
            
            # 构建请求参数
            params = {
                "appid": self.app_id,
                "mchid": self.mch_id,
                "description": f"Molyte配额充值 - {order.cpu_hours}核时",
                "out_trade_no": order.order_no,
                "notify_url": self.notify_url,
                "amount": {
                    "total": int(order.amount * 100),  # 单位：分
                    "currency": "CNY"
                },
                "time_expire": (datetime.now() + timedelta(hours=2)).strftime("%Y-%m-%dT%H:%M:%S+08:00"),
            }
            
            # 生成签名 (实际需要使用商户证书签名)
            # 这里简化处理，生产环境需要完整实现
            headers = self._build_auth_header("POST", "/v3/pay/transactions/native", json.dumps(params))
            
            async with httpx.AsyncClient() as client:
                response = await client.post(url, json=params, headers=headers, timeout=30)
                result = response.json()
                
                if response.status_code == 200 and "code_url" in result:
                    return True, {
                        "qr_code": result["code_url"],  # 二维码链接
                        "order_no": order.order_no,
                        "expire_time": 7200,  # 2小时有效
                    }
                else:
                    return False, {"error": result.get("message", "创建支付失败")}
                    
        except Exception as e:
            logger.error(f"微信支付创建失败: {e}")
            return False, {"error": str(e)}
    
    def _build_auth_header(self, method: str, url: str, body: str) -> Dict[str, str]:
        """构建认证头（简化版，生产环境需完整实现）"""
        timestamp = str(int(time.time()))
        nonce = uuid.uuid4().hex
        
        # 实际需要使用私钥签名
        # signature = self._sign(method, url, timestamp, nonce, body)
        
        return {
            "Content-Type": "application/json",
            "Accept": "application/json",
            # "Authorization": f'WECHATPAY2-SHA256-RSA2048 mchid="{self.mch_id}",nonce_str="{nonce}",timestamp="{timestamp}",serial_no="{self.serial_no}",signature="{signature}"'
        }
    
    async def query_payment(self, order_no: str) -> Tuple[bool, Dict[str, Any]]:
        """查询微信支付状态"""
        try:
            url = f"https://api.mch.weixin.qq.com/v3/pay/transactions/out-trade-no/{order_no}?mchid={self.mch_id}"
            headers = self._build_auth_header("GET", f"/v3/pay/transactions/out-trade-no/{order_no}?mchid={self.mch_id}", "")
            
            async with httpx.AsyncClient() as client:
                response = await client.get(url, headers=headers, timeout=30)
                result = response.json()
                
                if response.status_code == 200:
                    trade_state = result.get("trade_state")
                    return True, {
                        "paid": trade_state == "SUCCESS",
                        "trade_state": trade_state,
                        "transaction_id": result.get("transaction_id"),
                    }
                return False, {"error": result.get("message", "查询失败")}
                
        except Exception as e:
            logger.error(f"微信支付查询失败: {e}")
            return False, {"error": str(e)}
    
    def verify_callback(self, data: Dict) -> Tuple[bool, str]:
        """验证微信支付回调签名"""
        # 实际需要验证签名
        return True, data.get("out_trade_no", "")


class AlipayProvider(PaymentProvider):
    """支付宝支付"""
    
    def __init__(self):
        self.app_id = getattr(settings, 'ALIPAY_APP_ID', '')
        self.private_key = getattr(settings, 'ALIPAY_PRIVATE_KEY', '')
        self.alipay_public_key = getattr(settings, 'ALIPAY_PUBLIC_KEY', '')
        self.notify_url = getattr(settings, 'ALIPAY_NOTIFY_URL', '')
        self.return_url = getattr(settings, 'ALIPAY_RETURN_URL', '')
        self.gateway = "https://openapi.alipay.com/gateway.do"
    
    async def create_payment(self, order: RechargeOrder) -> Tuple[bool, Dict[str, Any]]:
        """创建支付宝订单（电脑网站支付）"""
        try:
            # 构建业务参数
            biz_content = {
                "out_trade_no": order.order_no,
                "total_amount": str(order.amount),
                "subject": f"Molyte配额充值 - {order.cpu_hours}核时",
                "product_code": "FAST_INSTANT_TRADE_PAY",
            }

            # 公共参数
            params = {
                "app_id": self.app_id,
                "method": "alipay.trade.page.pay",
                "charset": "utf-8",
                "sign_type": "RSA2",
                "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "version": "1.0",
                "notify_url": self.notify_url,
                "return_url": self.return_url,
                "biz_content": json.dumps(biz_content),
            }

            # 生成签名（实际需要使用私钥签名）
            # params["sign"] = self._sign(params)

            # 构建支付链接
            query_string = "&".join([f"{k}={v}" for k, v in sorted(params.items())])
            pay_url = f"{self.gateway}?{query_string}"

            return True, {
                "pay_url": pay_url,
                "order_no": order.order_no,
                "expire_time": 7200,
            }

        except Exception as e:
            logger.error(f"支付宝支付创建失败: {e}")
            return False, {"error": str(e)}

    async def query_payment(self, order_no: str) -> Tuple[bool, Dict[str, Any]]:
        """查询支付宝支付状态"""
        try:
            biz_content = {"out_trade_no": order_no}
            params = {
                "app_id": self.app_id,
                "method": "alipay.trade.query",
                "charset": "utf-8",
                "sign_type": "RSA2",
                "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "version": "1.0",
                "biz_content": json.dumps(biz_content),
            }

            async with httpx.AsyncClient() as client:
                response = await client.post(self.gateway, data=params, timeout=30)
                result = response.json()

                trade_status = result.get("alipay_trade_query_response", {}).get("trade_status")
                return True, {
                    "paid": trade_status == "TRADE_SUCCESS",
                    "trade_status": trade_status,
                }

        except Exception as e:
            logger.error(f"支付宝查询失败: {e}")
            return False, {"error": str(e)}

    def verify_callback(self, data: Dict) -> Tuple[bool, str]:
        """验证支付宝回调签名"""
        # 实际需要验证签名
        return True, data.get("out_trade_no", "")


class MockPaymentProvider(PaymentProvider):
    """模拟支付（开发测试用）"""

    async def create_payment(self, order: RechargeOrder) -> Tuple[bool, Dict[str, Any]]:
        """创建模拟支付"""
        return True, {
            "mock": True,
            "order_no": order.order_no,
            "amount": order.amount,
            "message": "这是模拟支付，点击确认支付按钮完成支付",
        }

    async def query_payment(self, order_no: str) -> Tuple[bool, Dict[str, Any]]:
        """模拟查询"""
        return True, {"paid": False, "mock": True}

    def verify_callback(self, data: Dict) -> Tuple[bool, str]:
        return True, data.get("order_no", "")


class PaymentService:
    """支付服务管理器"""

    def __init__(self):
        self.providers = {
            "wechat": WechatPayProvider(),
            "alipay": AlipayProvider(),
            "simulated": MockPaymentProvider(),
        }

    def get_provider(self, method: str) -> PaymentProvider:
        """获取支付提供商"""
        return self.providers.get(method, self.providers["simulated"])

    async def create_payment(self, order: RechargeOrder) -> Tuple[bool, Dict[str, Any]]:
        """创建支付"""
        provider = self.get_provider(order.payment_method)
        return await provider.create_payment(order)

    async def query_payment(self, order: RechargeOrder) -> Tuple[bool, Dict[str, Any]]:
        """查询支付状态"""
        provider = self.get_provider(order.payment_method)
        return await provider.query_payment(order.order_no)

    def verify_callback(self, method: str, data: Dict) -> Tuple[bool, str]:
        """验证支付回调"""
        provider = self.get_provider(method)
        return provider.verify_callback(data)


# 单例实例
payment_service = PaymentService()

