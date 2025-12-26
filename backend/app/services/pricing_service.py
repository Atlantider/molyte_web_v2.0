"""
定价服务 - 任务类型差异化定价和用户折扣
"""
from typing import Optional, Tuple
from sqlalchemy.orm import Session
from sqlalchemy import func, text

from app.models.user import User, UserType
from app.core.logger import logger


class TaskTypePricingModel:
    """任务类型定价模型(临时,直到创建正式模型)"""
    pass


class UserTypeDiscountModel:
    """用户类型折扣模型(临时,直到创建正式模型)"""
    pass


class PricingService:
    """定价服务类"""
    
    # 默认定价(当数据库中没有配置时使用)
    DEFAULT_TASK_TYPE_PRICES = {
        'MD': 1.0,
        'QC': 2.0,
        'POSTPROCESS': 5.0,
        'REACTION_NETWORK': 10.0,
        'FORCEFIELD': 100.0,
    }
    
    DEFAULT_USER_TYPE_DISCOUNTS = {
        'GUEST': 1.0,      # 访客 - 标准价格
        'USER': 0.9,       # 普通用户 - 9折
        'PREMIUM': 0.7,    # 高级用户 - 7折
        'ADMIN': 0.0,      # 管理员 - 免费
    }
    
    @staticmethod
    def get_task_type_price(db: Session, task_type: str) -> float:
        """
        获取任务类型单价
        
        Args:
            db: 数据库会话
            task_type: 任务类型 (MD/QC/POSTPROCESS/REACTION_NETWORK/FORCEFIELD)
            
        Returns:
            float: 核时单价(元/核时)
        """
        try:
            # 从数据库查询
            result = db.execute(
                text("SELECT price_per_hour FROM task_type_pricing WHERE task_type = :task_type AND is_active = TRUE"),
                {"task_type": task_type}
            ).first()
            
            if result:
                return float(result[0])
        except Exception as e:
            logger.warning(f"Failed to get task type price from DB: {e}")
        
        # 使用默认值
        return PricingService.DEFAULT_TASK_TYPE_PRICES.get(task_type, 1.0)
    
    @staticmethod
    def get_user_discount_rate(db: Session, user_id: int) -> float:
        """
        获取用户折扣率
        
        优先级:
        1. 用户自定义折扣 (users.custom_discount_rate)
        2. 用户类型折扣 (user_type_discounts表)
        3. 无折扣 (1.0)
        
        Args:
            db: 数据库会话
            user_id: 用户ID
            
        Returns:
            float: 折扣率 (0.9表示9折, 1.0表示无折扣)
        """
        user = db.query(User).filter(User.id == user_id).first()
        if not user:
            return 1.0
        
        # 1. 检查用户自定义折扣
        if user.custom_discount_rate is not None:
            return float(user.custom_discount_rate)
        
        # 2. 检查用户类型折扣
        try:
            user_type_value = user.user_type.value if hasattr(user.user_type, 'value') else str(user.user_type)
            
            result = db.execute(
                text("SELECT discount_rate FROM user_type_discounts WHERE user_type = :user_type AND is_active = TRUE"),
                {"user_type": user_type_value}
            ).first()
            
            if result:
                return float(result[0])
        except Exception as e:
            logger.warning(f"Failed to get user type discount from DB: {e}")
            # 使用默认值
            user_type_value = user.user_type.value if hasattr(user.user_type, 'value') else str(user.user_type)
            return PricingService.DEFAULT_USER_TYPE_DISCOUNTS.get(user_type_value, 1.0)
        
        # 3. 无折扣
        return 1.0
    
    @staticmethod
    def calculate_final_price(db: Session, task_type: str, cpu_hours: float, user_id: int) -> float:
        """
        计算最终价格
        
        公式: 最终价格 = 任务类型单价 × 核时 × 用户折扣率
        
        Args:
            db: 数据库会话
            task_type: 任务类型
            cpu_hours: 核时数
            user_id: 用户ID
            
        Returns:
            float: 最终价格(元)
        """
        task_price = PricingService.get_task_type_price(db, task_type)
        discount_rate = PricingService.get_user_discount_rate(db, user_id)
        
        final_price = task_price * cpu_hours * discount_rate
        
        logger.info(
            f"Price calculation: task_type={task_type}, cpu_hours={cpu_hours:.2f}, "
            f"task_price={task_price:.2f}, discount={discount_rate:.2f}, final={final_price:.2f}"
        )
        
        return final_price
    
    @staticmethod
    def admin_update_task_type_price(
        db: Session,
        task_type: str,
        price: float,
        admin_id: int
    ) -> Tuple[bool, str]:
        """
        管理员更新任务类型单价
        
        Args:
            db: 数据库会话
            task_type: 任务类型
            price: 新单价
            admin_id: 管理员ID
            
        Returns:
            Tuple[bool, str]: (成功标志, 消息)
        """
        try:
            # 检查是否存在
            result = db.execute(
                text("SELECT id FROM task_type_pricing WHERE task_type = :task_type"),
                {"task_type": task_type}
            ).first()
            
            if result:
                # 更新
                db.execute(
                    text("UPDATE task_type_pricing SET price_per_hour = :price, updated_by = :admin_id, updated_at = NOW() "
                    "WHERE task_type = :task_type"),
                    {"price": price, "admin_id": admin_id, "task_type": task_type}
                )
            else:
                # 插入
                db.execute(
                    text("INSERT INTO task_type_pricing (task_type, price_per_hour, updated_by) "
                    "VALUES (:task_type, :price, :admin_id)"),
                    {"task_type": task_type, "price": price, "admin_id": admin_id}
                )
            
            db.commit()
            logger.info(f"Admin {admin_id} updated {task_type} price to {price}")
            return True, f"已更新 {task_type} 单价为 ¥{price:.2f}/核时"
            
        except Exception as e:
            db.rollback()
            logger.error(f"Failed to update task type price: {e}")
            return False, f"更新失败: {str(e)}"
    
    @staticmethod
    def admin_update_user_type_price(
        db: Session,
        user_type: str,
        core_hour_price: float,
        admin_id: int
    ) -> Tuple[bool, str]:
        """
        管理员更新用户类型核时单价
        
        Args:
            db: 数据库会话
            user_type: 用户类型
            core_hour_price: 核时单价
            admin_id: 管理员ID
            
        Returns:
            Tuple[bool, str]: (成功标志, 消息)
        """
        try:
            # 检查是否存在
            result = db.execute(
                text("SELECT id FROM user_type_prices WHERE user_type = :user_type"),
                {"user_type": user_type}
            ).first()
            
            if result:
                # 更新
                db.execute(
                    text("UPDATE user_type_prices SET core_hour_price = :price, updated_by = :admin_id, updated_at = NOW() "
                    "WHERE user_type = :user_type"),
                    {"price": core_hour_price, "admin_id": admin_id, "user_type": user_type}
                )
            else:
                # 插入  
                db.execute(
                    text("INSERT INTO user_type_prices (user_type, core_hour_price, updated_by) "
                    "VALUES (:user_type, :price, :admin_id)"),
                    {"user_type": user_type, "price": core_hour_price, "admin_id": admin_id}
                )
            
            db.commit()
            logger.info(f"Admin {admin_id} updated {user_type} core hour price to {core_hour_price}")
            return True, f"已更新 {user_type} 核时单价为 ¥{core_hour_price:.2f}/核时"
            
        except Exception as e:
            db.rollback()
            logger.error(f"Failed to update user type price: {e}")
            return False, f"更新失败: {str(e)}"
    
    @staticmethod
    def admin_set_user_custom_discount(
        db: Session,
        user_id: int,
        discount_rate: Optional[float],
        admin_id: int
    ) -> Tuple[bool, str]:
        """
        管理员设置特定用户自定义折扣
        
        Args:
            db: 数据库会话
            user_id: 用户ID
            discount_rate: 折扣率 (None表示删除自定义折扣)
            admin_id: 管理员ID
            
        Returns:
            Tuple[bool, str]: (成功标志, 消息)
        """
        user = db.query(User).filter(User.id == user_id).first()
        if not user:
            return False, f"用户 {user_id} 不存在"
        
        user.custom_discount_rate = discount_rate
        db.commit()
        
        if discount_rate is None:
            logger.info(f"Admin {admin_id} removed custom discount for user {user_id}")
            return True, f"已删除用户 {user.username} 的自定义折扣"
        else:
            discount_percent = int(discount_rate * 100)
            logger.info(f"Admin {admin_id} set custom discount {discount_rate} for user {user_id}")
            return True, f"已设置用户 {user.username} 自定义折扣为 {discount_percent}%"
