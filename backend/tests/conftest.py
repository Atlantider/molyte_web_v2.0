"""
Pytest 配置文件

配置测试环境、数据库、fixtures 等
"""

import pytest
import os
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, Session
from sqlalchemy.pool import StaticPool

# 设置测试环境变量
os.environ['TESTING'] = 'true'


@pytest.fixture(scope="session")
def test_db_engine():
    """创建测试数据库引擎"""
    # 尝试连接到 PostgreSQL 测试数据库
    # 如果失败，则跳过数据库测试
    try:
        from sqlalchemy import text

        # 尝试连接到 PostgreSQL
        engine = create_engine(
            "postgresql://molyte_user:molyte_password@localhost/molyte_db_test",
            echo=False
        )

        # 测试连接
        with engine.connect() as conn:
            conn.execute(text("SELECT 1"))

        # 创建所有表
        from app.database import Base
        Base.metadata.drop_all(bind=engine)
        Base.metadata.create_all(bind=engine)

        yield engine

        # 清理
        Base.metadata.drop_all(bind=engine)

    except Exception as e:
        print(f"无法连接到 PostgreSQL: {e}")
        print("跳过数据库测试")
        pytest.skip("PostgreSQL 数据库不可用")


@pytest.fixture(scope="function")
def db_session(test_db_engine):
    """创建测试数据库会话"""
    connection = test_db_engine.connect()
    transaction = connection.begin()
    session = sessionmaker(autocommit=False, autoflush=False, bind=connection)()
    
    yield session
    
    session.close()
    transaction.rollback()
    connection.close()


@pytest.fixture(scope="function")
def test_user(db_session):
    """创建测试用户"""
    from app.models import User
    
    user = User(
        username="testuser",
        email="test@example.com",
        password_hash="hashed_password",
        is_active=True,
        is_admin=False,
        total_cpu_hours=100.0,
        balance_cpu_hours=100.0,
        frozen_cpu_hours=0.0,
        used_cpu_hours=0.0
    )
    db_session.add(user)
    db_session.commit()
    db_session.refresh(user)
    return user


@pytest.fixture(scope="function")
def test_admin(db_session):
    """创建测试管理员"""
    from app.models import User
    
    admin = User(
        username="admin",
        email="admin@example.com",
        password_hash="hashed_password",
        is_active=True,
        is_admin=True,
        total_cpu_hours=0.0,
        balance_cpu_hours=0.0,
        frozen_cpu_hours=0.0,
        used_cpu_hours=0.0
    )
    db_session.add(admin)
    db_session.commit()
    db_session.refresh(admin)
    return admin



@pytest.fixture(scope="function")
def test_master_account(db_session, test_user, test_organization):
    """创建测试主账号"""
    from app.models import MasterAccount
    from datetime import datetime
    
    master = MasterAccount(
        user_id=test_user.id,
        organization_id=test_organization.id,
        total_cpu_hours=500.0,
        balance_cpu_hours=500.0,
        frozen_cpu_hours=0.0,
        used_cpu_hours=0.0,
        max_sub_accounts=10,
        current_sub_accounts=0,
        is_active=True,
        created_at=datetime.utcnow(),
        updated_at=datetime.utcnow()
    )
    db_session.add(master)
    db_session.commit()
    db_session.refresh(master)
    return master


@pytest.fixture(scope="function")
def test_sub_account(db_session, test_master_account):
    """创建测试子账号"""
    from app.models import User, SubAccount
    from datetime import datetime
    
    # 创建子账号用户
    sub_user = User(
        username="subuser",
        email="sub@example.com",
        password_hash="hashed_password",
        is_active=True,
        is_admin=False,
        total_cpu_hours=0.0,
        balance_cpu_hours=0.0,
        frozen_cpu_hours=0.0,
        used_cpu_hours=0.0
    )
    db_session.add(sub_user)
    db_session.commit()
    db_session.refresh(sub_user)
    
    # 创建子账号
    sub_account = SubAccount(
        user_id=sub_user.id,
        master_account_id=test_master_account.id,
        personal_quota=100.0,
        personal_used=0.0,
        is_active=True,
        created_at=datetime.utcnow(),
        updated_at=datetime.utcnow()
    )
    db_session.add(sub_account)
    db_session.commit()
    db_session.refresh(sub_account)
    return sub_account

