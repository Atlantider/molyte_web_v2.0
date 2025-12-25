"""
Quota System Unification Migration

添加enable_resource_limits字段并更新字段默认值

Revision ID: quota_unification_v1
Created: 2025-12-25
"""

def upgrade():
    """升级数据库"""
    import sqlalchemy as sa
    from alembic import op
    
    # 1. 添加资源保护开关字段
    op.add_column('users',
        sa.Column('enable_resource_limits', sa.Boolean(),
                  nullable=False, server_default='true'))
    
    # 2. 更新管理员用户的配置
    # 管理员不受资源限制
    op.execute("""
        UPDATE users
        SET enable_resource_limits = false,
            concurrent_job_limit = 0,
            daily_job_limit = 0,
            storage_quota_gb = 0
        WHERE role = 'ADMIN';
    """)
    
    # 3. 更新高级用户的配置
    # PREMIUM用户获得更高限制
    op.execute("""
        UPDATE users
        SET concurrent_job_limit = 10,
            daily_job_limit = 50,
            storage_quota_gb = 200.0
        WHERE role = 'PREMIUM';
    """)
    
    # 4. 更新普通用户的配置
    # USER角色使用标准限制
    op.execute("""
        UPDATE users
        SET concurrent_job_limit = 5,
            daily_job_limit = 20,
            storage_quota_gb = 50.0
        WHERE role = 'USER' AND concurrent_job_limit < 5;
    """)
    
    #  5. GUEST用户限制更严格
    op.execute("""
        UPDATE users
        SET concurrent_job_limit = 2,
            daily_job_limit = 5,
            storage_quota_gb = 10.0
        WHERE role = 'GUEST';
    """)
    
    print("✅ Quota system unification migration completed")
    print("  - ADMIN: 无限制")
    print("  - PREMIUM: 10并发/50每日/200GB")
    print("  - USER: 5并发/20每日/50GB")
    print("  - GUEST: 2并发/5每日/10GB")


def downgrade():
    """回滚数据库"""
    from alembic import op
    
    # 移除新增字段
    op.drop_column('users', 'enable_resource_limits')
    
    print("⚠️  Quota system unification migration rolled back")
