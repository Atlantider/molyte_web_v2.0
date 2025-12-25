"""Add can_use_gaussian permission to users

Revision ID: add_gaussian_permission
Revises: (previous_revision)
Create Date: 2025-12-24

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'add_gaussian_permission'
down_revision = None  # 需要手动填入上一个迁移的ID
branch_labels = None
depends_on = None


def upgrade():
    # 添加can_use_gaussian字段
    # 默认False,需要管理员手动授权
    op.add_column('users', sa.Column('can_use_gaussian', sa.Boolean(), nullable=False, server_default='false'))
    
    # 可选: 给所有ADMIN用户自动授权Gaussian使用权限
    op.execute("""
        UPDATE users 
        SET can_use_gaussian = true 
        WHERE role = 'ADMIN'
    """)


def downgrade():
    # 移除can_use_gaussian字段
    op.drop_column('users', 'can_use_gaussian')
