"""Add qc_engine and retry fields to QCJob

Revision ID: add_qc_engine_retry
Revises: add_gaussian_permission
Create Date: 2025-12-24

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'add_qc_engine_retry'
down_revision = 'add_gaussian_permission'
branch_labels = None
depends_on = None


def upgrade():
    # 添加qc_engine字段，默认gaussian
    op.add_column('qc_jobs', sa.Column('qc_engine', sa.String(), nullable=False, server_default='gaussian'))
    
    # 添加重试相关字段
    op.add_column('qc_jobs', sa.Column('retry_count', sa.Integer(), nullable=False, server_default='0'))
    op.add_column('qc_jobs', sa.Column('max_retries', sa.Integer(), nullable=False, server_default='3'))


def downgrade():
    # 移除添加的字段
    op.drop_column('qc_jobs', 'max_retries')
    op.drop_column('qc_jobs', 'retry_count')
    op.drop_column('qc_jobs', 'qc_engine')
