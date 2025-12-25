"""
Add labels column to electrolyte_systems table

Revision ID: add_electrolyte_labels
Create Date: 2025-12-25
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql

# revision identifiers
revision = 'add_electrolyte_labels'
down_revision = None  # 请根据实际情况设置前一个migration的revision ID
branch_labels = None
depends_on = None


def upgrade():
    """添加 labels 列到 electrolyte_systems 表"""
    op.add_column(
        'electrolyte_systems',
        sa.Column('labels', postgresql.JSONB(), nullable=True, server_default='{}')
    )
    
    # 为 labels 添加 GIN 索引以支持 JSON 查询
    op.create_index(
        'idx_electrolyte_labels',
        'electrolyte_systems',
        ['labels'],
        postgresql_using='gin'
    )


def downgrade():
    """删除 labels 列"""
    op.drop_index('idx_electrolyte_labels', table_name='electrolyte_systems')
    op.drop_column('electrolyte_systems', 'labels')
