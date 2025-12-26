"""
Database migration for Reaction Network module
Create reaction_network_jobs, molecules, and reactions tables

Revision ID: add_reaction_network_tables
Revises: previous_migration_id
Create Date: 2025-12-25
"""

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql

# revision identifiers, used by Alembic
revision = 'add_reaction_network_tables'
down_revision = None  # Update this with the actual previous migration ID
branch_labels = None
depends_on = None


def upgrade():
    # Create reaction_network_jobs table
    op.create_table(
        'reaction_network_jobs',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('user_id', sa.Integer(), nullable=False),
        sa.Column('job_name', sa.String(length=255), nullable=False),
        sa.Column('description', sa.Text(), nullable=True),
        sa.Column('status', sa.Enum('CREATED', 'QUEUED', 'RUNNING', 'POSTPROCESSING', 'COMPLETED', 'FAILED', 'CANCELLED', name='reactionnetworkjobstatus'), nullable=False),
        sa.Column('progress', sa.Float(), nullable=True),
        sa.Column('error_message', sa.Text(), nullable=True),
        sa.Column('initial_smiles', postgresql.JSONB(astext_type=sa.Text()), nullable=False),
        sa.Column('temperature', sa.Float(), nullable=True),
        sa.Column('electrode_type', sa.Enum('anode', 'cathode', name='electrodetype'), nullable=True),
        sa.Column('voltage', sa.Float(), nullable=True),
        sa.Column('max_generations', sa.Integer(), nullable=True),
        sa.Column('max_species', sa.Integer(), nullable=True),
        sa.Column('energy_cutoff', sa.Float(), nullable=True),
        sa.Column('config', postgresql.JSONB(astext_type=sa.Text()), nullable=True),
        sa.Column('slurm_job_id', sa.String(length=64), nullable=True),
        sa.Column('work_dir', sa.String(length=512), nullable=True),
        sa.Column('slurm_partition', sa.String(length=64), nullable=True),
        sa.Column('slurm_cpus', sa.Integer(), nullable=True),
        sa.Column('slurm_time', sa.Integer(), nullable=True),
        sa.Column('actual_cpu_hours', sa.Float(), nullable=True),
        sa.Column('num_molecules', sa.Integer(), nullable=True),
        sa.Column('num_reactions', sa.Integer(), nullable=True),
        sa.Column('max_generation_reached', sa.Integer(), nullable=True),
        sa.Column('network_json_path', sa.String(length=512), nullable=True),
        sa.Column('visualization_png_path', sa.String(length=512), nullable=True),
        sa.Column('visualization_html_path', sa.String(length=512), nullable=True),
        sa.Column('created_at', sa.DateTime(), nullable=False),
        sa.Column('updated_at', sa.DateTime(), nullable=True),
        sa.Column('started_at', sa.DateTime(), nullable=True),
        sa.Column('finished_at', sa.DateTime(), nullable=True),
        sa.Column('is_deleted', sa.Boolean(), nullable=True),
        sa.Column('deleted_at', sa.DateTime(), nullable=True),
        sa.Column('deleted_by', sa.Integer(), nullable=True),
        sa.Column('delete_reason', sa.String(length=512), nullable=True),
        sa.ForeignKeyConstraint(['user_id'], ['users.id'], ),
        sa.ForeignKeyConstraint(['deleted_by'], ['users.id'], ),
        sa.PrimaryKeyConstraint('id')
    )
    op.create_index(op.f('ix_reaction_network_jobs_id'), 'reaction_network_jobs', ['id'], unique=False)
    op.create_index(op.f('ix_reaction_network_jobs_user_id'), 'reaction_network_jobs', ['user_id'], unique=False)
    op.create_index(op.f('ix_reaction_network_jobs_status'), 'reaction_network_jobs', ['status'], unique=False)
    op.create_index(op.f('ix_reaction_network_jobs_slurm_job_id'), 'reaction_network_jobs', ['slurm_job_id'], unique=False)

    # Create reaction_network_molecules table
    op.create_table(
        'reaction_network_molecules',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('job_id', sa.Integer(), nullable=False),
        sa.Column('name', sa.String(length=255), nullable=False),
        sa.Column('smiles', sa.String(length=512), nullable=False),
        sa.Column('generation', sa.Integer(), nullable=True),
        sa.Column('energy_kcal', sa.Float(), nullable=True),
        sa.Column('molecular_weight', sa.Float(), nullable=True),
        sa.Column('num_atoms', sa.Integer(), nullable=True),
        sa.Column('num_heavy_atoms', sa.Integer(), nullable=True),
        sa.Column('formal_charge', sa.Integer(), nullable=True),
        sa.Column('num_rings', sa.Integer(), nullable=True),
        sa.Column('morgan_fingerprint', sa.String(length=512), nullable=True),
        sa.Column('properties', postgresql.JSONB(astext_type=sa.Text()), nullable=True),
        sa.Column('created_at', sa.DateTime(), nullable=True),
        sa.ForeignKeyConstraint(['job_id'], ['reaction_network_jobs.id'], ondelete='CASCADE'),
        sa.PrimaryKeyConstraint('id')
    )
    op.create_index(op.f('ix_reaction_network_molecules_id'), 'reaction_network_molecules', ['id'], unique=False)
    op.create_index(op.f('ix_reaction_network_molecules_job_id'), 'reaction_network_molecules', ['job_id'], unique=False)
    op.create_index(op.f('ix_reaction_network_molecules_smiles'), 'reaction_network_molecules', ['smiles'], unique=False)
    op.create_index(op.f('ix_reaction_network_molecules_generation'), 'reaction_network_molecules', ['generation'], unique=False)

    # Create reaction_network_reactions table
    op.create_table(
        'reaction_network_reactions',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('job_id', sa.Integer(), nullable=False),
        sa.Column('reactant_smiles', postgresql.JSONB(astext_type=sa.Text()), nullable=False),
        sa.Column('product_smiles', postgresql.JSONB(astext_type=sa.Text()), nullable=False),
        sa.Column('operator_name', sa.String(length=128), nullable=True),
        sa.Column('reaction_type', sa.String(length=128), nullable=True),
        sa.Column('reaction_energy', sa.Float(), nullable=True),
        sa.Column('activation_energy', sa.Float(), nullable=True),
        sa.Column('driving_force', sa.String(length=128), nullable=True),
        sa.Column('is_reversible', sa.Boolean(), nullable=True),
        sa.Column('equilibrium_constant', sa.Float(), nullable=True),
        sa.Column('metadata', postgresql.JSONB(astext_type=sa.Text()), nullable=True),
        sa.Column('created_at', sa.DateTime(), nullable=True),
        sa.ForeignKeyConstraint(['job_id'], ['reaction_network_jobs.id'], ondelete='CASCADE'),
        sa.PrimaryKeyConstraint('id')
    )
    op.create_index(op.f('ix_reaction_network_reactions_id'), 'reaction_network_reactions', ['id'], unique=False)
    op.create_index(op.f('ix_reaction_network_reactions_job_id'), 'reaction_network_reactions', ['job_id'], unique=False)
    op.create_index(op.f('ix_reaction_network_reactions_operator_name'), 'reaction_network_reactions', ['operator_name'], unique=False)
    op.create_index(op.f('ix_reaction_network_reactions_reaction_type'), 'reaction_network_reactions', ['reaction_type'], unique=False)


def downgrade():
    op.drop_index(op.f('ix_reaction_network_reactions_reaction_type'), table_name='reaction_network_reactions')
    op.drop_index(op.f('ix_reaction_network_reactions_operator_name'), table_name='reaction_network_reactions')
    op.drop_index(op.f('ix_reaction_network_reactions_job_id'), table_name='reaction_network_reactions')
    op.drop_index(op.f('ix_reaction_network_reactions_id'), table_name='reaction_network_reactions')
    op.drop_table('reaction_network_reactions')
    
    op.drop_index(op.f('ix_reaction_network_molecules_generation'), table_name='reaction_network_molecules')
    op.drop_index(op.f('ix_reaction_network_molecules_smiles'), table_name='reaction_network_molecules')
    op.drop_index(op.f('ix_reaction_network_molecules_job_id'), table_name='reaction_network_molecules')
    op.drop_index(op.f('ix_reaction_network_molecules_id'), table_name='reaction_network_molecules')
    op.drop_table('reaction_network_molecules')
    
    op.drop_index(op.f('ix_reaction_network_jobs_slurm_job_id'), table_name='reaction_network_jobs')
    op.drop_index(op.f('ix_reaction_network_jobs_status'), table_name='reaction_network_jobs')
    op.drop_index(op.f('ix_reaction_network_jobs_user_id'), table_name='reaction_network_jobs')
    op.drop_index(op.f('ix_reaction_network_jobs_id'), table_name='reaction_network_jobs')
    op.drop_table('reaction_network_jobs')
    
    op.execute('DROP TYPE reactionnetworkjobstatus')
    op.execute('DROP TYPE electrodetype')
