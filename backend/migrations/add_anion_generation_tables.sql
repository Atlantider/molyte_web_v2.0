-- Add anion generation and library tables
-- Created: 2025-12-08
-- Note: This migration adds support for anion force field auto-generation

-- Drop old tables if they exist (for clean migration)
DROP TABLE IF EXISTS anion_generation_jobs CASCADE;
DROP TABLE IF EXISTS anion_library CASCADE;

-- Create anion_generation_jobs table
CREATE TABLE anion_generation_jobs (
    id SERIAL PRIMARY KEY,
    job_id VARCHAR(36) UNIQUE NOT NULL,
    user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,

    -- Anion information
    anion_name VARCHAR(50) NOT NULL,
    display_name VARCHAR(255),
    charge INTEGER DEFAULT -1,

    -- Input information
    identifier_type VARCHAR(20) NOT NULL,  -- "smiles" or "cas"
    identifier_value VARCHAR(500) NOT NULL,

    -- Job status
    status VARCHAR(20) NOT NULL DEFAULT 'pending',  -- pending, running, success, failed
    message TEXT,

    -- Generated files
    lt_path VARCHAR(500),
    pdb_path VARCHAR(500),

    -- Intermediate files
    work_dir VARCHAR(500),
    gaussian_log VARCHAR(500),
    mol2_file VARCHAR(500),
    gromacs_top VARCHAR(500),
    sob_output VARCHAR(500),

    -- Related QC job
    qc_job_id INTEGER REFERENCES qc_jobs(id) ON DELETE SET NULL,

    -- Metadata
    config JSONB DEFAULT '{}',

    -- Timestamps
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    started_at TIMESTAMP WITH TIME ZONE,
    finished_at TIMESTAMP WITH TIME ZONE,

    -- Soft delete
    is_deleted BOOLEAN DEFAULT FALSE,
    deleted_at TIMESTAMP WITH TIME ZONE,
    deleted_by INTEGER REFERENCES users(id) ON DELETE SET NULL,
    delete_reason VARCHAR(500)
);

-- Create indexes for anion_generation_jobs
CREATE INDEX IF NOT EXISTS idx_anion_gen_user_id ON anion_generation_jobs(user_id);
CREATE INDEX IF NOT EXISTS idx_anion_gen_status ON anion_generation_jobs(status);
CREATE INDEX IF NOT EXISTS idx_anion_gen_anion_name ON anion_generation_jobs(anion_name);
CREATE INDEX IF NOT EXISTS idx_anion_gen_created_at ON anion_generation_jobs(created_at);
CREATE INDEX IF NOT EXISTS idx_anion_gen_job_id ON anion_generation_jobs(job_id);

-- Create anion_library table
CREATE TABLE IF NOT EXISTS anion_library (
    id SERIAL PRIMARY KEY,
    anion_name VARCHAR(50) UNIQUE NOT NULL,
    display_name VARCHAR(255) NOT NULL,
    charge INTEGER DEFAULT -1,
    
    -- File paths
    lt_path VARCHAR(500) NOT NULL,
    pdb_path VARCHAR(500) NOT NULL,
    
    -- Source information
    source VARCHAR(100) DEFAULT 'manual',  -- "manual" or "auto_generated_sob_gaussian"
    generation_job_id INTEGER REFERENCES anion_generation_jobs(id) ON DELETE SET NULL,
    
    -- Metadata
    description TEXT,
    config JSONB DEFAULT '{}',
    
    -- Timestamps
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    created_by INTEGER REFERENCES users(id) ON DELETE SET NULL,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    
    -- Soft delete
    is_deleted BOOLEAN DEFAULT FALSE,
    deleted_at TIMESTAMP WITH TIME ZONE
);

-- Create indexes for anion_library
CREATE INDEX IF NOT EXISTS idx_anion_lib_name ON anion_library(anion_name);
CREATE INDEX IF NOT EXISTS idx_anion_lib_created_at ON anion_library(created_at);
CREATE INDEX IF NOT EXISTS idx_anion_lib_is_deleted ON anion_library(is_deleted);

