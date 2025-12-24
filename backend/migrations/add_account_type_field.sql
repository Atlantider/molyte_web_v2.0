-- Migration: Add account_type field to users table
-- Purpose: Support Scheme B refactoring - account type classification
-- Date: 2025-12-17

-- Step 1: Create ENUM type for account_type
DO $$ BEGIN
    CREATE TYPE account_type_enum AS ENUM (
        'personal',
        'master_account',
        'sub_account'
    );
EXCEPTION
    WHEN duplicate_object THEN NULL;
END $$;

-- Step 2: Add account_type column to users table
ALTER TABLE users
ADD COLUMN IF NOT EXISTS account_type account_type_enum DEFAULT 'personal' NOT NULL;

-- Step 3: Create index on account_type for faster queries
CREATE INDEX IF NOT EXISTS idx_users_account_type ON users(account_type);

-- Step 4: Verify migration
-- SELECT COUNT(*) as total_users,
--        account_type,
--        COUNT(*) as count
-- FROM users
-- GROUP BY account_type;

