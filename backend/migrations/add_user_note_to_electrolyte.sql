-- Add user_note field to electrolyte_systems table
-- This field stores the user's custom name/note, separate from the system-generated name

ALTER TABLE electrolyte_systems
ADD COLUMN user_note VARCHAR(255) NULL;

-- Create index for user_note if needed
CREATE INDEX idx_electrolyte_systems_user_note ON electrolyte_systems(user_note);

