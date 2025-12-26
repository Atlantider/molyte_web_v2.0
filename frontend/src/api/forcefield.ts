/**
 * Force field and anion generation API client
 */
import client from './client';

export interface AnionGenerationRequest {
  anion_name: string;
  display_name?: string;
  charge?: number;
  identifier_type: 'smiles' | 'cas';
  identifier_value: string;
}

export interface AnionGenerationResponse {
  job_id: string;
  status: string;
}

export interface AnionGenerationStatusResponse {
  job_id: string;
  status: 'pending' | 'running' | 'success' | 'failed';
  message?: string;
  anion_key?: string;
  files?: {
    lt_path: string;
    pdb_path: string;
  };
  created_at: string;
  updated_at: string;
  started_at?: string;
  finished_at?: string;
}

export interface AnionLibraryEntry {
  anion_name: string;
  display_name: string;
  charge: number;
  lt_path: string;
  pdb_path: string;
  source: string;
  description?: string;
  created_at: string;
}

export interface AnionLibraryListResponse {
  anions: AnionLibraryEntry[];
  total: number;
}

/**
 * Submit a job to auto-generate anion force field
 */
export async function submitAnionGeneration(
  request: AnionGenerationRequest
): Promise<AnionGenerationResponse> {
  const response = await client.post<AnionGenerationResponse>(
    '/forcefield/anions/auto-generate',
    request
  );
  return response.data;
}

/**
 * Get the status of an anion generation job
 */
export async function getAnionGenerationStatus(
  jobId: string
): Promise<AnionGenerationStatusResponse> {
  const response = await client.get<AnionGenerationStatusResponse>(
    `/forcefield/anions/auto-generate/${jobId}`
  );
  return response.data;
}

/**
 * Get list of available anions in the library
 */
export async function getAnionLibrary(): Promise<AnionLibraryListResponse> {
  const response = await client.get<AnionLibraryListResponse>(
    '/forcefield/anions/library'
  );
  return response.data;
}

