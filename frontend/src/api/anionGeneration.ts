/**
 * 阴离子力场自动生成 API
 */
import client from './client';

export interface AnionGenerationRequest {
  anion_name: string;
  display_name: string;
  charge: number;
  identifier_type: 'smiles' | 'cas';
  identifier_value: string;
}

export interface AnionGenerationResponse {
  job_id: string;
  status: 'pending' | 'running' | 'success' | 'failed';
  message?: string;
}

export interface AnionGenerationStatusResponse {
  job_id: string;
  status: 'pending' | 'running' | 'success' | 'failed';
  message?: string;
  current_step?: string;
  anion_key?: string;
  files?: {
    lt_path: string;
    pdb_path: string;
  };
  created_at: string;
  updated_at: string;
  completed_at?: string;
}

/**
 * 提交阴离子力场自动生成任务
 */
export const submitAnionGeneration = async (
  request: AnionGenerationRequest
): Promise<AnionGenerationResponse> => {
  const response = await client.post('/api/v1/forcefield/anions/auto-generate', request);
  return response.data;
};

/**
 * 查询阴离子生成任务状态
 */
export const getAnionGenerationStatus = async (
  jobId: string
): Promise<AnionGenerationStatusResponse> => {
  const response = await client.get(`/api/v1/forcefield/anions/auto-generate/${jobId}`);
  return response.data;
};

