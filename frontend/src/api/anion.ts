/**
 * 阴离子生成 API
 */
import client from './client';

export interface AnionGenerationJob {
  id: number;
  job_id: string;
  anion_name: string;
  display_name: string;
  charge: number;
  identifier_type: string;
  identifier_value: string;
  status: string;
  message: string;
  lt_path?: string;
  pdb_path?: string;
  created_at: string;
  started_at?: string;
  finished_at?: string;
}

export interface AnionGenerationRequest {
  anion_name: string;
  display_name: string;
  charge: number;
  identifier_type: 'smiles' | 'cas';
  identifier_value: string;
}

/**
 * 获取用户的阴离子生成任务列表
 */
export const getAnionGenerationJobs = async (): Promise<{ jobs: AnionGenerationJob[] }> => {
  const response = await client.get('/forcefield/anions/jobs');
  return response.data;
};

/**
 * 提交新的阴离子生成任务
 */
export const submitAnionGenerationJob = async (data: AnionGenerationRequest): Promise<any> => {
  const response = await client.post('/forcefield/anions/auto-generate', data);
  return response.data;
};

/**
 * 获取阴离子生成任务的状态
 */
export const getAnionGenerationStatus = async (jobId: string): Promise<any> => {
  const response = await client.get(`/forcefield/anions/auto-generate/${jobId}`);
  return response.data;
};

/**
 * 删除阴离子生成任务
 */
export const deleteAnionGenerationJob = async (jobId: number): Promise<any> => {
  const response = await client.delete(`/forcefield/anions/jobs/${jobId}`);
  return response.data;
};

/**
 * 获取阴离子库
 */
export const getAnionLibrary = async (): Promise<any> => {
  const response = await client.get('/forcefield/anions/library');
  return response.data;
};

