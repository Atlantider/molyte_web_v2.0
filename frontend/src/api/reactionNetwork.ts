/**
 * Reaction Network API Client
 * 反应网络API客户端
 */

import apiClient from './client';

// ============================================================================
// Types
// ============================================================================

export interface ReactionNetworkJobCreate {
    job_name: string;
    description?: string;
    initial_smiles: string[];
    temperature?: number;
    electrode_type?: 'anode' | 'cathode';
    voltage?: number;
    max_generations?: number;
    max_species?: number;
    energy_cutoff?: number;
    slurm_partition?: string;
    slurm_cpus?: number;
    slurm_time?: number;
}

export interface ReactionNetworkJob {
    id: number;
    user_id: number;
    job_name: string;
    description?: string;
    status: 'CREATED' | 'QUEUED' | 'RUNNING' | 'POSTPROCESSING' | 'COMPLETED' | 'FAILED' | 'CANCELLED';
    progress: number;
    initial_smiles: string[];
    temperature: number;
    electrode_type: string;
    voltage: number;
    max_generations: number;
    max_species: number;
    energy_cutoff: number;
    slurm_job_id?: string;
    work_dir?: string;
    num_molecules?: number;
    num_reactions?: number;
    max_generation_reached?: number;
    network_json_path?: string;
    visualization_png_path?: string;
    visualization_html_path?: string;
    created_at: string;
    updated_at?: string;
    started_at?: string;
    finished_at?: string;
    error_message?: string;
}

export interface Molecule {
    id: number;
    job_id: number;
    name: string;
    smiles: string;
    generation: number;
    energy_kcal?: number;
    molecular_weight?: number;
    num_atoms?: number;
    formal_charge?: number;
    properties: Record<string, any>;
    created_at: string;
}

export interface Reaction {
    id: number;
    job_id: number;
    reactant_smiles: string[];
    product_smiles: string[];
    operator_name?: string;
    reaction_type?: string;
    reaction_energy?: number;
    activation_energy?: number;
    driving_force?: string;
    metadata: Record<string, any>;
    created_at: string;
}

export interface NetworkNode {
    id: string;
    label: string;
    generation: number;
    energy?: number;
    properties: Record<string, any>;
}

export interface NetworkEdge {
    source: string;
    target: string;
    label?: string;
    operator?: string;
    energy?: number;
    properties: Record<string, any>;
}

export interface NetworkVisualizationData {
    nodes: NetworkNode[];
    edges: NetworkEdge[];
    statistics: Record<string, any>;
}

// ============================================================================
// API Functions
// ============================================================================

/**
 * 创建反应网络任务
 */
export async function createReactionNetworkJob(
    data: ReactionNetworkJobCreate
): Promise<ReactionNetworkJob> {
    const response = await apiClient.post('/reaction-network/jobs', data);
    return response.data;
}

/**
 * 获取任务列表
 */
export async function getReactionNetworkJobs(params?: {
    status?: string;
    skip?: number;
    limit?: number;
    sort_by?: string;
    sort_desc?: boolean;
}): Promise<{ total: number; jobs: ReactionNetworkJob[] }> {
    const response = await apiClient.get('/reaction-network/jobs', { params });
    return response.data;
}

/**
 * 获取任务详情
 */
export async function getReactionNetworkJob(jobId: number): Promise<ReactionNetworkJob & {
    molecules: Molecule[];
    reactions: Reaction[];
}> {
    const response = await apiClient.get(`/reaction-network/jobs/${jobId}`);
    return response.data;
}

/**
 * 提交任务到集群
 */
export async function submitReactionNetworkJob(jobId: number): Promise<{
    success: boolean;
    message: string;
    job_id: number;
    slurm_job_id?: string;
}> {
    const response = await apiClient.post(`/reaction-network/jobs/${jobId}/submit`);
    return response.data;
}

/**
 * 删除任务
 */
export async function deleteReactionNetworkJob(jobId: number): Promise<{
    message: string;
    job_id: number;
}> {
    const response = await apiClient.delete(`/reaction-network/jobs/${jobId}`);
    return response.data;
}

/**
 * 获取分子列表
 */
export async function getJobMolecules(jobId: number, params?: {
    generation?: number;
    min_energy?: number;
    max_energy?: number;
    skip?: number;
    limit?: number;
}): Promise<{ total: number; molecules: Molecule[] }> {
    const response = await apiClient.get(`/reaction-network/jobs/${jobId}/molecules`, { params });
    return response.data;
}

/**
 * 获取反应列表
 */
export async function getJobReactions(jobId: number, params?: {
    operator?: string;
    reaction_type?: string;
    min_energy?: number;
    max_energy?: number;
    skip?: number;
    limit?: number;
}): Promise<{ total: number; reactions: Reaction[] }> {
    const response = await apiClient.get(`/reaction-network/jobs/${jobId}/reactions`, { params });
    return response.data;
}

/**
 * 获取网络可视化数据
 */
export async function getNetworkVisualizationData(jobId: number, params?: {
    max_generation?: number;
}): Promise<NetworkVisualizationData> {
    const response = await apiClient.get(`/reaction-network/jobs/${jobId}/network`, { params });
    return response.data;
}

/**
 * 激活的算符信息
 */
export interface ActivatedOperator {
    name: string;
    description: string;
    weight: number;
    conditions: string[];
    activation_reasons: string[];
    required_drives: string[];
    enhancing_drives: string[];
    molecular_checks: string[];
}

export interface OperatorsResponse {
    job_id: number;
    num_operators: number;
    operators: ActivatedOperator[];
    environment: {
        temperature: number;
        voltage: number;
        electrode_type: string;
        active_drives: string[];
    };
}

/**
 * 获取激活的算符信息
 */
export async function getActivatedOperators(jobId: number): Promise<OperatorsResponse> {
    const response = await apiClient.get(`/reaction-network/jobs/${jobId}/operators`);
    return response.data;
}

