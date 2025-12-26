/**
 * JobSubmit 主组件类型定义
 */

import type { MDJob, ElectrolyteSystem } from '../../types';
import type { PartitionInfo } from '../../api/slurm';
import type { DuplicateCheckResponse } from '../../api/qc';

export interface QCMoleculeParams {
    functional: string;
    basis_set: string;
    solvent_model: string;
}

export interface GlobalQCParams {
    accuracy_level: string;
    solvent_model: string;
    solvent_name: string;
    use_recommended_params: boolean;
    custom_eps?: number;
    custom_eps_inf?: number;
    custom_solvent_name?: string;
}

export interface JobSubmitState {
    loading: boolean;
    submitting: boolean;
    editMode: boolean;
    isSubmittedJob: boolean;

    // Data
    job: MDJob | null;
    electrolyte: ElectrolyteSystem | null;
    partitions: PartitionInfo[];

    // QC Parameters
    moleculeParams: Record<string, QCMoleculeParams>;
    editingMolecule: string | null;
    editingGlobalQC: boolean;
    globalQCParams: GlobalQCParams | null;

    // Duplicate Check
    duplicateCheckResult: DuplicateCheckResponse | null;
    checkingDuplicates: boolean;
}

export type JobSubmitAction =
    | { type: 'SET_LOADING'; loading: boolean }
    | { type: 'SET_SUBMITTING'; submitting: boolean }
    | { type: 'SET_JOB'; job: MDJob }
    | { type: 'SET_ELECTROLYTE'; electrolyte: ElectrolyteSystem }
    | { type: 'SET_PARTITIONS'; partitions: PartitionInfo[] }
    | { type: 'SET_MOLECULE_PARAMS'; key: string; params: QCMoleculeParams }
    | { type: 'SET_DUPLICATE_CHECK_RESULT'; result: DuplicateCheckResponse | null }
    | { type: 'SET_CHECKING_DUPLICATES'; checking: boolean };
