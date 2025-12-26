/**
 * 任务提交主逻辑Hook
 * 
 * 管理任务加载、提交、保存等核心业务逻辑
 */
import { useState, useEffect, useCallback } from 'react';
import { useNavigate } from 'react-router-dom';
import { message, Modal } from 'antd';
import { WalletOutlined } from '@ant-design/icons';
import type { MDJob, ElectrolyteSystem, MDJobCreate } from '../../../types';
import type { PartitionInfo } from '../../../api/slurm';
import { getMDJob, updateMDJobConfig, submitJobToCluster, createMDJob } from '../../../api/jobs';
import { getElectrolyte } from '../../../api/electrolytes';
import { checkCanSubmit } from '../../../api/billing';
import { getPartitions } from '../../../api/slurm';

interface UseJobSubmitOptions {
    jobId: number;
    onSubmitSuccess?: () => void;
}

export function useJobSubmit({ jobId, onSubmitSuccess }: UseJobSubmitOptions) {
    const navigate = useNavigate();

    // State
    const [loading, setLoading] = useState(false);
    const [submitting, setSubmitting] = useState(false);
    const [job, setJob] = useState<MDJob | null>(null);
    const [electrolyte, setElectrolyte] = useState<ElectrolyteSystem | null>(null);
    const [partitions, setPartitions] = useState<PartitionInfo[]>([]);
    const [isSubmittedJob, setIsSubmittedJob] = useState(false);

    // 加载任务数据
    const loadJobData = useCallback(async () => {
        setLoading(true);
        try {
            const jobData = await getMDJob(jobId);
            setJob(jobData);

            // 检查任务状态
            setIsSubmittedJob(jobData.status !== 'CREATED' && jobData.status !== 'CANCELLED');

            // 加载配方数据
            const electrolyteData = await getElectrolyte(jobData.system_id);
            setElectrolyte(electrolyteData);

            // 加载Slurm分区
            try {
                const partitionsData = await getPartitions();
                setPartitions(partitionsData);
            } catch (err) {
                console.error('加载分区信息失败:', err);
                setPartitions([{
                    name: 'cpu',
                    state: 'up',
                    total_nodes: 0,
                    available_nodes: 0,
                    total_cpus: 0,
                    available_cpus: 0
                }]);
            }

            return { job: jobData, electrolyte: electrolyteData };
        } catch (error: any) {
            message.error('加载任务信息失败: ' + (error.response?.data?.detail || error.message));
            navigate('/workspace/liquid-electrolyte/md');
            throw error;
        } finally {
            setLoading(false);
        }
    }, [jobId, navigate]);

    // 提交到集群
    const submitToCluster = useCallback(async () => {
        // 检查余额
        try {
            const canSubmitResult = await checkCanSubmit();
            if (!canSubmitResult.can_submit) {
                Modal.confirm({
                    title: '余额不足',
                    icon: <WalletOutlined style={{ color: '#faad14' }} />,
                    content: (
                        <div>
                            <p>{canSubmitResult.reason} </p>
                            < p > 是否前往充值？</p>
                        </div>
                    ),
                    okText: '前往充值',
                    cancelText: '取消',
                    onOk: () => navigate('/workspace/recharge'),
                });
                return false;
            }
        } catch (error) {
            console.error('检查余额失败:', error);
        }

        // 提交
        try {
            setSubmitting(true);
            await submitJobToCluster(jobId);
            message.success('任务已提交到集群！');
            onSubmitSuccess?.();
            navigate(`/workspace/liquid-electrolyte/md/${jobId}`);
            return true;
        } catch (error: any) {
            const detail = error.response?.data?.detail || error.message;

            if (error.response?.status === 402) {
                Modal.confirm({
                    title: '余额不足',
                    icon: <WalletOutlined style={{ color: '#faad14' }} />,
                    content: (
                        <div>
                            <p>{detail} </p>
                            < p > 是否前往充值？</p>
                        </div>
                    ),
                    okText: '前往充值',
                    cancelText: '取消',
                    onOk: () => navigate('/workspace/recharge'),
                });
            } else {
                message.error('提交失败: ' + detail);
            }
            return false;
        } finally {
            setSubmitting(false);
        }
    }, [jobId, navigate, onSubmitSuccess]);

    // 保存配置
    const saveConfig = useCallback(async (formValues: any) => {
        try {
            setSubmitting(true);

            // 构建QC配置
            const qcOptions = job?.config?.qc_enabled ? {
                enabled: true,
                accuracy_level: formValues.qc_accuracy_level || job.config?.qc_accuracy_level || 'standard',
                basis_set: formValues.qc_basis_set || job.config?.qc_basis_set || '6-31++g(d,p)',
                functional: formValues.qc_functional || job.config?.qc_functional || 'B3LYP',
                solvent_model: formValues.qc_solvent_model || job.config?.qc_solvent_model || 'pcm',
                solvent_name: formValues.qc_solvent_name || job.config?.qc_solvent_name || 'water',
                use_recommended_params: formValues.qc_use_recommended_params !== undefined
                    ? formValues.qc_use_recommended_params
                    : (job.config?.qc_use_recommended_params !== false),
            } : undefined;

            if (isSubmittedJob) {
                // 已提交任务：创建新任务
                const generateCopyName = (originalName: string) => {
                    const copyMatch = originalName.match(/-copy(-(\d+))?$/);
                    if (copyMatch) {
                        const copyNumber = copyMatch[2] ? parseInt(copyMatch[2]) + 1 : 2;
                        return originalName.replace(/-copy(-\d+)?$/, `-copy-${copyNumber}`);
                    }
                    return `${originalName}-copy`;
                };

                const originalName = job!.config?.job_name || '';
                const newJobName = generateCopyName(originalName);

                const newJobData: MDJobCreate = {
                    system_id: job!.system_id,
                    job_name: formValues.user_note || undefined,
                    nsteps_npt: formValues.nsteps_npt,
                    nsteps_nvt: formValues.nsteps_nvt,
                    timestep: formValues.timestep,
                    temperature: formValues.temperature,
                    pressure: formValues.pressure,
                    freq_trj_npt: formValues.freq_trj_npt,
                    freq_trj_nvt: formValues.freq_trj_nvt,
                    thermo_freq: formValues.thermo_freq,
                    slurm_partition: formValues.slurm_partition,
                    slurm_nodes: formValues.slurm_nodes,
                    slurm_ntasks: formValues.slurm_ntasks,
                    slurm_cpus_per_task: formValues.slurm_cpus_per_task,
                    slurm_time: formValues.slurm_time,
                    submit_to_cluster: false,
                    qc_options: qcOptions,
                };

                const newJob = await createMDJob(newJobData);
                message.success(`已创建新任务：${newJobName}`);
                navigate(`/workspace/liquid-electrolyte/md/${newJob.id}/submit`);
            } else {
                // 未提交任务：直接更新
                const updateData = {
                    ...formValues,
                    slurm_partition: formValues.slurm_partition,
                    slurm_nodes: formValues.slurm_nodes,
                    slurm_ntasks: formValues.slurm_ntasks,
                    slurm_cpus_per_task: formValues.slurm_cpus_per_task,
                    slurm_time: formValues.slurm_time,
                    ...(job?.config?.qc_enabled && {
                        qc_accuracy_level: formValues.qc_accuracy_level || job.config?.qc_accuracy_level,
                        qc_basis_set: formValues.qc_basis_set || job.config?.qc_basis_set,
                        qc_functional: formValues.qc_functional || job.config?.qc_functional,
                        qc_solvent_model: formValues.qc_solvent_model || job.config?.qc_solvent_model,
                        qc_solvent_name: formValues.qc_solvent_name || job.config?.qc_solvent_name,
                        qc_use_recommended_params: formValues.qc_use_recommended_params !== undefined
                            ? formValues.qc_use_recommended_params
                            : job.config?.qc_use_recommended_params,
                    }),
                };

                const updatedJob = await updateMDJobConfig(jobId, updateData);
                setJob(updatedJob);
                message.success('配置已更新！');
            }

            return true;
        } catch (error: any) {
            if (error.errorFields) {
                message.error('请检查表单填写');
            } else {
                message.error('保存失败: ' + (error.response?.data?.detail || error.message));
            }
            return false;
        } finally {
            setSubmitting(false);
        }
    }, [job, jobId, isSubmittedJob, navigate]);

    // 初始化加载
    useEffect(() => {
        if (jobId) {
            loadJobData();
        }
    }, [jobId, loadJobData]);

    return {
        // State
        loading,
        submitting,
        job,
        electrolyte,
        partitions,
        isSubmittedJob,

        // Actions
        loadJobData,
        submitToCluster,
        saveConfig,
    };
}
