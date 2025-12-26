/**
 * JobSubmit 主组件 - 重构版本
 * 
 * 将原1,251行的巨大文件拆分为模块化组件和hooks
 * 改进了可维护性、可测试性和代码可读性
 */
import React, { useState } from 'react';
import { useNavigate, useParams } from 'react-router-dom';
import {
    Button,
    Space,
    Typography,
    Alert,
    Form,
    Modal,
    theme,
    message,
} from 'antd';
import { ArrowLeftOutlined, ThunderboltOutlined, CheckCircleOutlined } from '@ant-design/icons';

// Hooks
import { useJobSubmit } from './hooks/useJobSubmit';
import { useDuplicateCheck } from './hooks/useDuplicateCheck';

// Components
import { JobInfoCard } from './components/JobInfoCard';
import { MDParamsForm } from './components/MDParamsForm';
import { SlurmConfigPanel } from './components/SlurmConfigPanel';
import { DuplicateCheckAlert } from './components/DuplicateCheckAlert';

const { Title, Text } = Typography;

export default function JobSubmit() {
    const navigate = useNavigate();
    const { jobId } = useParams<{ jobId: string }>();
    const { token } = theme.useToken();
    const [form] = Form.useForm();

    // Local state
    const [editMode, setEditMode] = useState(false);
    const [moleculeParams, setMoleculeParams] = useState<Record<string, any>>({});

    // Business logic hooks
    const {
        loading,
        submitting,
        job,
        electrolyte,
        partitions,
        isSubmittedJob,
        submitToCluster,
        saveConfig,
    } = useJobSubmit({
        jobId: Number(jobId),
        onSubmitSuccess: () => {
            message.success('任务提交成功！');
        },
    });

    const {
        checking: checkingDuplicates,
        result: duplicateCheckResult,
        checkDuplicates,
    } = useDuplicateCheck({
        job,
        electrolyte,
        moleculeParams,
    });

    // 处理编辑按钮点击
    const handleEditClick = () => {
        if (isSubmittedJob) {
            Modal.confirm({
                title: '创建新任务',
                content: '该任务已提交，无法直接修改。是否要基于当前配置创建一个新任务？',
                okText: '创建新任务',
                cancelText: '取消',
                onOk: () => {
                    setEditMode(true);
                },
            });
        } else {
            setEditMode(true);
        }
    };

    // 处理保存
    const handleSave = async () => {
        try {
            const values = await form.validateFields();
            const success = await saveConfig(values);
            if (success && !isSubmittedJob) {
                setEditMode(false);
            }
        } catch (error) {
            console.error('保存失败:', error);
        }
    };

    // 处理提交
    const handleSubmit = async () => {
        // 检查QC重复计算
        let duplicateInfo = duplicateCheckResult;
        if (job?.config?.qc_enabled && !duplicateInfo) {
            duplicateInfo = await checkDuplicates();
        }

        // 构建确认内容
        let confirmContent: React.ReactNode = '确定要将此任务提交到 Slurm 集群执行吗？提交后将开始计算。';

        if (duplicateInfo && duplicateInfo.existing_count > 0) {
            confirmContent = (
                <div>
                    <p>确定要将此任务提交到 Slurm 集群执行吗？</p>
                    <Alert
                        type="success"
                        showIcon
                        icon={<CheckCircleOutlined />}
                        style={{ marginTop: 12 }}
                        message={
                            <span>
                                检测到 <strong>{duplicateInfo.existing_count}</strong> 个分子已有计算结果，
                                将直接复用，无需重复计算！
                            </span>
                        }
                        description={
                            duplicateInfo.new_count > 0 ? (
                                <span>另外 {duplicateInfo.new_count} 个分子将执行新计算。</span>
                            ) : (
                                <span>所有QC计算都将复用已有结果，节省计算时间和资源！</span>
                            )
                        }
                    />
                </div>
            );
        }

        Modal.confirm({
            title: '确认提交任务到集群',
            content: confirmContent,
            okText: '确定提交',
            cancelText: '取消',
            onOk: async () => {
                await submitToCluster();
            },
        });
    };

    // 初始化表单值
    React.useEffect(() => {
        if (job?.config) {
            const { job_name, ...configWithoutJobName } = job.config;
            form.setFieldsValue(configWithoutJobName);
        }
    }, [job, form]);

    // 加载状态
    if (loading || !job || !electrolyte) {
        return (
            <div style={{
                display: 'flex',
                justifyContent: 'center',
                alignItems: 'center',
                height: 'calc(100vh - 64px)',
                background: token.colorBgLayout,
            }}>
                加载中...
            </div>
        );
    }

    return (
        <div style={{ padding: 24, background: token.colorBgLayout, minHeight: 'calc(100vh - 64px)' }}>
            {/* 页面头部 */}
            <div style={{ marginBottom: 24 }}>
                <Space style={{ marginBottom: 16 }}>
                    <Button
                        icon={<ArrowLeftOutlined />}
                        onClick={() => navigate('/workspace/liquid-electrolyte/md')}
                        style={{ borderRadius: 8 }}
                    >
                        返回任务列表
                    </Button>
                </Space>
                <Title level={2} style={{ margin: 0, marginBottom: 8 }}>
                    <ThunderboltOutlined style={{ marginRight: 12, color: '#1677ff' }} />
                    {isSubmittedJob ? '查看任务配置' : '提交计算任务'}
                </Title>
                <Text type="secondary">
                    {isSubmittedJob ? '查看已提交任务的配置参数' : '检查并确认计算参数后提交到集群'}
                </Text>
            </div>

            {/* 状态提示 */}
            {isSubmittedJob ? (
                <Alert
                    message="任务已提交"
                    description="该任务已提交到集群，无法直接修改。如需修改参数，可以基于当前配置创建新任务。"
                    type="warning"
                    showIcon
                    style={{ marginBottom: 24, borderRadius: 8 }}
                />
            ) : (
                <Alert
                    message="检查计算参数"
                    description="请仔细检查以下计算参数，确认无误后提交到集群执行"
                    type="info"
                    showIcon
                    style={{ marginBottom: 24, borderRadius: 8 }}
                />
            )}

            {/* 重复计算提示 */}
            <DuplicateCheckAlert result={duplicateCheckResult} />

            {/* 任务和配方信息 */}
            <JobInfoCard job={job} electrolyte={electrolyte} />

            {/* 表单区域 */}
            <Form form={form} layout="vertical">
                {/* MD参数配置 */}
                <MDParamsForm
                    job={job}
                    editMode={editMode}
                    onEditClick={!editMode ? handleEditClick : undefined}
                />

                {/* Slurm配置 */}
                <SlurmConfigPanel
                    job={job}
                    partitions={partitions}
                    editMode={editMode}
                />
            </Form>

            {/* 底部操作按钮 */}
            <div style={{
                position: 'sticky',
                bottom: 0,
                padding: '16px 0',
                background: token.colorBgLayout,
                borderTop: `1px solid ${token.colorBorder}`,
                marginTop: 24
            }}>
                <Space>
                    {editMode ? (
                        <>
                            <Button onClick={() => setEditMode(false)}>取消</Button>
                            <Button type="primary" onClick={handleSave} loading={submitting}>
                                {isSubmittedJob ? '创建新任务' : '保存配置'}
                            </Button>
                        </>
                    ) : (
                        <>
                            {!isSubmittedJob && (
                                <Button
                                    type="primary"
                                    size="large"
                                    onClick={handleSubmit}
                                    loading={submitting || checkingDuplicates}
                                    style={{ borderRadius: 8 }}
                                >
                                    提交到集群
                                </Button>
                            )}
                        </>
                    )}
                </Space>
            </div>
        </div>
    );
}
