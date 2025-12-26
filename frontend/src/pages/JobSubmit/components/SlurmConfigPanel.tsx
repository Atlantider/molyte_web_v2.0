/**
 * Slurm配置面板组件
 */
import React from 'react';
import { Card, Form, Select, InputNumber, Row, Col, Space, Tag, Button, Tooltip } from 'antd';
import { BulbOutlined, SyncOutlined } from '@ant-design/icons';
import type { PartitionInfo } from '../../api/slurm';
import type { MDJob } from '../../types';

interface SlurmConfigPanelProps {
    job: MDJob;
    partitions: PartitionInfo[];
    editMode: boolean;
    onGetSuggestion?: () => void;
}

export function SlurmConfigPanel({ job, partitions, editMode, onGetSuggestion }: SlurmConfigPanelProps) {
    if (!editMode) {
        // 只读模式 - 显示当前配置
        return (
            <Card
                title="Slurm 集群资源配置"
                style={{
                    marginBottom: 24,
                    borderRadius: 12,
                    boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
                    border: 'none'
                }}
            >
                <Space direction="vertical" style={{ width: '100%' }}>
                    <div><strong>分区：</strong>{job.config?.slurm_partition || 'cpu'}</div>
                    <div><strong>节点数：</strong>{job.config?.slurm_nodes || 1}</div>
                    <div><strong>任务数：</strong>{job.config?.slurm_ntasks || 32}</div>
                    <div><strong>每任务CPU数：</strong>{job.config?.slurm_cpus_per_task || 1}</div>
                    <div><strong>时间限制：</strong>{job.config?.slurm_time || '24:00:00'}</div>
                </Space>
            </Card>
        );
    }

    // 编辑模式
    return (
        <Card
            title={
                <Space>
                    <span>Slurm 集群资源配置</span>
                    {onGetSuggestion && (
                        <Tooltip title="获取推荐配置">
                            <Button
                                type="link"
                                size="small"
                                icon={<BulbOutlined />}
                                onClick={onGetSuggestion}
                            >
                                智能推荐
                            </Button>
                        </Tooltip>
                    )}
                </Space>
            }
            style={{
                marginBottom: 24,
                borderRadius: 12,
                boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
                border: 'none'
            }}
        >
            <Row gutter={16}>
                <Col span={12}>
                    <Form.Item
                        label="分区 (Partition)"
                        name="slurm_partition"
                        tooltip="选择运行任务的计算分区"
                        rules={[{ required: true, message: '请选择分区' }]}
                    >
                        <Select placeholder="选择分区">
                            {partitions.map(p => (
                                <Select.Option key={p.name} value={p.name}>
                                    <Space>
                                        {p.name}
                                        <Tag color={p.state === 'up' ? 'success' : 'error'}>
                                            {p.available_cpus}/{p.total_cpus} CPUs
                                        </Tag>
                                    </Space>
                                </Select.Option>
                            ))}
                        </Select>
                    </Form.Item>
                </Col>

                <Col span={12}>
                    <Form.Item
                        label="节点数 (Nodes)"
                        name="slurm_nodes"
                        tooltip="使用的计算节点数量"
                        rules={[{ required: true, message: '请输入节点数' }]}
                    >
                        <InputNumber min={1} max={10} style={{ width: '100%' }} />
                    </Form.Item>
                </Col>
            </Row>

            <Row gutter={16}>
                <Col span={12}>
                    <Form.Item
                        label="任务数 (Tasks)"
                        name="slurm_ntasks"
                        tooltip="并行任务数，通常设置为总CPU核数"
                        rules={[{ required: true, message: '请输入任务数' }]}
                    >
                        <InputNumber min={1} max={128} style={{ width: '100%' }} />
                    </Form.Item>
                </Col>

                <Col span={12}>
                    <Form.Item
                        label="每任务CPU数"
                        name="slurm_cpus_per_task"
                        tooltip="每个任务使用的CPU核数"
                        rules={[{ required: true, message: '请输入CPU数' }]}
                    >
                        <InputNumber min={1} max={32} style={{ width: '100%' }} />
                    </Form.Item>
                </Col>
            </Row>

            <Form.Item
                label="时间限制 (Time Limit)"
                name="slurm_time"
                tooltip="任务运行的最大时间，格式：HH:MM:SS 或 D-HH:MM:SS"
                rules={[
                    { required: true, message: '请输入时间限制' },
                    {
                        pattern: /^(\d+-)?(\d{1,2}):(\d{2}):(\d{2})$/,
                        message: '格式错误，请使用 HH:MM:SS 或 D-HH:MM:SS'
                    }
                ]}
            >
                <Select placeholder="选择时间限制">
                    <Select.Option value="06:00:00">6小时</Select.Option>
                    <Select.Option value="12:00:00">12小时</Select.Option>
                    <Select.Option value="24:00:00">24小时</Select.Option>
                    <Select.Option value="48:00:00">48小时</Select.Option>
                    <Select.Option value="3-00:00:00">3天</Select.Option>
                    <Select.Option value="7-00:00:00">7天</Select.Option>
                </Select>
            </Form.Item>
        </Card>
    );
}
