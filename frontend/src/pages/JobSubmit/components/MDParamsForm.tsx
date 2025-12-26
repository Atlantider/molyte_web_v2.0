/**
 * MD参数配置表单组件
 */
import React from 'react';
import { Form, InputNumber, Select, Row, Col, Typography, Card, Space, Alert, Tag, Button, Input } from 'antd';
import { InfoCircleOutlined, EditOutlined } from '@ant-design/icons';
import AccuracyLevelSelector from '../../../components/AccuracyLevelSelector';
import type { MDJob, MDJobCreate } from '../../../types';

interface MDParamsFormProps {
    job: MDJob;
    editMode: boolean;
    onEditClick?: () => void;
}

export function MDParamsForm({ job, editMode, onEditClick }: MDParamsFormProps) {
    if (!editMode) {
        // 只读模式 - 显示当前配置
        return (
            <Card
                title={
                    <Space>
                        <span>MD计算参数配置</span>
                        {job?.config?.qc_enabled && <Tag color="blue">MD</Tag>}
                    </Space>
                }
                extra={
                    onEditClick && (
                        <Button icon={<EditOutlined />} onClick={onEditClick} style={{ borderRadius: 8 }}>
                            修改参数
                        </Button>
                    )
                }
                style={{
                    marginBottom: 24,
                    borderRadius: 12,
                    boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
                    border: 'none'
                }}
            >
                <Space direction="vertical" style={{ width: '100%' }} size="middle">
                    <Row gutter={16}>
                        <Col span={8}>
                            <div><strong>NPT步数：</strong>{job.config?.nsteps_npt?.toLocaleString() || 0}</div>
                        </Col>
                        <Col span={8}>
                            <div><strong>NVT步数：</strong>{job.config?.nsteps_nvt?.toLocaleString() || 0}</div>
                        </Col>
                        <Col span={8}>
                            <div><strong>时间步长：</strong>{job.config?.timestep || 0} fs</div>
                        </Col>
                    </Row>
                    <Row gutter={16}>
                        <Col span={8}>
                            <div><strong>温度：</strong>{job.config?.temperature || 0} K</div>
                        </Col>
                        <Col span={8}>
                            <div><strong>压力：</strong>{job.config?.pressure || 0} atm</div>
                        </Col>
                        <Col span={8}>
                            <div><strong>轨迹频率：</strong>{job.config?.freq_trj_npt || 0}</div>
                        </Col>
                    </Row>
                </Space>
            </Card>
        );
    }

    // 编辑模式
    return (
        <Card
            title={
                <Space>
                    <span>MD计算参数配置</span>
                    {job?.config?.qc_enabled && <Tag color="blue">MD</Tag>}
                </Space>
            }
            style={{
                marginBottom: 24,
                borderRadius: 12,
                boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
                border: 'none'
            }}
        >
            <Form.Item
                label="备注信息（可选）"
                name="user_note"
                tooltip="可选的备注信息，用于记录任务目的或特殊说明"
            >
                <Input placeholder="可选备注（如：高温测试、对照组等）" allowClear />
            </Form.Item>

            <Row gutter={16}>
                <Col span={8}>
                    <Form.Item
                        label="NPT步数"
                        name="nsteps_npt"
                        tooltip="NPT系综模拟步数"
                        rules={[{ required: true, message: '请输入NPT步数' }]}
                    >
                        <InputNumber min={0} step={100000} style={{ width: '100%' }} />
                    </Form.Item>
                </Col>
                <Col span={8}>
                    <Form.Item
                        label="NVT步数"
                        name="nsteps_nvt"
                        tooltip="NVT系综模拟步数"
                        rules={[{ required: true, message: '请输入NVT步数' }]}
                    >
                        <InputNumber min={0} step={100000} style={{ width: '100%' }} />
                    </Form.Item>
                </Col>
                <Col span={8}>
                    <Form.Item
                        label="时间步长 (fs)"
                        name="timestep"
                        tooltip="每步的时间增量"
                        rules={[{ required: true, message: '请输入时间步长' }]}
                    >
                        <InputNumber min={0.1} max={10} step={0.1} style={{ width: '100%' }} />
                    </Form.Item>
                </Col>
            </Row>

            <Row gutter={16}>
                <Col span={8}>
                    <Form.Item
                        label="温度 (K)"
                        name="temperature"
                        tooltip="模拟温度"
                        rules={[{ required: true, message: '请输入温度' }]}
                    >
                        <InputNumber min={0} max={1000} style={{ width: '100%' }} />
                    </Form.Item>
                </Col>
                <Col span={8}>
                    <Form.Item
                        label="压力 (atm)"
                        name="pressure"
                        tooltip="模拟压力"
                        rules={[{ required: true, message: '请输入压力' }]}
                    >
                        <InputNumber min={0} max={1000} style={{ width: '100%' }} />
                    </Form.Item>
                </Col>
                <Col span={8}>
                    <Form.Item
                        label="热力学输出频率"
                        name="thermo_freq"
                        tooltip="热力学数据输出频率"
                        rules={[{ required: true, message: '请输入输出频率' }]}
                    >
                        <InputNumber min={1} style={{ width: '100%' }} />
                    </Form.Item>
                </Col>
            </Row>

            <Row gutter={16}>
                <Col span={12}>
                    <Form.Item
                        label="NPT轨迹输出频率"
                        name="freq_trj_npt"
                        tooltip="NPT阶段轨迹输出频率"
                        rules={[{ required: true, message: '请输入轨迹频率' }]}
                    >
                        <InputNumber min={1} style={{ width: '100%' }} />
                    </Form.Item>
                </Col>
                <Col span={12}>
                    <Form.Item
                        label="NVT轨迹输出频率"
                        name="freq_trj_nvt"
                        tooltip="NVT阶段轨迹输出频率"
                        rules={[{ required: true, message: '请输入轨迹频率' }]}
                    >
                        <InputNumber min={1} style={{ width: '100%' }} />
                    </Form.Item>
                </Col>
            </Row>
        </Card>
    );
}
