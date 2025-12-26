/**
 * 任务和配方信息展示组件
 */
import { Descriptions, Tag, Space, Typography, Tooltip, Card, theme } from 'antd';
import { InfoCircleOutlined, ThunderboltOutlined, ExperimentOutlined } from '@ant-design/icons';
import type { MDJob, ElectrolyteSystem } from '../../../types';
import type { PartitionInfo } from '../../../api/slurm';

interface JobInfoCardProps {
    job: MDJob;
    electrolyte: ElectrolyteSystem;
}

export function JobInfoCard({ job, electrolyte }: JobInfoCardProps) {
    const { token } = theme.useToken();

    return (
        <>
            {/* 配方信息 */}
            <Card
                title="电解质配方信息"
                style={{
                    marginBottom: 24,
                    borderRadius: 12,
                    boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
                    border: 'none'
                }}
            >
                <Descriptions bordered column={2}>
                    <Descriptions.Item label="配方名称">{electrolyte.name}</Descriptions.Item>
                    <Descriptions.Item label="配方 ID">#{electrolyte.id}</Descriptions.Item>
                    <Descriptions.Item label="温度">{electrolyte.temperature} K</Descriptions.Item>
                    <Descriptions.Item label="压力">{electrolyte.pressure} atm</Descriptions.Item>
                    <Descriptions.Item label="盒子大小">
                        {electrolyte.box_size ? Number(electrolyte.box_size).toFixed(1) : '-'} Å
                    </Descriptions.Item>
                    <Descriptions.Item label="力场">{electrolyte.force_field || 'OPLS-AA'}</Descriptions.Item>
                </Descriptions>
            </Card>

            {/* 组分信息 */}
            <Card
                title="配方组分"
                style={{
                    marginBottom: 24,
                    borderRadius: 12,
                    boxShadow: '0 2px 8px rgba(0,0,0,0.06)',
                    border: 'none'
                }}
            >
                <Space direction="vertical" style={{ width: '100%' }} size="middle">
                    {/* 溶剂 */}
                    {electrolyte.solvents && electrolyte.solvents.length > 0 && (
                        <div>
                            <div style={{ marginBottom: 8, fontWeight: 500 }}>溶剂：</div>
                            <Space wrap>
                                {electrolyte.solvents.map((s: any, idx: number) => (
                                    <Tag key={idx} color="blue">
                                        {s.name} ({s.mole_fraction || 0}%)
                                    </Tag>
                                ))}
                            </Space>
                        </div>
                    )}

                    {/* 阳离子 */}
                    {electrolyte.cations && electrolyte.cations.length > 0 && (
                        <div>
                            <div style={{ marginBottom: 8, fontWeight: 500 }}>阳离子：</div>
                            <Space wrap>
                                {electrolyte.cations.map((c: any, idx: number) => (
                                    <Tag key={idx} color="green">
                                        {c.name}
                                    </Tag>
                                ))}
                            </Space>
                        </div>
                    )}

                    {/* 阴离子 */}
                    {electrolyte.anions && electrolyte.anions.length > 0 && (
                        <div>
                            <div style={{ marginBottom: 8, fontWeight: 500 }}>阴离子：</div>
                            <Space wrap>
                                {electrolyte.anions.map((a: any, idx: number) => (
                                    <Tag key={idx} color="orange">
                                        {a.name}
                                    </Tag>
                                ))}
                            </Space>
                        </div>
                    )}
                </Space>
            </Card>
        </>
    );
}
