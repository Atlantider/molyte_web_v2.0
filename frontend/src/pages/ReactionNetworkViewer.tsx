/**
 * Reaction Network Visualization Component
 * 反应网络交互式可视化 - 使用react-force-graph
 */

import React, { useState, useEffect, useRef } from 'react';
import {
    Card,
    Button,
    Space,
    Select,
    Slider,
    Switch,
    Row,
    Col,
    Statistic,
    Tag,
    Drawer,
    Descriptions,
    message,
    Tabs
} from 'antd';
import {
    FullscreenOutlined,
    DownloadOutlined,
    ReloadOutlined,
    ZoomInOutlined,
    ZoomOutOutlined,
    InfoCircleOutlined
} from '@ant-design/icons';
import ForceGraph2D from 'react-force-graph-2d';
import { useParams } from 'react-router-dom';
import { getNetworkVisualizationData, type NetworkNode, type NetworkEdge } from '../api/reactionNetwork';

const { Option } = Select;
const { TabPane } = Tabs;

interface NetworkData {
    nodes: NetworkNode[];
    links: NetworkEdge[];
    statistics: Record<string, any>;
}

const ReactionNetworkViewer: React.FC = () => {
    const { id } = useParams<{ id: string }>();
    const graphRef = useRef<any>();

    const [networkData, setNetworkData] = useState<NetworkData | null>(null);
    const [loading, setLoading] = useState(true);
    const [selectedNode, setSelectedNode] = useState<NetworkNode | null>(null);
    const [selectedLink, setSelectedLink] = useState<NetworkEdge | null>(null);
    const [drawerVisible, setDrawerVisible] = useState(false);

    // 可视化控制
    const [filterGeneration, setFilterGeneration] = useState<number | null>(null);
    const [colorScheme, setColorScheme] = useState<'generation' | 'energy'>('generation');
    const [showLabels, setShowLabels] = useState(true);
    const [linkDistance, setLinkDistance] = useState(150);
    const [chargeStrength, setChargeStrength] = useState(-300);

    useEffect(() => {
        loadNetworkData();
    }, [id, filterGeneration]);

    const loadNetworkData = async () => {
        setLoading(true);
        try {
            const data = await getNetworkVisualizationData(Number(id), {
                max_generation: filterGeneration || undefined
            });

            // 转换数据格式为react-force-graph需要的格式
            const nodes = data.nodes.map(node => ({
                id: node.id,
                label: node.label,
                generation: node.generation,
                energy: node.energy,
                properties: node.properties,
                // 颜色配置
                color: getNodeColor(node, colorScheme),
                size: 8 + (data.nodes.length > 50 ? 0 : 4) // 大网络时节点小一点
            }));

            const links = data.edges.map(edge => ({
                source: edge.source,
                target: edge.target,
                label: edge.label,
                operator: edge.operator,
                energy: edge.energy,
                properties: edge.properties,
                color: getLinkColor(edge)
            }));

            setNetworkData({ nodes, links, statistics: data.statistics });
        } catch (error) {
            message.error('加载网络数据失败');
        } finally {
            setLoading(false);
        }
    };

    const getNodeColor = (node: NetworkNode, scheme: 'generation' | 'energy') => {
        if (scheme === 'generation') {
            // 按代数着色
            const colors = ['#1890ff', '#52c41a', '#faad14', '#f5222d', '#722ed1', '#eb2f96'];
            return colors[node.generation % colors.length];
        } else {
            // 按能量着色
            if (!node.energy) return '#d9d9d9';
            const energy = node.energy;
            if (energy < -50) return '#0050b3';
            if (energy < 0) return '#1890ff';
            if (energy < 50) return '#52c41a';
            return '#f5222d';
        }
    };

    const getLinkColor = (edge: NetworkEdge) => {
        if (!edge.energy) return '#d9d9d9';
        const energy = edge.energy;
        if (energy < 0) return '#52c41a'; // 放能反应
        if (energy < 50) return '#faad14';
        return '#f5222d'; // 高能反应
    };

    const handleNodeClick = (node: any) => {
        setSelectedNode(node);
        setSelectedLink(null);
        setDrawerVisible(true);
    };

    const handleLinkClick = (link: any) => {
        setSelectedLink(link);
        setSelectedNode(null);
        setDrawerVisible(true);
    };

    const handleDownload = () => {
        if (!graphRef.current) return;

        // 截图下载
        const canvas = graphRef.current.renderer().domElement;
        canvas.toBlob((blob: Blob) => {
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = `reaction_network_${id}.png`;
            a.click();
            URL.revokeObjectURL(url);
        });

        message.success('网络图已保存');
    };

    const handleExportJSON = () => {
        if (!networkData) return;

        const dataStr = JSON.stringify(networkData, null, 2);
        const blob = new Blob([dataStr], { type: 'application/json' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `network_data_${id}.json`;
        a.click();
        URL.revokeObjectURL(url);

        message.success('数据已导出');
    };

    const handleZoomToFit = () => {
        if (graphRef.current) {
            graphRef.current.zoomToFit(400);
        }
    };

    if (loading || !networkData) {
        return <div style={{ padding: '100px', textAlign: 'center' }}>加载中...</div>;
    }

    return (
        <div style={{ height: '100vh', display: 'flex', flexDirection: 'column' }}>
            {/* 控制面板 */}
            <Card bodyStyle={{ padding: '12px 24px' }}>
                <Row gutter={16} align="middle">
                    <Col flex="auto">
                        <Space size="large">
                            <div>
                                <span style={{ marginRight: 8 }}>代数筛选:</span>
                                <Select
                                    value={filterGeneration}
                                    onChange={setFilterGeneration}
                                    style={{ width: 120 }}
                                    allowClear
                                    placeholder="全部"
                                >
                                    {[0, 1, 2, 3, 4, 5].map(gen => (
                                        <Option key={gen} value={gen}>第 {gen} 代</Option>
                                    ))}
                                </Select>
                            </div>

                            <div>
                                <span style={{ marginRight: 8 }}>配色方案:</span>
                                <Select
                                    value={colorScheme}
                                    onChange={setColorScheme}
                                    style={{ width: 120 }}
                                >
                                    <Option value="generation">按代数</Option>
                                    <Option value="energy">按能量</Option>
                                </Select>
                            </div>

                            <div>
                                <span style={{ marginRight: 8 }}>显示标签:</span>
                                <Switch checked={showLabels} onChange={setShowLabels} />
                            </div>
                        </Space>
                    </Col>

                    <Col>
                        <Space>
                            <Button icon={<ZoomInOutlined />} onClick={() => graphRef.current?.zoom(1.5, 400)}>
                                放大
                            </Button>
                            <Button icon={<ZoomOutOutlined />} onClick={() => graphRef.current?.zoom(0.75, 400)}>
                                缩小
                            </Button>
                            <Button icon={<FullscreenOutlined />} onClick={handleZoomToFit}>
                                适应窗口
                            </Button>
                            <Button icon={<DownloadOutlined />} onClick={handleDownload}>
                                下载图片
                            </Button>
                            <Button icon={<DownloadOutlined />} onClick={handleExportJSON}>
                                导出JSON
                            </Button>
                            <Button icon={<ReloadOutlined />} onClick={loadNetworkData}>
                                刷新
                            </Button>
                        </Space>
                    </Col>
                </Row>

                {/* 统计信息 */}
                <Row gutter={16} style={{ marginTop: '16px' }}>
                    <Col span={6}>
                        <Statistic
                            title="节点数"
                            value={networkData.statistics.num_nodes}
                            prefix={<InfoCircleOutlined />}
                        />
                    </Col>
                    <Col span={6}>
                        <Statistic
                            title="边数"
                            value={networkData.statistics.num_edges}
                            valueStyle={{ color: '#52c41a' }}
                        />
                    </Col>
                    <Col span={6}>
                        <Statistic
                            title="最大代数"
                            value={networkData.statistics.max_generation}
                            valueStyle={{ color: '#faad14' }}
                        />
                    </Col>
                    <Col span={6}>
                        <div>
                            <div style={{ fontSize: '14px', color: '#666', marginBottom: '4px' }}>代数分布</div>
                            <Space wrap>
                                {Object.entries(networkData.statistics.generation_distribution || {}).map(([gen, count]) => (
                                    <Tag key={gen} color="blue">Gen {gen}: {count as number}</Tag>
                                ))}
                            </Space>
                        </div>
                    </Col>
                </Row>

                {/* 高级控制 */}
                <Row gutter={16} style={{ marginTop: '16px' }}>
                    <Col span={12}>
                        <div>
                            <span style={{ marginRight: 8 }}>连接距离: {linkDistance}</span>
                            <Slider
                                min={50}
                                max={300}
                                value={linkDistance}
                                onChange={setLinkDistance}
                                style={{ width: 200 }}
                            />
                        </div>
                    </Col>
                    <Col span={12}>
                        <div>
                            <span style={{ marginRight: 8 }}>排斥力: {chargeStrength}</span>
                            <Slider
                                min={-500}
                                max={-100}
                                value={chargeStrength}
                                onChange={setChargeStrength}
                                style={{ width: 200 }}
                            />
                        </div>
                    </Col>
                </Row>
            </Card>

            {/* 力导向图 */}
            <div style={{ flex: 1, background: '#f0f2f5' }}>
                <ForceGraph2D
                    ref={graphRef}
                    graphData={{ nodes: networkData.nodes, links: networkData.links }}
                    nodeLabel={(node: any) => showLabels ? node.label : ''}
                    nodeColor={(node: any) => node.color}
                    nodeRelSize={6}
                    nodeVal={(node: any) => node.size}
                    linkLabel={(link: any) => link.label || link.operator}
                    linkColor={(link: any) => link.color}
                    linkDirectionalArrowLength={6}
                    linkDirectionalArrowRelPos={1}
                    linkWidth={2}
                    linkCurvature={0.2}
                    onNodeClick={handleNodeClick}
                    onLinkClick={handleLinkClick}
                    d3AlphaDecay={0.02}
                    d3VelocityDecay={0.3}
                    warmupTicks={100}
                    cooldownTicks={0}
                    linkDistance={linkDistance}
                    d3Force={{
                        charge: { strength: chargeStrength },
                        center: { x: 0, y: 0 }
                    }}
                />
            </div>

            {/* 详情抽屉 */}
            <Drawer
                title={selectedNode ? '分子详情' : '反应详情'}
                placement="right"
                width={450}
                onClose={() => setDrawerVisible(false)}
                visible={drawerVisible}
            >
                {selectedNode && (
                    <Tabs defaultActiveKey="basic">
                        <TabPane tab="基本信息" key="basic">
                            <Descriptions column={1} bordered>
                                <Descriptions.Item label="名称">{selectedNode.label}</Descriptions.Item>
                                <Descriptions.Item label="SMILES">
                                    <code style={{ fontSize: '12px', wordBreak: 'break-all' }}>{selectedNode.id}</code>
                                </Descriptions.Item>
                                <Descriptions.Item label="代数">
                                    <Tag color="blue">第 {selectedNode.generation} 代</Tag>
                                </Descriptions.Item>
                                <Descriptions.Item label="能量">
                                    {selectedNode.energy ? `${selectedNode.energy.toFixed(2)} kcal/mol` : '-'}
                                </Descriptions.Item>
                            </Descriptions>
                        </TabPane>
                        <TabPane tab="分子属性" key="properties">
                            <Descriptions column={1} bordered>
                                {Object.entries(selectedNode.properties || {}).map(([key, value]) => (
                                    <Descriptions.Item key={key} label={key}>
                                        {String(value)}
                                    </Descriptions.Item>
                                ))}
                            </Descriptions>
                        </TabPane>
                    </Tabs>
                )}

                {selectedLink && (
                    <Descriptions column={1} bordered>
                        <Descriptions.Item label="反应物">{selectedLink.source}</Descriptions.Item>
                        <Descriptions.Item label="产物">{selectedLink.target}</Descriptions.Item>
                        <Descriptions.Item label="算符">
                            <Tag color="orange">{selectedLink.operator}</Tag>
                        </Descriptions.Item>
                        <Descriptions.Item label="反应能">
                            {selectedLink.energy ? `${selectedLink.energy.toFixed(2)} kcal/mol` : '-'}
                        </Descriptions.Item>
                        {Object.entries(selectedLink.properties || {}).map(([key, value]) => (
                            <Descriptions.Item key={key} label={key}>
                                {String(value)}
                            </Descriptions.Item>
                        ))}
                    </Descriptions>
                )}
            </Drawer>
        </div>
    );
};

export default ReactionNetworkViewer;
