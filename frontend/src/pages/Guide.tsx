/**
 * Molyte 科学计算指南
 * 专注于计算原理、功能意义及操作最佳实践
 */
import { useNavigate } from 'react-router-dom';
import {
  Typography,
  Card,
  Row,
  Col,
  Collapse,
  Button,
  Divider,
  Space,
  Tabs,
  Tag,
  Steps,
  Alert,
  Tooltip,
  Timeline,
  Image,
  Table
} from 'antd';
import {
  HomeOutlined,
  ThunderboltOutlined,
  ExperimentOutlined,
  RocketOutlined,
  LineChartOutlined,
  BulbOutlined,
  DatabaseOutlined,
  FireOutlined,
  ClusterOutlined,
  ApartmentOutlined,
  CheckCircleOutlined,
  InfoCircleOutlined,
  SyncOutlined,
  SettingOutlined,
  DashboardOutlined
} from '@ant-design/icons';
import { useAuthStore } from '../stores/authStore';
import './Guide.css';

const { Title, Paragraph, Text } = Typography;
const { TabPane } = Tabs;
const { Panel } = Collapse;

export default function Guide() {
  const navigate = useNavigate();
  const { isAuthenticated } = useAuthStore();

  return (
    <div className="guide-container">
      {/* Background Elements from Home */}
      <div className="guide-background">
        <div className="bg-particles"></div>
        <div className="bg-gradient"></div>
        <div className="bg-grid"></div>
        <div className="floating-shapes">
          <div className="shape shape-1"></div>
          <div className="shape shape-2"></div>
          <div className="shape shape-3"></div>
        </div>
      </div>

      {/* Header */}
      <header className="guide-header">
        <div className="guide-header-content">
          <div className="guide-logo" onClick={() => navigate('/')}>
            <div className="guide-logo-icon">
              <ThunderboltOutlined />
            </div>
            <span className="guide-logo-text">Molyte 计算指南</span>
          </div>
          <nav className="guide-nav">
            <a href="#md-simulation">MD 溶剂化</a>
            <a href="#reaction-network">SEI 反应网络</a>
            <a href="#md-simulation">MD 溶剂化</a>
            <a href="#reaction-network">SEI 反应网络</a>
            <a href="#qc-calc">量子化学</a>
            <a href="#post-process">数据后处理</a>
            <a href="#best-practices">最佳实践</a>
          </nav>
          <div className="guide-header-actions">
            <Button type="text" onClick={() => navigate('/')} className="guide-back-btn">
              返回首页
            </Button>
            {isAuthenticated ? (
              <Button type="primary" onClick={() => navigate('/workspace/dashboard')} className="guide-action-btn">
                开始计算
              </Button>
            ) : (
              <Button type="primary" onClick={() => navigate('/login')} className="guide-action-btn">
                登录平台
              </Button>
            )}
          </div>
        </div>
      </header>

      <main className="guide-main">
        {/* Intro */}
        <section className="guide-hero">
          <div className="guide-hero-content">
            <Title level={1} className="guide-hero-title">电池材料计算背后的科学</Title>
            <Paragraph className="guide-hero-subtitle">
              从微观尺度理解电解液性能：原理、方法与结果解读
            </Paragraph>
          </div>
        </section>

        {/* Module 1: Liquid Electrolyte MD */}
        <section id="md-simulation" className="guide-section">
          <div className="section-header">
            <div className="section-icon"><ExperimentOutlined /></div>
            <Title level={2}>溶元调配与分子动力学 (MD)</Title>
          </div>

          <Card className="guide-card">
            <Tabs defaultActiveKey="principle" size="large">
              <TabPane tab={<span><BulbOutlined /> 计算原理</span>} key="principle">
                <Row gutter={24}>
                  <Col span={14}>
                    <Title level={4}>经典力场模拟 (OPLS-AA)</Title>
                    <Paragraph>
                      本模块的核心是通过分子动力学（MD）模拟电解液在微观尺度的热运动。我们采用 <Text strong>OPLS-AA (Optimized Potentials for Liquid Simulations)</Text> 全原子力场，
                      它包含键合相互作用（键长、键角、二面角）和非键合相互作用（范德华力、库仑力）。
                    </Paragraph>
                    <div className="formula-box">
                      E_total = E_bond + E_angle + E_torsion + E_vdw + E_coulomb
                    </div>
                    <Paragraph>
                      <ul>
                        <li><Text strong>E_vdw</Text>: 采用 Lennard-Jones 12-6 势能函数描述原子间的排斥与吸引。</li>
                        <li><Text strong>E_coulomb</Text>: 通过 CM1A/LBCC 电荷模型计算静电相互作用，通过 PME 算法处理长程作用。</li>
                      </ul>
                    </Paragraph>
                    <Alert
                      message="为什么要算 MD？"
                      description="实验只能测得宏观性质（如电导率），但无法直接观测锂离子的微观溶剂化结构（Solvation Structure）。MD 模拟如同显微镜，让我们'看到'锂离子被哪些溶剂分子包围，这直接决定了去溶剂化能（影响低温性能）和界面成膜机制。"
                      type="info"
                      showIcon
                    />
                  </Col>
                  <Col span={10}>
                    <div className="theory-card">
                      <Title level={5}>模拟流程详解</Title>
                      <Steps direction="vertical" size="small" current={1} items={[
                        { title: 'Packmol 建模', description: '根据设定的摩尔比，在模拟盒子中随机填充阳离子、阴离子和溶剂分子。' },
                        { title: '能量最小化', description: '消除原子重叠，防止初始构型能量过高导致体系崩溃 ("爆炸")。' },
                        { title: 'NPT 平衡 (等温等压)', description: '关键步骤！调整盒子体积以达到实验密度。通常需要 1ns 以上。' },
                        { title: 'NVT 采样 (等温等容)', description: '在平衡密度下收集轨迹数据，用于计算扩散系数和 RDF。' }
                      ]} />
                    </div>
                  </Col>
                </Row>
              </TabPane>

              <TabPane tab={<span><LineChartOutlined /> 结果解读</span>} key="results">
                <Row gutter={24}>
                  <Col span={12}>
                    <Card title="1. 径向分布函数 (RDF, g(r))" bordered={false} className="feature-item">
                      <Row gutter={24} align="middle">
                        <Col span={12}>
                          <Paragraph>
                            RDF 反映了以 Li+ 为中心，距离 r 处出现氧原子（来自溶剂或阴离子）的概率密度。
                          </Paragraph>
                          <ul>
                            <li><Text strong>第一峰位置</Text>: 代表 Li-O 键长（通常约 1.8-2.1 Å）。</li>
                            <li><Text strong>配位数 (CN)</Text>: 对第一峰积分。CN ≈ 4 说明是典型的四配位结构；CN &gt; 5 可能存在溶剂共享或聚集体。</li>
                            <li><Text strong>意义</Text>: 判断溶剂的竞争能力。如果阴离子的峰很高，说明形成了接触离子对 (CIP) 或聚集体 (AGG)。</li>
                          </ul>
                        </Col>
                        <Col span={12} style={{ textAlign: 'center' }}>
                          <Image
                            src="/assets/rdf_chart_blue.png"
                            alt="RDF Chart Visualization"
                            style={{ borderRadius: 8, border: '1px solid rgba(255,255,255,0.1)', maxWidth: '100%', maxHeight: 300 }}
                          />
                        </Col>
                      </Row>
                    </Card>
                  </Col>
                  <Col span={12}>
                    <Card title="2. 均方位移 (MSD) 与自扩散系数" bordered={false} className="feature-item">
                      <Row gutter={24} align="middle">
                        <Col span={12}>
                          <Paragraph>
                            MSD 描述粒子随时间扩散的距离平方。其斜率直接正比于扩散系数 D。
                          </Paragraph>
                          <div className="formula-box" style={{ margin: '16px 0' }}>D = MSD / (6 * t)</div>
                          <ul>
                            <li><Text strong>斜率越陡</Text>: 扩散越快，离子电导率越高。</li>
                            <li><Text strong>非线性区</Text>: 模拟初期（弹道区）的数据不可用，必须取 MSD 线性段计算斜率。</li>
                            <li><Text strong>机理</Text>: 若锂离子的 D 与溶剂的 D 相近，说明是 Vehicular 机制。</li>
                          </ul>
                        </Col>
                        <Col span={12} style={{ textAlign: 'center' }}>
                          <Image
                            src="/assets/msd_chart_blue.png"
                            alt="MSD Chart Visualization"
                            style={{ borderRadius: 8, border: '1px solid rgba(255,255,255,0.1)', maxWidth: '100%', maxHeight: 300 }}
                          />
                        </Col>
                      </Row>
                    </Card>
                  </Col>
                </Row>
                <Row gutter={24} style={{ marginTop: 24 }}>
                  <Col span={24}>
                    <Card title="3. 溶剂化结构可视化 (3D Visualization)" bordered={false} className="feature-item">
                      <Row gutter={24} align="middle">
                        <Col span={10}>
                          <Paragraph>
                            高精度的 3D 渲染图帮助直观理解离子的配位环境。紫色为中心锂离子，周围环绕着 EC/DMC 溶剂分子，形成第一配位壳层。
                          </Paragraph>
                          <ul>
                            <li><Text strong>配位壳层</Text>: 观察第一配位层的溶剂分子数量和排列方式</li>
                            <li><Text strong>离子对</Text>: 识别阴离子是否进入配位层 (CIP/AGG 结构)</li>
                            <li><Text strong>空间结构</Text>: 理解离子-溶剂、离子-阴离子的空间相互作用模式</li>
                          </ul>
                        </Col>
                        <Col span={14} style={{ textAlign: 'center' }}>
                          <Image
                            src="/assets/solvation_structure.png"
                            alt="Li+ Solvation Structure"
                            style={{ borderRadius: 8, border: '1px solid rgba(255,255,255,0.1)', maxWidth: '100%', maxHeight: 350 }}
                          />
                        </Col>
                      </Row>
                    </Card>
                  </Col>
                </Row>
              </TabPane>

              <TabPane tab={<span><SettingOutlined /> 操作指南</span>} key="guide">
                <Alert
                  message="完整操作流程"
                  description="从配方设计到结果分析的全流程指导，帮助您快速上手并获得高质量结果。"
                  type="success"
                  showIcon
                  style={{ marginBottom: 24 }}
                />
                <Collapse accordion defaultActiveKey={['1']}>
                  <Panel header="📋 Step 1: 配方设计与参数设置" key="1">
                    <Title level={5}>1.1 组分选择</Title>
                    <Paragraph>
                      <ul>
                        <li><Text strong>阳离子</Text>: 通常选择 Li+（锂离子电池）或 Na+（钠离子电池）</li>
                        <li><Text strong>阴离子</Text>:
                          <br />• PF6⁻, TFSI⁻ (常规电解液)
                          <br />• FSI⁻, BOB⁻ (低温/高压专用)
                          <br />• 建议盐浓度: 1.0-1.2 M (标准), 3-5 M (高浓度电解液 HCE)
                        </li>
                        <li><Text strong>溶剂</Text>:
                          <br />• EC+DMC (经典配方，循环稳定)
                          <br />• EC+EMC (低温性能优)
                          <br />• FEC 添加剂 (5-10%, 改善 SEI)
                        </li>
                      </ul>
                    </Paragraph>
                    <Title level={5}>1.2 温度与时长策略</Title>
                    <Table
                      size="small"
                      dataSource={[
                        { key: '1', scenario: '快速筛选', temp: '298 K', npt: '0.5 ns', nvt: '1 ns', purpose: '定性对比配位结构' },
                        { key: '2', scenario: '标准计算', temp: '298 K', npt: '2 ns', nvt: '5 ns', purpose: '发表级数据质量' },
                        { key: '3', scenario: '高温加速', temp: '333-353 K', npt: '1 ns', nvt: '3 ns', purpose: '高粘度体系快速收敛' },
                        { key: '4', scenario: '低温性能', temp: '253 K', npt: '5 ns', nvt: '20 ns', purpose: '低温电池设计' },
                      ]}
                      columns={[
                        { title: '应用场景', dataIndex: 'scenario', key: 'scenario' },
                        { title: '温度', dataIndex: 'temp', key: 'temp' },
                        { title: 'NPT 平衡', dataIndex: 'npt', key: 'npt' },
                        { title: 'NVT 采样', dataIndex: 'nvt', key: 'nvt' },
                        { title: '用途', dataIndex: 'purpose', key: 'purpose' },
                      ]}
                      pagination={false}
                    />
                  </Panel>

                  <Panel header="⚙️ Step 2: 任务提交与监控" key="2">
                    <Paragraph>
                      <Text strong>提交前检查清单:</Text>
                      <ul>
                        <li>✓ 确认盒子尺寸: 建议 &gt;40 Å (避免周期性边界效应)</li>
                        <li>✓ 分子数量: 总原子数 3000-8000 (平衡精度与速度)</li>
                        <li>✓ 时间步长: 默认 1 fs (刚性键) 或 0.5 fs (柔性全原子)</li>
                      </ul>
                    </Paragraph>
                    <Paragraph>
                      <Text strong>任务监控:</Text>
                      <br />• 在 "任务列表" 页面实时查看进度
                      <br />• NPT 阶段: 关注密度收敛曲线 (目标: ±0.02 g/cm³)
                      <br />• NVT 阶段: 检查温度波动 (应 &lt;5 K)
                    </Paragraph>
                  </Panel>

                  <Panel header="📊 Step 3: 结果分析与判断" key="3">
                    <Title level={5}>3.1 RDF 分析</Title>
                    <Paragraph>
                      <Text strong>关键指标:</Text>
                      <ul>
                        <li>第一峰位置: Li-O 距离应在 1.9-2.1 Å</li>
                        <li>配位数 CN: 积分第一峰至第一谷底
                          <br />• CN = 4: 典型四配位 (EC/DMC 体系)
                          <br />• CN = 5-6: 高浓度或强配位溶剂
                          <br />• CN &lt; 4: 可能存在离子对 (CIP/AGG)
                        </li>
                        <li>阴离子峰强度: 峰高 &gt;1.5 说明阴离子参与配位 (有利于 SEI)</li>
                      </ul>
                    </Paragraph>
                    <Title level={5}>3.2 MSD 与扩散系数</Title>
                    <Paragraph>
                      <Text strong>数据处理:</Text>
                      <ul>
                        <li>舍弃前 20% 数据 (弹道区)</li>
                        <li>线性拟合 MSD vs t 曲线 (R² &gt; 0.98)</li>
                        <li>计算 D = slope / 6</li>
                        <li>典型值参考:
                          <br />• Li+ 在 EC/DMC (1M LiPF6): D ≈ 2-5 × 10⁻⁶ cm²/s
                          <br />• 若 D &lt; 1 × 10⁻⁶: 体系过于粘稠或未平衡
                        </li>
                      </ul>
                    </Paragraph>
                  </Panel>

                  <Panel header="⚠️ 常见问题与解决" key="4">
                    <Paragraph>
                      <Text strong type="danger">问题 1: 密度不收敛</Text>
                      <br />→ 延长 NPT 时间至 5 ns
                      <br />→ 检查初始构型是否合理 (Packmol 可能生成重叠)
                    </Paragraph>
                    <Paragraph>
                      <Text strong type="danger">问题 2: MSD 曲线不是直线</Text>
                      <br />→ NVT 时间不足，延长至 10 ns
                      <br />→ 体系可能存在相分离或结晶
                    </Paragraph>
                    <Paragraph>
                      <Text strong type="danger">问题 3: 扩散系数异常大/小</Text>
                      <br />→ 检查温度设置是否正确
                      <br />→ 验证力场参数 (OPLS-AA 对某些新型溶剂可能不准)
                    </Paragraph>
                  </Panel>
                </Collapse>
              </TabPane>
            </Tabs>
          </Card>
        </section>

        {/* Module 2: Reaction Network */}
        <section id="reaction-network" className="guide-section">
          <div className="section-header">
            <div className="section-icon"><FireOutlined /></div>
            <Title level={2}>溶盐反应网络 (Reaction Network)</Title>
          </div>

          <Alert
            message="预测 SEI/CEI 膜成分的关键技术"
            description="本模块通过模拟分子在电压、温度和活性自由基驱动下的化学演化，预测界面膜的组成成分。"
            type="warning"
            showIcon
            style={{ marginBottom: 24 }}
          />

          <Card className="guide-card">
            <Row gutter={[24, 24]}>
              <Col xs={24} lg={12}>
                <div className="theory-card" style={{ height: '100%' }}>
                  <Title level={4}>RSNet 算法原理</Title>
                  <Paragraph>
                    RSNet (Reaction Space Network) 是一种基于规则的反应生成算法。不同于昂贵的 AIMD（从头算分子动力学），
                    RSNet 通过预定义的<Text strong>反应算符 (Operators)</Text> 快速遍历化学空间。
                  </Paragraph>
                  <Divider />
                  <Title level={5}>算符激活机制</Title>
                  <Paragraph>
                    算符是否“激活”取决于环境驱动力：
                  </Paragraph>
                  <ul>
                    <li><Tag color="red">电压驱动</Tag>: 当电势低于 LUMO（负极还原）或高于 HOMO（正极氧化）时，激活电子转移算符。</li>
                    <li><Tag color="orange">热驱动</Tag>: 高温可激活高能垒的键断裂反应（如 C-F 键断裂）。</li>
                    <li><Tag color="blue">自由基驱动</Tag>: 一旦生成高活性自由基（如 H•, F•），立即激活链式反应算符。</li>
                  </ul>
                </div>
              </Col>

              <Col xs={24} lg={12}>
                <div className="feature-item" style={{ height: '100%' }}>
                  <Title level={4}>结果应用：该看什么？</Title>
                  <Timeline
                    items={[
                      {
                        color: 'green',
                        children: (
                          <>
                            <Text strong>第一代产物 (Gen 1):</Text>
                            <p>通常是直接的氧化/还原产物。例如 EC 得到双电子还原生成 LEDC (Li2CO3前体)。这决定了 SEI 的主要骨架。</p>
                          </>
                        ),
                      },
                      {
                        color: 'red',
                        children: (
                          <>
                            <Text strong>终端产物 (Thermodynamic Sinks):</Text>
                            <p>只有热力学极度稳定的产物才会累积，如 <Text strong>LiF, Li2CO3, Li2O</Text>。这些无机成分含量越高，通常意味着界面膜越致密、更稳定。</p>
                          </>
                        ),
                      },
                      {
                        color: 'gray',
                        children: (
                          <>
                            <Text strong>气体释放:</Text>
                            <p>关注路径中是否生成 CO2, C2H4, H2 等小分子气体。这是电池胀气和安全性评估的重要指标。</p>
                          </>
                        ),
                      },
                    ]}
                  />
                </div>
              </Col>
            </Row>

            <div style={{ marginTop: 24 }}>
              <Title level={4}><SettingOutlined /> 完整操作流程</Title>
              <Collapse accordion>
                <Panel header="🎯 Step 1: 任务配置" key="1">
                  <Title level={5}>1.1 电极类型选择</Title>
                  <Table
                    size="small"
                    dataSource={[
                      { key: '1', electrode: 'Anode (负极)', voltage: '0.01-0.5 V', env: 'Li+ + e⁻ 还原', target: 'SEI 膜成分预测' },
                      { key: '2', electrode: 'Cathode (正极)', voltage: '4.3-4.8 V', env: '氧化环境 + 金属离子催化', target: 'CEI 膜、气体释放' },
                      { key: '3', electrode: 'Bulk (本体)', voltage: '任意', env: '热分解', target: '高温稳定性' },
                    ]}
                    columns={[
                      { title: '电极', dataIndex: 'electrode', key: 'electrode' },
                      { title: '电压范围', dataIndex: 'voltage', key: 'voltage' },
                      { title: '反应环境', dataIndex: 'env', key: 'env' },
                      { title: '研究目标', dataIndex: 'target', key: 'target' },
                    ]}
                    pagination={false}
                  />
                  <Title level={5} style={{ marginTop: 16 }}>1.2 关键参数设置</Title>
                  <Paragraph>
                    <ul>
                      <li><Text strong>代数 (Generations)</Text>:
                        <br />• Gen 1-2: 初级产物 (单体、自由基)
                        <br />• Gen 3: 二聚体、环状产物 (推荐)
                        <br />• Gen 4+: 计算量爆炸，仅用于特殊研究
                      </li>
                      <li><Text strong>温度</Text>:
                        <br />• 298 K: 常温循环
                        <br />• 333 K: 高温存储 (加速老化)
                        <br />• 353-373 K: 热失控预警
                      </li>
                      <li><Text strong>溶剂/盐选择</Text>:
                        <br />• 必须与 MD 模拟保持一致
                        <br />• 可添加添加剂 (如 VC, FEC) 研究协同效应
                      </li>
                    </ul>
                  </Paragraph>
                </Panel>

                <Panel header="🔬 Step 2: 结果解读策略" key="2">
                  <Title level={5}>2.1 反应网络图分析</Title>
                  <Paragraph>
                    <Text strong>优先关注:</Text>
                    <ul>
                      <li>节点大小 = 产物浓度 → 找最大的几个节点</li>
                      <li>边的粗细 = 反应速率 → 主要反应路径</li>
                      <li>颜色编码:
                        <br />• 绿色: 稳定产物 (LiF, Li2CO3)
                        <br />• 红色: 气体 (CO2, C2H4) - 安全隐患
                        <br />• 黄色: 中间体 (可能继续反应)
                      </li>
                    </ul>
                  </Paragraph>
                  <Title level={5}>2.2 关键产物识别</Title>
                  <Paragraph>
                    <Text strong>SEI 膜优质指标:</Text>
                    <ul>
                      <li>✓ LiF 含量 &gt;30%: 高机械强度</li>
                      <li>✓ Li2CO3 / LEDC: 离子传导性</li>
                      <li>✓ 有机聚合物 (PVDF 类): 柔韧性</li>
                      <li>✗ 气体产率 &gt;10%: 胀气风险</li>
                      <li>✗ 可溶性产物: 持续消耗电解液</li>
                    </ul>
                  </Paragraph>
                </Panel>

                <Panel header="📈 Step 3: 数据导出与应用" key="3">
                  <Paragraph>
                    <Text strong>可导出数据:</Text>
                    <ul>
                      <li>反应网络 JSON (用于机器学习)</li>
                      <li>产物分布 CSV (绘制柱状图)</li>
                      <li>SMILES 列表 (后续 QC 计算)</li>
                    </ul>
                  </Paragraph>
                  <Paragraph>
                    <Text strong>典型应用场景:</Text>
                    <ul>
                      <li>对比不同电解液配方的 SEI 成分差异</li>
                      <li>筛选低产气配方 (电动汽车安全)</li>
                      <li>预测添加剂的作用机制</li>
                    </ul>
                  </Paragraph>
                </Panel>

                <Panel header="⚠️ 常见错误与规避" key="4">
                  <Paragraph>
                    <Text type="danger" strong>错误 1: 电压设置不当</Text>
                    <br />→ SEI 研究必须 &lt;0.8 V (否则无还原反应)
                    <br />→ CEI 研究必须 &gt;4.2 V (否则无氧化反应)
                  </Paragraph>
                  <Paragraph>
                    <Text type="danger" strong>错误 2: 代数过高</Text>
                    <br />→ Gen 5+ 会产生数千种产物，难以分析
                    <br />→ 且高代产物通常浓度极低 (&lt;0.1%)
                  </Paragraph>
                  <Paragraph>
                    <Text type="danger" strong>错误 3: 忽略溶剂效应</Text>
                    <br />→ 不同溶剂会显著改变反应路径
                    <br />→ 必须用实际电解液配方，不能只算单一溶剂
                  </Paragraph>
                </Panel>
              </Collapse>
            </div>
          </Card>
        </section>

        {/* Module 3: Quantum Chemistry */}
        <section id="qc-calc" className="guide-section">
          <div className="section-header">
            <div className="section-icon"><ThunderboltOutlined /></div>
            <Title level={2}>量子化学计算 (QC)</Title>
          </div>
          <Card className="guide-card">
            <Paragraph>
              通过求解薛定谔方程（基于 DFT 近似），获得分子的电子结构。这是理解电解液氧化还原稳定性的基石。
            </Paragraph>
            <Row gutter={24}>
              <Col span={8}>
                <Card title="HOMO & LUMO" size="small" className="feature-item">
                  <div style={{ textAlign: 'center', marginBottom: 16 }}>
                    <div style={{ height: 4, width: '100%', background: '#ff4d4f', marginBottom: 4 }}></div>
                    <Text type="secondary">LUMO (易得电子)</Text>
                    <div style={{ margin: '12px 0', fontSize: 24 }}>⚡ Energy Gap</div>
                    <Text type="secondary">HOMO (易失电子)</Text>
                    <div style={{ height: 4, width: '100%', background: '#1890ff', marginTop: 4 }}></div>
                  </div>
                  <Paragraph>
                    <Text strong>HOMO 能级低</Text> &rarr; 抗氧化能力强 &rarr; 适合高压正极。
                    <br />
                    <Text strong>LUMO 能级高</Text> &rarr; 抗还原能力强 &rarr; 在负极稳定。
                  </Paragraph>
                </Card>
              </Col>
              <Col span={8}>
                <Card title="静电势 (ESP)" size="small" className="feature-item">
                  <Paragraph>
                    ESP 映射图显示了分子表面的电荷分布。
                    <br /><br />
                    <Text strong>红色区域 (负电)</Text>: 易吸引 Li+，是配位位点（如 EC 的羰基氧）。
                    <br />
                    <Text strong>蓝色区域 (正电)</Text>: 易受亲核试剂攻击（如 F- 进攻）。
                  </Paragraph>
                </Card>
              </Col>
              <Col span={8}>
                <Card title="计算建议" size="small" className="theory-card">
                  <Title level={5}>方法选择</Title>
                  <p><Text code>B3LYP/6-31G*</Text>:
                    <br />性价比之王，适合快速筛选几百个分子。</p>
                  <p><Text code>M06-2X/def2-TZVP</Text>:
                    <br />描述弱相互作用（如 π-π 堆积）更准确，计算反应能时推荐使用。</p>
                  <p><Text code>wB97X-D</Text>:
                    <br />包含色散校正，精度最高，但耗时较长。</p>
                </Card>
              </Col>
            </Row>
          </Card>
        </section>

        {/* Module 4: Post-Process */}
        <section id="post-process" className="guide-section">
          <div className="section-header">
            <div className="section-icon"><ClusterOutlined /></div>
            <Title level={2}>数据后处理 (Post-Processing)</Title>
          </div>
          <Card className="guide-card">
            <Row gutter={[24, 24]}>
              <Col xs={24} lg={12}>
                <Title level={4}>去溶剂化能 (Desolvation Energy)</Title>
                <div style={{ marginBottom: 16 }}>
                  <Image
                    src="/assets/desolvation_process.png"
                    fallback="/assets/rdf_chart_blue.png"
                    alt="Desolvation Process Diagram"
                    style={{ borderRadius: 8, border: '1px solid rgba(255,255,255,0.1)', maxWidth: '100%', maxHeight: 200, objectFit: 'cover' }}
                  />
                </div>
                <Paragraph>
                  <Text strong>物理意义:</Text>
                  锂离子在嵌入电极（如石墨或硅负极）之前，必须先摆脱其溶剂化壳层的束缚。这个"脱衣"过程的能量代价即为去溶剂化能。
                  <br />
                  去溶剂化通常是低温下电池内阻的主要来源（占总阻抗的 60% 以上）。去溶剂化能越低，低温充电性能越好。
                </Paragraph>
                <div className="theory-card">
                  <Title level={5}>两种计算模式详解</Title>
                  <Paragraph>
                    <Text strong>1. Stepwise (逐步脱溶剂机制)</Text>:
                    <br />
                    模拟锂离子逐个丢弃溶剂分子的过程。
                    <div className="formula-box" style={{ margin: '8px 0', fontSize: 12 }}>
                      Li(Sol)₄ &rarr; Li(Sol)₃ + Sol (&Delta;G₁) <br />
                      Li(Sol)₃ &rarr; Li(Sol)₂ + Sol (&Delta;G₂)
                    </div>
                    通常最后一步（如配位数 2&rarr;1 或 3&rarr;2）能垒最高，是决速步。
                  </Paragraph>
                  <Paragraph>
                    <Text strong>2. Full (完全脱溶剂)</Text>:
                    <br />
                    Li<sup>+</sup>(Solvent)<sub>n</sub> &rarr; Li<sup>+</sup> + n &times; Solvent
                    <br />
                    计算将离子完全剥离裸露所需的总能量。这反映了溶剂与离子的本征结合强度，常用于筛选新型高压溶剂（强结合=高氧化稳定性）。
                  </Paragraph>
                </div>
              </Col>
              <Col xs={24} lg={12}>
                <Card title="氧化还原电位与重组能 (Redox & Reorganization)" bordered={false} className="feature-item">
                  <div style={{ marginBottom: 16, textAlign: 'center' }}>
                    <Image
                      src="/assets/redox_energy_diagram.png"
                      alt="Marcus Theory Energy Diagram"
                      style={{ borderRadius: 8, border: '1px solid rgba(255,255,255,0.1)', maxWidth: '100%', maxHeight: 200, objectFit: 'contain' }}
                    />
                  </div>
                  <div className="theory-card">
                    <Paragraph>
                      <Text strong>1. 氧化/还原电位 (Redox Potential)</Text>:
                      <br />
                      我们采用最高精度的 <Text strong>G4MP2</Text> 或 <Text strong>M06-2X</Text> 方法计算分子的绝热氧化/还原电位。
                      <div className="formula-box" style={{ margin: '8px 0', fontSize: 13 }}>
                        E_ox = (G(M⁺) - G(M)) / F - 1.4V
                      </div>
                      相比于简单的 HOMO/LUMO 能级，该方法考虑了溶剂化效应和构型弛豫，预测值与 CV 实验吻合度更高（误差 &lt; 0.2V）。
                    </Paragraph>
                    <Paragraph>
                      <Text strong>2. 重组能 (Reorganization Energy, &lambda;)</Text>:
                      <br />
                      定义为电子转移过程中，分子因几何构型变化而消耗的能量。它是 <Text strong>Marcus 电子转移理论</Text> 的核心参数。
                      <br />
                      <ul>
                        <li><Text strong>&lambda; 大 (如 EC, &gt;0.8 eV)</Text>: 结构变化剧烈，电子转移慢，阻抗高。</li>
                        <li><Text strong>&lambda; 小 (如 芳香族分子, &lt;0.2 eV)</Text>: 结构刚性，电子转移极快，适合做添加剂或快充溶剂。</li>
                      </ul>
                    </Paragraph>
                  </div>
                </Card>
              </Col>
            </Row>
          </Card>
        </section>

        {/* Footer */}
        <footer className="guide-footer" style={{ textAlign: 'center', marginTop: 80, paddingBottom: 40 }}>
          Molyte Team &copy; 2025
        </footer>
      </main>
    </div>
  );
}
