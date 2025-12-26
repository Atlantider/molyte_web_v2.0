# RSNet Simple API - 快速入门指南

## 简介

RSNet Simple API 是一个简化的反应网络生成接口，让您可以轻松地：
- 输入分子的SMILES表示
- 自动发现可能的反应类型
- 计算反应能量
- 生成可视化的反应网络图

## 安装

确保已安装所有依赖：

```bash
cd d:\codex\rsnet-main\rsnet-main
pip install -r requirements.txt
```

## 快速开始

### 方法1：Python脚本

```python
from rsnet_simple_api import generate_reaction_network

# 生成反应网络
result = generate_reaction_network(
    smiles_list=['C1COC(=O)O1', '[Li+]', 'F[P-](F)(F)(F)(F)F'],
    temperature=300.0,
    max_generations=3
)

# 查看结果
print(f"生成了 {len(result['molecules'])} 个分子")
print(f"发现了 {len(result['reactions'])} 个反应")
print(f"可视化保存在: {result['visualization_path']}")
```

### 方法2：命令行

```bash
# 基本用法
python rsnet_simple_api.py --smiles "CCO" "C=C"

# 指定参数
python rsnet_simple_api.py --smiles "C1COC(=O)O1" "[Li+]" \
    --temperature 300 \
    --electrode anode \
    --max-generations 3

# 从文件读取SMILES
python rsnet_simple_api.py --input molecules.txt --output ./results
```

### 方法3：交互式示例

```bash
python examples_simple_api.py
```

## API参数说明

### `generate_reaction_network()` 函数

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `smiles_list` | List[str] | **必需** | SMILES字符串列表 |
| `temperature` | float | 300.0 | 温度 (K) |
| `electrode_type` | str | 'anode' | 电极类型 ('anode' 或 'cathode') |
| `voltage` | float | 0.1 | 电压 (V) |
| `max_generations` | int | 3 | 最大代数 |
| `max_species` | int | 50 | 最大分子数 |
| `energy_cutoff` | float | 80.0 | 能量截断值 (kcal/mol) |
| `visualize` | bool | True | 是否生成可视化 |
| `save_results` | bool | True | 是否保存结果到文件 |
| `output_dir` | str | './rsnet_output' | 输出目录 |

### 返回值

返回一个字典，包含：

```python
{
    'network': ReactionNetwork,          # 网络对象
    'molecules': List[Molecule],         # 所有分子列表
    'reactions': List[Reaction],         # 所有反应列表
    'statistics': {                      # 统计信息
        'num_molecules': int,
        'num_reactions': int,
        'max_generation': int
    },
    'generation_time': float,            # 生成时间（秒）
    'environment': {...},                # 环境参数
    'visualization_path': str,           # 可视化文件路径
    'json_path': str                     # JSON结果文件路径
}
```

## 使用示例

### 示例1：锂离子电池电解液

```python
from rsnet_simple_api import generate_reaction_network

result = generate_reaction_network(
    smiles_list=[
        'C1COC(=O)O1',           # EC (碳酸乙烯酯)
        '[Li+]',                  # 锂离子
        'F[P-](F)(F)(F)(F)F'     # PF6- 阴离子
    ],
    temperature=300.0,
    electrode_type='anode',
    voltage=0.1,
    max_generations=3
)

print(f"发现 {len(result['reactions'])} 个反应")
```

### 示例2：有机分子反应

```python
result = generate_reaction_network(
    smiles_list=['CCO', 'C=C'],  # 乙醇 + 乙烯
    temperature=400.0,
    max_generations=3
)

# 查看反应
for rxn in result['reactions']:
    reactants = ' + '.join([r.smiles for r in rxn.reactants])
    products = ' + '.join([p.smiles for p in rxn.products])
    print(f"{reactants} → {products}")
```

### 示例3：程序化访问

```python
result = generate_reaction_network(
    smiles_list=['C1CC1'],  # 环丙烷
    temperature=500.0,
    visualize=False,
    save_results=False
)

# 直接访问网络对象
network = result['network']

# 分析分子
for mol in result['molecules']:
    print(f"{mol.name}: {mol.smiles}")

# 分析反应类型
reaction_types = {}
for rxn in result['reactions']:
    op = rxn.operator_name if hasattr(rxn, 'operator_name') else 'unknown'
    reaction_types[op] = reaction_types.get(op, 0) + 1

print("反应类型分布:", reaction_types)
```

## 输出文件

### 可视化图片 (PNG)

自动生成的网络可视化图，包含：
- 左侧：反应网络拓扑图
  - 节点：分子（按代数着色）
  - 边：反应（箭头表示方向）
- 右侧：统计信息面板

文件名格式：`network_YYYYMMDD_HHMMSS.png`

### JSON结果文件

包含完整的网络信息：

```json
{
  "statistics": {
    "num_molecules": 15,
    "num_reactions": 8,
    "max_generation": 3
  },
  "molecules": [
    {
      "name": "mol_1",
      "smiles": "CCO",
      "generation": 0
    }
  ],
  "reactions": [
    {
      "reactants": ["CCO"],
      "products": ["CC=O"],
      "operator": "dehydrogenation",
      "energy": -10.5
    }
  ]
}
```

文件名格式：`results_YYYYMMDD_HHMMSS.json`

## 常见问题

### Q: 如何提高反应发现的数量？

A: 可以调整以下参数：
- 增加 `max_generations`（更多代数）
- 增加 `max_species`（允许更多分子）
- 增加 `energy_cutoff`（允许更高能量的反应）
- 提高 `temperature`（更多热激活反应）

### Q: 如何加快生成速度？

A: 可以：
- 减少 `max_generations`
- 减少 `max_species`
- 降低 `energy_cutoff`

### Q: 支持哪些类型的分子？

A: 支持任何有效的SMILES表示的分子，包括：
- 有机分子
- 离子（如 `[Li+]`, `[Cl-]`）
- 自由基
- 金属配合物

### Q: 能量计算准确吗？

A: 系统使用半经验方法（xTB）进行能量估算，适合快速筛选。如需高精度，建议使用DFT等方法进一步验证。

## 高级用法

### 自定义环境条件

```python
from rsnet.core.environment import Environment

# 创建自定义环境
env = Environment(
    temperature=350.0,
    electrode_type='cathode',
    voltage=4.2,
    li_activity=0.5,
    solvent='DMC'
)

# 使用自定义环境需要直接使用API类
from rsnet_simple_api import RSNetSimpleAPI

api = RSNetSimpleAPI()
result = api.generate_reaction_network(
    smiles_list=['CCO'],
    temperature=env.temperature,
    electrode_type=env.electrode_type,
    voltage=env.voltage
)
```

### 批量处理

```python
import pandas as pd
from rsnet_simple_api import generate_reaction_network

# 从CSV读取分子列表
df = pd.read_csv('molecules.csv')

results = []
for idx, row in df.iterrows():
    result = generate_reaction_network(
        smiles_list=[row['smiles']],
        temperature=row['temperature'],
        max_generations=2,
        output_dir=f'./results/mol_{idx}'
    )
    results.append({
        'name': row['name'],
        'num_reactions': len(result['reactions']),
        'num_molecules': len(result['molecules'])
    })

# 保存汇总结果
summary_df = pd.DataFrame(results)
summary_df.to_csv('summary.csv', index=False)
```

## 运行测试

```bash
# 运行所有测试
python -m pytest test_simple_api.py -v

# 运行特定测试
python -m pytest test_simple_api.py::TestRSNetSimpleAPI::test_basic_network_generation -v
```

## 故障排除

### 问题：导入错误

```
ModuleNotFoundError: No module named 'rsnet'
```

**解决方案**：确保在正确的目录运行，并且rsnet包在Python路径中：

```bash
cd d:\codex\rsnet-main\rsnet-main
export PYTHONPATH=$PYTHONPATH:$(pwd)  # Linux/Mac
$env:PYTHONPATH += ";$(pwd)"          # Windows PowerShell
```

### 问题：可视化不显示中文

**解决方案**：安装支持中文的字体，或在代码中指定字体：

```python
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['SimHei']  # Windows
```

## 更多资源

- 完整示例：运行 `python examples_simple_api.py`
- API文档：查看 `rsnet_simple_api.py` 中的docstring
- 测试用例：参考 `test_simple_api.py`

## 联系与反馈

如有问题或建议，请提交Issue或Pull Request。
