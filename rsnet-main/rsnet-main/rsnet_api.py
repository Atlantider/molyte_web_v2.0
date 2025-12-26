#!/usr/bin/env python3
"""
RSNet API - 通用反应网络生成API
可以作为Web服务或命令行工具使用
"""

import sys
import json
import argparse
from pathlib import Path
from typing import Dict, List, Any, Optional
import time

# 导入简化的RSNet应用程序
from simple_rsnet_app import SimpleRSNetApp


class RSNetAPI:
    """RSNet API接口"""
    
    def __init__(self):
        self.app = SimpleRSNetApp()
    
    def generate_network(self, request: Dict[str, Any]) -> Dict[str, Any]:
        """
        生成反应网络的主要API接口
        
        Args:
            request: 请求字典
                {
                    "molecules": [
                        {"smiles": "C1COC(=O)O1", "name": "EC"},
                        {"smiles": "[Li+]", "name": "Li_ion"},
                        ...
                    ],
                    "environment": {
                        "temperature": 300.0,
                        "electrode_type": "anode",  # "anode" or "cathode"
                        "voltage": 0.1,
                        "li_activity": 1.0,
                        "interface_type": "SEI"
                    },
                    "options": {
                        "max_generations": 4,
                        "max_molecules": 50,
                        "energy_cutoff": 100.0
                    }
                }
        
        Returns:
            response: 响应字典
                {
                    "success": true,
                    "execution_time": 0.123,
                    "network": {
                        "total_molecules": 15,
                        "total_reactions": 8,
                        "generations": 3,
                        "molecules_by_generation": {...},
                        "reactions": [...],
                        "energy_statistics": {...}
                    },
                    "message": "Network generation completed successfully"
                }
        """
        
        try:
            # 验证输入
            validation_result = self._validate_request(request)
            if not validation_result['valid']:
                return {
                    'success': False,
                    'error': validation_result['error'],
                    'message': 'Invalid input parameters'
                }
            
            # 重置应用程序状态
            self.app.molecules = {}
            self.app.reactions = []
            self.app.generation_stats = []
            
            # 运行网络生成
            results = self.app.run(request)
            
            # 格式化响应
            response = {
                'success': True,
                'execution_time': results['execution_time'],
                'network': {
                    'total_molecules': results['total_molecules'],
                    'total_reactions': results['total_reactions'],
                    'generations': results['generations'],
                    'generation_stats': results['generation_stats'],
                    'molecules_by_generation': results['molecules_by_generation'],
                    'reactions': results['reactions'],
                    'reaction_types': results['reaction_types'],
                    'energy_statistics': results['energy_statistics']
                },
                'message': 'Network generation completed successfully'
            }
            
            return response
            
        except Exception as e:
            return {
                'success': False,
                'error': str(e),
                'message': 'Network generation failed'
            }
    
    def _validate_request(self, request: Dict[str, Any]) -> Dict[str, Any]:
        """验证请求参数"""
        
        # 检查必需字段
        if 'molecules' not in request:
            return {'valid': False, 'error': 'Missing molecules field'}
        
        if 'environment' not in request:
            return {'valid': False, 'error': 'Missing environment field'}
        
        # 验证分子
        molecules = request['molecules']
        if not isinstance(molecules, list) or len(molecules) == 0:
            return {'valid': False, 'error': 'Molecules must be a non-empty list'}
        
        for i, mol in enumerate(molecules):
            if not isinstance(mol, dict):
                return {'valid': False, 'error': f'Molecule {i} must be a dictionary'}
            
            if 'smiles' not in mol or 'name' not in mol:
                return {'valid': False, 'error': f'Molecule {i} must have smiles and name fields'}
            
            if not isinstance(mol['smiles'], str) or not isinstance(mol['name'], str):
                return {'valid': False, 'error': f'Molecule {i} smiles and name must be strings'}
        
        # 验证环境
        environment = request['environment']
        if not isinstance(environment, dict):
            return {'valid': False, 'error': 'Environment must be a dictionary'}
        
        required_env_fields = ['temperature', 'electrode_type', 'voltage']
        for field in required_env_fields:
            if field not in environment:
                return {'valid': False, 'error': f'Missing environment field: {field}'}
        
        # 验证环境参数范围
        if not (200 <= environment['temperature'] <= 500):
            return {'valid': False, 'error': 'Temperature must be between 200-500 K'}
        
        if environment['electrode_type'] not in ['anode', 'cathode']:
            return {'valid': False, 'error': 'Electrode type must be anode or cathode'}
        
        if not (-2.0 <= environment['voltage'] <= 5.0):
            return {'valid': False, 'error': 'Voltage must be between -2.0 and 5.0 V'}
        
        return {'valid': True}
    
    def get_supported_molecules(self) -> Dict[str, Any]:
        """获取支持的分子类型"""
        return {
            'common_solvents': [
                {'name': 'EC', 'smiles': 'C1COC(=O)O1', 'description': 'Ethylene Carbonate'},
                {'name': 'DMC', 'smiles': 'COC(=O)OC', 'description': 'Dimethyl Carbonate'},
                {'name': 'PC', 'smiles': 'CC1COC(=O)O1', 'description': 'Propylene Carbonate'},
                {'name': 'EMC', 'smiles': 'CCOC(=O)OC', 'description': 'Ethyl Methyl Carbonate'}
            ],
            'common_salts': [
                {'name': 'LiPF6', 'cation': '[Li+]', 'anion': 'F[P-](F)(F)(F)(F)F', 'description': 'Lithium Hexafluorophosphate'},
                {'name': 'LiTFSI', 'cation': '[Li+]', 'anion': 'O=S(=O)(N([S](=O)(=O)C(F)(F)F))[C](F)(F)F', 'description': 'Lithium Bis(trifluoromethanesulfonyl)imide'},
                {'name': 'NaClO4', 'cation': '[Na+]', 'anion': '[O-][Cl+3]([O-])([O-])[O-]', 'description': 'Sodium Perchlorate'}
            ],
            'common_additives': [
                {'name': 'VC', 'smiles': 'C1COC(=O)O1', 'description': 'Vinylene Carbonate'},
                {'name': 'FEC', 'smiles': 'FC1COC(=O)O1', 'description': 'Fluoroethylene Carbonate'}
            ]
        }
    
    def get_environment_presets(self) -> Dict[str, Any]:
        """获取环境预设"""
        return {
            'li_ion_anode_sei': {
                'temperature': 300.0,
                'electrode_type': 'anode',
                'voltage': 0.1,
                'li_activity': 1.0,
                'interface_type': 'SEI',
                'description': 'Li-ion battery anode SEI formation conditions'
            },
            'li_ion_cathode': {
                'temperature': 300.0,
                'electrode_type': 'cathode',
                'voltage': 4.2,
                'li_activity': 0.1,
                'interface_type': 'CEI',
                'description': 'Li-ion battery cathode conditions'
            },
            'na_ion_anode': {
                'temperature': 300.0,
                'electrode_type': 'anode',
                'voltage': 0.2,
                'li_activity': 0.0,
                'interface_type': 'SEI',
                'description': 'Na-ion battery anode conditions'
            },
            'high_temperature': {
                'temperature': 400.0,
                'electrode_type': 'anode',
                'voltage': 0.1,
                'li_activity': 1.0,
                'interface_type': 'SEI',
                'description': 'High temperature conditions'
            }
        }


def create_example_requests() -> List[Dict[str, Any]]:
    """创建示例请求"""
    
    examples = [
        {
            'name': 'Li+ PF6- EC System',
            'request': {
                'molecules': [
                    {'smiles': '[Li+]', 'name': 'Li_ion'},
                    {'smiles': 'F[P-](F)(F)(F)(F)F', 'name': 'PF6_anion'},
                    {'smiles': 'C1COC(=O)O1', 'name': 'EC'}
                ],
                'environment': {
                    'temperature': 300.0,
                    'electrode_type': 'anode',
                    'voltage': 0.1,
                    'li_activity': 1.0,
                    'interface_type': 'SEI'
                }
            }
        },
        {
            'name': 'Na+ ClO4- DMC System',
            'request': {
                'molecules': [
                    {'smiles': '[Na+]', 'name': 'Na_ion'},
                    {'smiles': '[O-][Cl+3]([O-])([O-])[O-]', 'name': 'ClO4_anion'},
                    {'smiles': 'COC(=O)OC', 'name': 'DMC'}
                ],
                'environment': {
                    'temperature': 300.0,
                    'electrode_type': 'anode',
                    'voltage': 0.2,
                    'li_activity': 0.0,
                    'interface_type': 'SEI'
                }
            }
        },
        {
            'name': 'Li+ TFSI- PC System',
            'request': {
                'molecules': [
                    {'smiles': '[Li+]', 'name': 'Li_ion'},
                    {'smiles': 'O=S(=O)(N([S](=O)(=O)C(F)(F)F))[C](F)(F)F', 'name': 'TFSI_anion'},
                    {'smiles': 'CC1COC(=O)O1', 'name': 'PC'}
                ],
                'environment': {
                    'temperature': 350.0,
                    'electrode_type': 'anode',
                    'voltage': 0.15,
                    'li_activity': 1.0,
                    'interface_type': 'SEI'
                }
            }
        }
    ]
    
    return examples


def main():
    """命令行接口"""
    
    parser = argparse.ArgumentParser(description='RSNet API - Universal Reaction Network Generator')
    parser.add_argument('--input', '-i', type=str, help='Input JSON file with request')
    parser.add_argument('--output', '-o', type=str, help='Output JSON file for results')
    parser.add_argument('--example', '-e', type=str, help='Run example (li_pf6_ec, na_clo4_dmc, li_tfsi_pc)')
    parser.add_argument('--list-molecules', action='store_true', help='List supported molecules')
    parser.add_argument('--list-environments', action='store_true', help='List environment presets')
    parser.add_argument('--create-examples', action='store_true', help='Create example input files')
    
    args = parser.parse_args()
    
    api = RSNetAPI()
    
    # 列出支持的分子
    if args.list_molecules:
        molecules = api.get_supported_molecules()
        print("Supported Molecules:")
        print(json.dumps(molecules, indent=2))
        return 0
    
    # 列出环境预设
    if args.list_environments:
        environments = api.get_environment_presets()
        print("Environment Presets:")
        print(json.dumps(environments, indent=2))
        return 0
    
    # 创建示例文件
    if args.create_examples:
        examples = create_example_requests()
        for example in examples:
            filename = example['name'].lower().replace(' ', '_').replace('+', '').replace('-', '') + '_request.json'
            with open(filename, 'w') as f:
                json.dump(example['request'], f, indent=2)
            print(f"Created example: {filename}")
        return 0
    
    # 运行示例
    if args.example:
        examples = create_example_requests()
        example_map = {
            'li_pf6_ec': 0,
            'na_clo4_dmc': 1,
            'li_tfsi_pc': 2
        }
        
        if args.example not in example_map:
            print(f"Unknown example: {args.example}")
            print(f"Available examples: {list(example_map.keys())}")
            return 1
        
        example = examples[example_map[args.example]]
        print(f"Running example: {example['name']}")
        
        response = api.generate_network(example['request'])
        
        if response['success']:
            print(f"✅ Success! Generated network with {response['network']['total_molecules']} molecules and {response['network']['total_reactions']} reactions")
            
            output_file = args.output or f"{args.example}_results.json"
            with open(output_file, 'w') as f:
                json.dump(response, f, indent=2, default=str)
            print(f"Results saved to: {output_file}")
        else:
            print(f"❌ Failed: {response['message']}")
            print(f"Error: {response['error']}")
        
        return 0
    
    # 从文件运行
    if args.input:
        if not Path(args.input).exists():
            print(f"Input file not found: {args.input}")
            return 1
        
        with open(args.input, 'r') as f:
            request = json.load(f)
        
        print(f"Processing request from: {args.input}")
        response = api.generate_network(request)
        
        if response['success']:
            print(f"✅ Success! Generated network with {response['network']['total_molecules']} molecules and {response['network']['total_reactions']} reactions")
            
            output_file = args.output or args.input.replace('.json', '_results.json')
            with open(output_file, 'w') as f:
                json.dump(response, f, indent=2, default=str)
            print(f"Results saved to: {output_file}")
        else:
            print(f"❌ Failed: {response['message']}")
            print(f"Error: {response['error']}")
        
        return 0
    
    # 默认：显示帮助
    parser.print_help()
    return 0


if __name__ == '__main__':
    exit(main())
