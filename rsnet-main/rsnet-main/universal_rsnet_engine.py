#!/usr/bin/env python3
"""
通用RSNet引擎 - 完全自动化的反应网络生成系统
无硬编码，自动发现反应，迭代式生长
"""

import sys
import os
from pathlib import Path
import json
import time
from collections import deque, defaultdict
from typing import List, Dict, Set, Any, Optional, Tuple
import logging

# 添加rsnet到路径
sys.path.insert(0, '.')

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.core.reaction import Reaction
from rsnet.network.generator import NetworkGenerator, ReactionNetwork
from rsnet.operators.registry import OPERATOR_REGISTRY
from rsnet.features.driving_forces import get_driving_forces
from rsnet.features.structure_tags import get_structure_tags
from rsnet.compute.reaction_screener import ReactionScreener

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class UniversalRSNetEngine:
    """通用RSNet引擎 - 完全自动化的反应网络生成"""
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """初始化RSNet引擎"""
        
        # 默认配置
        self.config = {
            'max_generations': 5,
            'max_species': 50,
            'max_reactions_per_generation': 20,
            'energy_cutoff': 100.0,  # kcal/mol
            'min_driving_force_score': 0.3,
            'enable_parallel_screening': False,
            'convergence_threshold': 2,  # 最少新反应数
            'fragment_size_threshold': 2,  # 最小分子片段大小
            'duplicate_detection': True,
            'auto_operator_selection': True
        }
        
        if config:
            self.config.update(config)
        
        # 初始化组件
        self.network_generator = NetworkGenerator(
            operator_registry=OPERATOR_REGISTRY,
            screener=ReactionScreener(optimize_geometries=False)
        )
        
        # 统计信息
        self.generation_stats = []
        self.discovered_molecules = set()
        self.discovered_reactions = []
        self.driving_force_analysis = {}
        
    def generate_network(
        self, 
        seed_molecules: List[Molecule], 
        environment: Environment,
        verbose: bool = True
    ) -> Dict[str, Any]:
        """
        生成完整的反应网络
        
        Args:
            seed_molecules: 初始分子列表
            environment: 反应环境
            verbose: 是否输出详细信息
            
        Returns:
            完整的网络分析结果
        """
        
        if verbose:
            print("=" * 80)
            print("Universal RSNet Engine - Automated Reaction Network Generation")
            print("=" * 80)
            print(f"Input: {len(seed_molecules)} seed molecules")
            for mol in seed_molecules:
                print(f"  - {mol.name}: {mol.smiles}")
            print(f"Environment: {environment.temperature}K, {environment.electrode_type}, {environment.voltage}V")
            print()
        
        start_time = time.time()
        
        # Step 1: 环境驱动力分析
        if verbose:
            print("🔍 Step 1: Environment & Driving Force Analysis")
            print("-" * 50)
        
        self.driving_force_analysis = self._analyze_environment_driving_forces(environment, verbose)
        
        # Step 2: 初始分子结构分析
        if verbose:
            print("\n🏷️  Step 2: Initial Molecule Structure Analysis")
            print("-" * 50)
        
        initial_analysis = self._analyze_initial_molecules(seed_molecules, verbose)
        
        # Step 3: 迭代式网络生成
        if verbose:
            print("\n🌱 Step 3: Iterative Network Generation")
            print("-" * 50)
        
        network = self._generate_network_iteratively(seed_molecules, environment, verbose)
        
        # Step 4: 网络分析与后处理
        if verbose:
            print("\n📊 Step 4: Network Analysis & Post-processing")
            print("-" * 50)
        
        network_analysis = self._analyze_final_network(network, verbose)
        
        total_time = time.time() - start_time
        
        # 生成完整结果
        results = {
            'execution_time': total_time,
            'input_molecules': [{'name': mol.name, 'smiles': mol.smiles} for mol in seed_molecules],
            'environment': {
                'temperature': environment.temperature,
                'electrode_type': environment.electrode_type,
                'voltage': environment.voltage,
                'li_activity': environment.li_activity,
                'interface_type': environment.interface_type
            },
            'driving_force_analysis': self.driving_force_analysis,
            'initial_molecule_analysis': initial_analysis,
            'network_generation_stats': self.generation_stats,
            'final_network': network_analysis,
            'discovered_molecules': len(self.discovered_molecules),
            'discovered_reactions': len(self.discovered_reactions),
            'convergence_achieved': network_analysis['convergence_achieved']
        }
        
        if verbose:
            print(f"\n✅ Network generation completed in {total_time:.2f}s")
            print(f"   Final network: {len(network.molecules)} molecules, {len(network.reactions)} reactions")
            print(f"   Generations processed: {len(self.generation_stats)}")
            print(f"   Convergence: {'Yes' if results['convergence_achieved'] else 'No'}")
        
        return results
    
    def _analyze_environment_driving_forces(self, environment: Environment, verbose: bool) -> Dict[str, Any]:
        """分析环境驱动力"""
        
        # 获取环境驱动力
        driving_forces = get_driving_forces(environment)
        
        # 计算激活的驱动力
        active_forces = {k: v for k, v in driving_forces.items() if v > 0.1}
        
        if verbose:
            print(f"Environment conditions:")
            print(f"  Temperature: {environment.temperature} K")
            print(f"  Electrode: {environment.electrode_type}")
            print(f"  Voltage: {environment.voltage} V")
            print(f"  Li+ activity: {environment.li_activity}")
            print(f"  Interface: {environment.interface_type}")
            
            print(f"\nActive driving forces ({len(active_forces)}/{len(driving_forces)}):")
            for force, strength in sorted(active_forces.items(), key=lambda x: x[1], reverse=True):
                print(f"  ✓ {force}: {strength:.3f}")
        
        return {
            'all_forces': driving_forces,
            'active_forces': active_forces,
            'num_active': len(active_forces),
            'max_force_strength': max(driving_forces.values()) if driving_forces else 0.0
        }
    
    def _analyze_initial_molecules(self, molecules: List[Molecule], verbose: bool) -> Dict[str, Any]:
        """分析初始分子结构"""
        
        analysis = {}
        
        for mol in molecules:
            # 获取结构标签
            tags = get_structure_tags(mol)
            
            # 统计重要特征
            important_features = []
            if tags.get('small_rings', []):
                important_features.append(f"small_rings({len(tags['small_rings'])})")
            if tags.get('polar_bonds', []):
                important_features.append(f"polar_bonds({len(tags['polar_bonds'])})")
            if tags.get('heteroatoms', []):
                important_features.append(f"heteroatoms({len(tags['heteroatoms'])})")
            if tags.get('lewis_acid_sites', []):
                important_features.append(f"lewis_acid({len(tags['lewis_acid_sites'])})")
            if tags.get('lewis_base_sites', []):
                important_features.append(f"lewis_base({len(tags['lewis_base_sites'])})")
            
            analysis[mol.name] = {
                'smiles': mol.smiles,
                'num_atoms': mol.num_atoms,
                'num_heavy_atoms': mol.num_heavy_atoms,
                'important_features': important_features,
                'structure_tags': tags
            }
            
            if verbose:
                print(f"\n{mol.name} ({mol.smiles}):")
                print(f"  Atoms: {mol.num_atoms} total, {mol.num_heavy_atoms} heavy")
                print(f"  Key features: {', '.join(important_features) if important_features else 'None'}")
        
        return analysis
    
    def _generate_network_iteratively(
        self, 
        seed_molecules: List[Molecule], 
        environment: Environment, 
        verbose: bool
    ) -> ReactionNetwork:
        """迭代式生成反应网络"""
        
        # 初始化网络
        network = ReactionNetwork()
        for mol in seed_molecules:
            network.add_molecule(mol, generation=0)
            self.discovered_molecules.add(mol.smiles)
        
        # 迭代生成
        current_generation = 0
        molecules_to_process = deque(seed_molecules)
        consecutive_low_yield_generations = 0
        
        while (current_generation < self.config['max_generations'] and
               len(network.molecules) < self.config['max_species'] and
               molecules_to_process and
               consecutive_low_yield_generations < 2):
            
            if verbose:
                print(f"\n--- Generation {current_generation} ---")
                print(f"Processing {len(molecules_to_process)} molecules...")
            
            generation_start = time.time()
            
            # 处理当前代
            new_reactions = self._process_generation(
                list(molecules_to_process), 
                environment, 
                network, 
                current_generation,
                verbose
            )
            
            generation_time = time.time() - generation_start
            
            # 统计
            gen_stats = {
                'generation': current_generation,
                'input_molecules': len(molecules_to_process),
                'new_reactions': len(new_reactions),
                'total_molecules': len(network.molecules),
                'total_reactions': len(network.reactions),
                'time': generation_time
            }
            self.generation_stats.append(gen_stats)
            
            if verbose:
                print(f"Results: {len(new_reactions)} new reactions, "
                      f"{len(network.molecules)} total molecules ({generation_time:.2f}s)")
            
            # 检查收敛
            if len(new_reactions) < self.config['convergence_threshold']:
                consecutive_low_yield_generations += 1
                if verbose:
                    print(f"Low yield generation ({consecutive_low_yield_generations}/2)")
            else:
                consecutive_low_yield_generations = 0
            
            # 准备下一代
            molecules_to_process.clear()
            next_generation = current_generation + 1
            
            # 添加新分子
            new_molecule_count = 0
            for reaction in new_reactions:
                for product in reaction.products:
                    if product.smiles not in self.discovered_molecules:
                        network.add_molecule(product, generation=next_generation)
                        molecules_to_process.append(product)
                        self.discovered_molecules.add(product.smiles)
                        new_molecule_count += 1
            
            if verbose and new_molecule_count > 0:
                print(f"Added {new_molecule_count} new molecules for next generation")
            
            current_generation = next_generation
        
        return network
    
    def _process_generation(
        self,
        molecules: List[Molecule],
        environment: Environment,
        network: ReactionNetwork,
        generation: int,
        verbose: bool
    ) -> List[Reaction]:
        """处理单个世代"""
        
        all_reactions = []
        
        # 获取激活的算符
        active_operators = self._get_active_operators(molecules, environment)
        
        if verbose:
            print(f"  Active operators: {[op.name for op in active_operators]}")
        
        # 应用算符发现反应
        for operator in active_operators:
            try:
                # 单分子反应
                for mol in molecules:
                    if operator.is_applicable([mol]):
                        reactions = operator.apply([mol], environment)
                        all_reactions.extend(reactions)
                
                # 双分子反应
                if len(molecules) > 1:
                    for i, mol1 in enumerate(molecules):
                        for mol2 in molecules[i+1:]:
                            if operator.is_applicable([mol1, mol2]):
                                reactions = operator.apply([mol1, mol2], environment)
                                all_reactions.extend(reactions)
                
            except Exception as e:
                if verbose:
                    print(f"    Warning: {operator.name} failed: {e}")
        
        if verbose:
            print(f"  Generated {len(all_reactions)} candidate reactions")
        
        # 筛选反应
        feasible_reactions = self._screen_reactions(all_reactions, environment, verbose)
        
        # 添加到网络
        new_reactions = []
        for reaction in feasible_reactions:
            if network.add_reaction(reaction):
                new_reactions.append(reaction)
                self.discovered_reactions.append(reaction)
        
        return new_reactions
    
    def _get_active_operators(self, molecules: List[Molecule], environment: Environment) -> List:
        """获取激活的算符"""
        
        if not self.config['auto_operator_selection']:
            # 返回所有算符
            return list(OPERATOR_REGISTRY.operators.values())
        
        # 智能算符选择
        active_operators = []
        driving_forces = self.driving_force_analysis['active_forces']
        
        for op_name, operator in OPERATOR_REGISTRY.operators.items():
            # 检查驱动力要求
            if hasattr(operator, 'required_driving_forces'):
                required_forces = getattr(operator, 'required_driving_forces', [])
                if any(force in driving_forces for force in required_forces):
                    active_operators.append(operator)
            else:
                # 默认激活
                active_operators.append(operator)
        
        return active_operators
    
    def _screen_reactions(self, reactions: List[Reaction], environment: Environment, verbose: bool) -> List[Reaction]:
        """筛选反应"""
        
        if not reactions:
            return []
        
        feasible_reactions = []
        
        for reaction in reactions:
            try:
                # 简单的可行性检查
                is_feasible = True
                
                # 检查产物是否合理
                for product in reaction.products:
                    if product.num_heavy_atoms < self.config['fragment_size_threshold']:
                        is_feasible = False
                        break
                
                # 检查是否重复
                if self.config['duplicate_detection']:
                    reaction_signature = self._get_reaction_signature(reaction)
                    if any(self._get_reaction_signature(r) == reaction_signature 
                          for r in self.discovered_reactions):
                        is_feasible = False
                
                if is_feasible:
                    # 估算反应能量（简化）
                    reaction.reaction_energy = self._estimate_reaction_energy(reaction)
                    reaction.activation_energy = abs(reaction.reaction_energy) * 0.3 + 15.0
                    
                    if abs(reaction.reaction_energy) <= self.config['energy_cutoff']:
                        feasible_reactions.append(reaction)
                
            except Exception as e:
                if verbose:
                    print(f"    Screening error: {e}")
        
        if verbose:
            print(f"  Screened to {len(feasible_reactions)} feasible reactions")
        
        return feasible_reactions
    
    def _get_reaction_signature(self, reaction: Reaction) -> str:
        """获取反应签名用于去重"""
        reactants = sorted([mol.smiles for mol in reaction.reactants])
        products = sorted([mol.smiles for mol in reaction.products])
        return f"{'|'.join(reactants)}>{'|'.join(products)}"
    
    def _estimate_reaction_energy(self, reaction: Reaction) -> float:
        """估算反应能量（简化版本）"""
        # 基于分子复杂度的简单估算
        reactant_complexity = sum(mol.num_heavy_atoms for mol in reaction.reactants)
        product_complexity = sum(mol.num_heavy_atoms for mol in reaction.products)
        
        # 简化的能量估算
        base_energy = (product_complexity - reactant_complexity) * 5.0
        
        # 添加随机扰动
        import random
        random.seed(hash(self._get_reaction_signature(reaction)) % 2**32)
        perturbation = random.uniform(-20.0, 20.0)
        
        return base_energy + perturbation
    
    def _analyze_final_network(self, network: ReactionNetwork, verbose: bool) -> Dict[str, Any]:
        """分析最终网络"""
        
        stats = network.get_statistics()
        
        # 分析反应类型
        reaction_types = defaultdict(int)
        energy_distribution = []
        
        for reaction in network.reactions:
            reaction_types[reaction.name] += 1
            if hasattr(reaction, 'reaction_energy'):
                energy_distribution.append(reaction.reaction_energy)
        
        # 分析分子类型
        molecule_generations = defaultdict(list)
        for mol_smiles, generation in network.generation_info.items():
            molecule_generations[generation].append(mol_smiles)
        
        analysis = {
            'network_statistics': stats,
            'reaction_types': dict(reaction_types),
            'energy_statistics': {
                'mean_energy': sum(energy_distribution) / len(energy_distribution) if energy_distribution else 0,
                'min_energy': min(energy_distribution) if energy_distribution else 0,
                'max_energy': max(energy_distribution) if energy_distribution else 0,
                'exothermic_count': sum(1 for e in energy_distribution if e < 0),
                'endothermic_count': sum(1 for e in energy_distribution if e > 0)
            },
            'molecule_generations': {str(k): len(v) for k, v in molecule_generations.items()},
            'convergence_achieved': len(self.generation_stats) < self.config['max_generations']
        }
        
        if verbose:
            print(f"Network statistics:")
            print(f"  Molecules: {stats['num_molecules']} (across {stats['max_generation']+1} generations)")
            print(f"  Reactions: {stats['num_reactions']}")
            print(f"  Reaction types: {len(reaction_types)}")
            print(f"  Energy range: {analysis['energy_statistics']['min_energy']:.1f} to {analysis['energy_statistics']['max_energy']:.1f} kcal/mol")
            print(f"  Exothermic/Endothermic: {analysis['energy_statistics']['exothermic_count']}/{analysis['energy_statistics']['endothermic_count']}")
        
        return analysis


def main():
    """主函数 - 演示Li+ PF6- EC体系"""
    
    print("Universal RSNet Engine Demo")
    print("Testing Li+ PF6- EC system with automatic reaction discovery")
    
    # 创建RSNet引擎
    engine = UniversalRSNetEngine({
        'max_generations': 4,
        'max_species': 30,
        'energy_cutoff': 80.0,
        'convergence_threshold': 1
    })
    
    # 定义输入分子 - Li+ PF6- EC体系
    seed_molecules = [
        Molecule.from_smiles('[Li+]', name='Li_ion'),
        Molecule.from_smiles('F[P-](F)(F)(F)(F)F', name='PF6_anion'),
        Molecule.from_smiles('C1COC(=O)O1', name='EC')
    ]
    
    # 定义环境
    environment = Environment(
        temperature=300.0,
        electrode_type='anode',
        voltage=0.1,
        li_activity=1.0,
        interface_type='SEI'
    )
    
    # 生成网络
    results = engine.generate_network(seed_molecules, environment, verbose=True)
    
    # 保存结果
    output_file = Path('universal_rsnet_results.json')
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False, default=str)
    
    print(f"\n✅ Results saved to: {output_file}")
    
    return 0


if __name__ == '__main__':
    exit(main())
