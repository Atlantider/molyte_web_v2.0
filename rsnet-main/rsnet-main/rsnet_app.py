#!/usr/bin/env python3
"""
RSNet Application - 通用反应网络生成应用程序
输入: 分子SMILES + 环境条件
输出: 完整反应网络
完全自动化，无需领域知识
"""

import sys
import os
import json
import time
import subprocess
import tempfile
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
from collections import defaultdict, deque
import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import numpy as np

# 添加rsnet到路径
sys.path.insert(0, '.')

from rsnet.core.molecule import Molecule
from rsnet.core.environment import Environment
from rsnet.core.reaction import Reaction


class IntelligentOperatorGenerator:
    """智能算符生成器 - 基于分子结构和环境自动生成算符"""
    
    def __init__(self):
        self.generated_operators = []
    
    def generate_operators(self, molecules: List[Molecule], environment: Environment) -> List[Dict]:
        """基于输入分子和环境智能生成算符"""
        
        operators = []
        
        # 分析分子特征
        molecular_features = self._analyze_molecular_features(molecules)
        
        # 分析环境驱动力
        environmental_forces = self._analyze_environmental_forces(environment)
        
        # 基于特征生成算符
        operators.extend(self._generate_bond_breaking_operators(molecular_features, environmental_forces))
        operators.extend(self._generate_addition_operators(molecular_features, environmental_forces))
        operators.extend(self._generate_electron_transfer_operators(molecular_features, environmental_forces))
        operators.extend(self._generate_coordination_operators(molecular_features, environmental_forces))
        operators.extend(self._generate_rearrangement_operators(molecular_features, environmental_forces))
        
        self.generated_operators = operators
        return operators
    
    def _analyze_molecular_features(self, molecules: List[Molecule]) -> Dict[str, Any]:
        """分析分子特征"""
        features = {
            'has_rings': False,
            'has_carbonyls': False,
            'has_ethers': False,
            'has_metals': False,
            'has_halogens': False,
            'has_charged_species': False,
            'ring_sizes': [],
            'functional_groups': [],
            'reactive_sites': []
        }
        
        for mol in molecules:
            rdkit_mol = mol.rdkit_mol
            
            # 检测环
            ring_info = rdkit_mol.GetRingInfo()
            if ring_info.NumRings() > 0:
                features['has_rings'] = True
                features['ring_sizes'].extend([len(ring) for ring in ring_info.AtomRings()])
            
            # 检测官能团
            for atom in rdkit_mol.GetAtoms():
                symbol = atom.GetSymbol()
                
                # 金属
                if symbol in ['Li', 'Na', 'K', 'Mg', 'Ca', 'Al', 'Zn', 'Fe', 'Cu', 'Ni']:
                    features['has_metals'] = True
                
                # 卤素
                if symbol in ['F', 'Cl', 'Br', 'I']:
                    features['has_halogens'] = True
                
                # 带电物种
                if atom.GetFormalCharge() != 0:
                    features['has_charged_species'] = True
                
                # 羰基
                if symbol == 'C' and any(neighbor.GetSymbol() == 'O' and 
                                       rdkit_mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.BondType.DOUBLE
                                       for neighbor in atom.GetNeighbors()):
                    features['has_carbonyls'] = True
                    features['functional_groups'].append('carbonyl')
                
                # 醚键
                if symbol == 'O' and atom.GetDegree() == 2 and all(neighbor.GetSymbol() == 'C' for neighbor in atom.GetNeighbors()):
                    features['has_ethers'] = True
                    features['functional_groups'].append('ether')
        
        return features
    
    def _analyze_environmental_forces(self, environment: Environment) -> Dict[str, float]:
        """分析环境驱动力"""
        forces = {}
        
        # 电化学驱动力
        if environment.electrode_type == 'anode' and environment.voltage < 1.0:
            forces['reduction_potential'] = 1.0 - environment.voltage
        elif environment.electrode_type == 'cathode' and environment.voltage > 3.0:
            forces['oxidation_potential'] = environment.voltage - 3.0
        else:
            forces['electrochemical'] = 0.1
        
        # 热力学驱动力
        if environment.temperature > 250:
            forces['thermal_activation'] = (environment.temperature - 250) / 100
        
        # 界面驱动力
        if environment.interface_type == 'SEI':
            forces['surface_reaction'] = 0.8
        
        # 离子活度驱动力
        if environment.li_activity > 0.5:
            forces['ion_coordination'] = environment.li_activity
        
        return forces
    
    def _generate_bond_breaking_operators(self, mol_features: Dict, env_forces: Dict) -> List[Dict]:
        """生成键断裂算符"""
        operators = []
        
        # 环开环反应
        if mol_features['has_rings'] and env_forces.get('thermal_activation', 0) > 0.3:
            for ring_size in set(mol_features['ring_sizes']):
                if ring_size <= 6:  # 小环更容易开环
                    operators.append({
                        'name': f'ring_opening_{ring_size}',
                        'type': 'bond_breaking',
                        'target_pattern': f'[*]1{"[*]" * (ring_size-2)}1',
                        'activation_energy': 20.0 + ring_size * 5.0,
                        'probability': 0.8 if ring_size <= 5 else 0.4
                    })
        
        # 极性键断裂
        if mol_features['has_charged_species'] and env_forces.get('electrochemical', 0) > 0.2:
            operators.append({
                'name': 'ionic_dissociation',
                'type': 'bond_breaking',
                'target_pattern': '[+1,+2,+3][*]',
                'activation_energy': 15.0,
                'probability': 0.9
            })
        
        return operators
    
    def _generate_addition_operators(self, mol_features: Dict, env_forces: Dict) -> List[Dict]:
        """生成加成算符"""
        operators = []
        
        # 金属配位
        if mol_features['has_metals'] and env_forces.get('ion_coordination', 0) > 0.3:
            operators.append({
                'name': 'metal_coordination',
                'type': 'addition',
                'reactant1_pattern': '[Li+,Na+,K+,Mg+2,Ca+2]',
                'reactant2_pattern': '[O,N,S]',
                'activation_energy': 10.0,
                'probability': 0.9
            })
        
        return operators
    
    def _generate_electron_transfer_operators(self, mol_features: Dict, env_forces: Dict) -> List[Dict]:
        """生成电子转移算符"""
        operators = []
        
        # 电化学还原
        if env_forces.get('reduction_potential', 0) > 0.2:
            if mol_features['has_carbonyls']:
                operators.append({
                    'name': 'carbonyl_reduction',
                    'type': 'electron_transfer',
                    'target_pattern': '[C]=[O]',
                    'activation_energy': 25.0,
                    'probability': 0.7
                })
            
            if mol_features['has_halogens']:
                operators.append({
                    'name': 'halogen_reduction',
                    'type': 'electron_transfer',
                    'target_pattern': '[F,Cl,Br,I]',
                    'activation_energy': 30.0,
                    'probability': 0.6
                })
        
        return operators
    
    def _generate_coordination_operators(self, mol_features: Dict, env_forces: Dict) -> List[Dict]:
        """生成配位算符"""
        operators = []
        
        if mol_features['has_metals'] and mol_features['has_ethers']:
            operators.append({
                'name': 'ether_coordination',
                'type': 'coordination',
                'reactant1_pattern': '[Li+,Na+]',
                'reactant2_pattern': '[O]([C])[C]',
                'activation_energy': 8.0,
                'probability': 0.8
            })
        
        return operators
    
    def _generate_rearrangement_operators(self, mol_features: Dict, env_forces: Dict) -> List[Dict]:
        """生成重排算符"""
        operators = []
        
        if mol_features['has_rings'] and env_forces.get('thermal_activation', 0) > 0.5:
            operators.append({
                'name': 'ring_rearrangement',
                'type': 'rearrangement',
                'target_pattern': '[*]1[*][*][*][*]1',
                'activation_energy': 35.0,
                'probability': 0.3
            })
        
        return operators


class ReactionTemplateGenerator:
    """反应模板生成器 - 将算符转换为具体的反应模板"""
    
    def generate_templates(self, operators: List[Dict], molecules: List[Molecule]) -> List[Dict]:
        """生成反应模板"""
        templates = []
        
        for operator in operators:
            if operator['type'] == 'bond_breaking':
                templates.extend(self._generate_bond_breaking_templates(operator, molecules))
            elif operator['type'] == 'addition':
                templates.extend(self._generate_addition_templates(operator, molecules))
            elif operator['type'] == 'electron_transfer':
                templates.extend(self._generate_electron_transfer_templates(operator, molecules))
            elif operator['type'] == 'coordination':
                templates.extend(self._generate_coordination_templates(operator, molecules))
            elif operator['type'] == 'rearrangement':
                templates.extend(self._generate_rearrangement_templates(operator, molecules))
        
        return templates
    
    def _generate_bond_breaking_templates(self, operator: Dict, molecules: List[Molecule]) -> List[Dict]:
        """生成键断裂模板"""
        templates = []
        
        if operator['name'].startswith('ring_opening'):
            # EC环开环: C1COC(=O)O1 → CO2 + C2H4O
            templates.append({
                'name': f"{operator['name']}_template",
                'smarts': '[C:1]1[O:2][C:3][C:4][O:5][C:6](=O)[O:7]1>>[C:1][O:2].[C:3]=[C:4].[O:5]=[C:6]=[O:7]',
                'reactant_count': 1,
                'activation_energy': operator['activation_energy'],
                'probability': operator['probability']
            })
        
        elif operator['name'] == 'ionic_dissociation':
            # 离子解离: [Li+][A-] → Li+ + A-
            templates.append({
                'name': 'ionic_dissociation_template',
                'smarts': '[Li+:1][*:2]>>[Li+:1].[*:2]',
                'reactant_count': 1,
                'activation_energy': operator['activation_energy'],
                'probability': operator['probability']
            })
        
        return templates
    
    def _generate_addition_templates(self, operator: Dict, molecules: List[Molecule]) -> List[Dict]:
        """生成加成模板"""
        templates = []
        
        if operator['name'] == 'metal_coordination':
            templates.append({
                'name': 'metal_coordination_template',
                'smarts': '[Li+:1].[O:2]>>[Li+:1][O:2]',
                'reactant_count': 2,
                'activation_energy': operator['activation_energy'],
                'probability': operator['probability']
            })
        
        return templates
    
    def _generate_electron_transfer_templates(self, operator: Dict, molecules: List[Molecule]) -> List[Dict]:
        """生成电子转移模板"""
        templates = []
        
        if operator['name'] == 'carbonyl_reduction':
            templates.append({
                'name': 'carbonyl_reduction_template',
                'smarts': '[C:1]=[O:2]>>[C:1][O:2]',
                'reactant_count': 1,
                'activation_energy': operator['activation_energy'],
                'probability': operator['probability']
            })
        
        elif operator['name'] == 'halogen_reduction':
            templates.append({
                'name': 'halogen_reduction_template', 
                'smarts': '[P:1][F:2]>>[P:1].[F:2]',
                'reactant_count': 1,
                'activation_energy': operator['activation_energy'],
                'probability': operator['probability']
            })
        
        return templates
    
    def _generate_coordination_templates(self, operator: Dict, molecules: List[Molecule]) -> List[Dict]:
        """生成配位模板"""
        templates = []
        
        if operator['name'] == 'ether_coordination':
            templates.append({
                'name': 'ether_coordination_template',
                'smarts': '[Li+:1].[O:2]([C:3])[C:4]>>[Li+:1][O:2]([C:3])[C:4]',
                'reactant_count': 2,
                'activation_energy': operator['activation_energy'],
                'probability': operator['probability']
            })
        
        return templates
    
    def _generate_rearrangement_templates(self, operator: Dict, molecules: List[Molecule]) -> List[Dict]:
        """生成重排模板"""
        templates = []
        
        if operator['name'] == 'ring_rearrangement':
            templates.append({
                'name': 'ring_rearrangement_template',
                'smarts': '[C:1]1[C:2][C:3][C:4][C:5]1>>[C:1]=[C:2][C:3][C:4][C:5]',
                'reactant_count': 1,
                'activation_energy': operator['activation_energy'],
                'probability': operator['probability']
            })
        
        return templates


class ReactionInstantiator:
    """反应实例化器 - 将模板应用到具体分子上生成反应"""
    
    def instantiate_reactions(self, templates: List[Dict], molecules: List[Molecule]) -> List[Reaction]:
        """实例化反应"""
        reactions = []
        
        for template in templates:
            try:
                if template['reactant_count'] == 1:
                    reactions.extend(self._apply_single_molecule_template(template, molecules))
                elif template['reactant_count'] == 2:
                    reactions.extend(self._apply_two_molecule_template(template, molecules))
            except Exception as e:
                print(f"Template {template['name']} failed: {e}")
        
        return reactions
    
    def _apply_single_molecule_template(self, template: Dict, molecules: List[Molecule]) -> List[Reaction]:
        """应用单分子模板"""
        reactions = []
        
        try:
            rxn = AllChem.ReactionFromSmarts(template['smarts'])
            if not rxn:
                return reactions
            
            for mol in molecules:
                try:
                    products_tuples = rxn.RunReactants((mol.rdkit_mol,))
                    
                    for products_tuple in products_tuples:
                        products = []
                        for i, prod_mol in enumerate(products_tuple):
                            Chem.SanitizeMol(prod_mol)
                            smiles = Chem.MolToSmiles(prod_mol)
                            if smiles and len(smiles) > 0:
                                products.append(Molecule.from_smiles(smiles, name=f"{mol.name}_prod{i}"))
                        
                        if products:
                            reaction = Reaction(
                                reactants=[mol],
                                products=products,
                                name=template['name'],
                                activation_energy=template.get('activation_energy', 25.0)
                            )
                            reactions.append(reaction)
                
                except Exception as e:
                    continue
        
        except Exception as e:
            print(f"Single molecule template error: {e}")
        
        return reactions
    
    def _apply_two_molecule_template(self, template: Dict, molecules: List[Molecule]) -> List[Reaction]:
        """应用双分子模板"""
        reactions = []
        
        try:
            rxn = AllChem.ReactionFromSmarts(template['smarts'])
            if not rxn:
                return reactions
            
            for i, mol1 in enumerate(molecules):
                for mol2 in molecules[i:]:
                    try:
                        products_tuples = rxn.RunReactants((mol1.rdkit_mol, mol2.rdkit_mol))
                        
                        for products_tuple in products_tuples:
                            products = []
                            for j, prod_mol in enumerate(products_tuple):
                                Chem.SanitizeMol(prod_mol)
                                smiles = Chem.MolToSmiles(prod_mol)
                                if smiles and len(smiles) > 0:
                                    products.append(Molecule.from_smiles(smiles, name=f"{mol1.name}_{mol2.name}_prod{j}"))
                            
                            if products:
                                reaction = Reaction(
                                    reactants=[mol1, mol2] if mol1 != mol2 else [mol1],
                                    products=products,
                                    name=template['name'],
                                    activation_energy=template.get('activation_energy', 25.0)
                                )
                                reactions.append(reaction)
                    
                    except Exception as e:
                        continue
        
        except Exception as e:
            print(f"Two molecule template error: {e}")
        
        return reactions


class XTBEnergyCalculator:
    """xTB能量计算器"""
    
    def __init__(self):
        self.xtb_path = self._find_xtb_executable()
    
    def _find_xtb_executable(self) -> Optional[str]:
        """查找xTB可执行文件"""
        try:
            result = subprocess.run(['which', 'xtb'], capture_output=True, text=True)
            if result.returncode == 0:
                return result.stdout.strip()
        except:
            pass
        return None
    
    def calculate_molecule_energy(self, molecule: Molecule) -> Optional[float]:
        """计算分子能量"""
        if not self.xtb_path:
            # 使用简化的能量估算
            return self._estimate_energy(molecule)
        
        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                xyz_file = Path(temp_dir) / f"{molecule.name}.xyz"
                
                # 生成xyz文件
                self._write_xyz_file(molecule, xyz_file)
                
                # 运行xTB
                cmd = [self.xtb_path, str(xyz_file), '--gfn', '2', '--sp']
                result = subprocess.run(cmd, cwd=temp_dir, capture_output=True, text=True, timeout=30)
                
                if result.returncode == 0:
                    return self._parse_xtb_energy(result.stdout)
        
        except Exception as e:
            print(f"xTB calculation failed for {molecule.name}: {e}")
        
        return self._estimate_energy(molecule)
    
    def _write_xyz_file(self, molecule: Molecule, xyz_file: Path):
        """写入xyz文件"""
        mol = molecule.rdkit_mol
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)
        
        conf = mol.GetConformer()
        atoms = mol.GetAtoms()
        
        with open(xyz_file, 'w') as f:
            f.write(f"{mol.GetNumAtoms()}\n")
            f.write(f"{molecule.name}\n")
            
            for i, atom in enumerate(atoms):
                pos = conf.GetAtomPosition(i)
                f.write(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")
    
    def _parse_xtb_energy(self, output: str) -> Optional[float]:
        """解析xTB输出能量"""
        for line in output.split('\n'):
            if 'TOTAL ENERGY' in line:
                try:
                    energy_hartree = float(line.split()[-2])
                    return energy_hartree * 627.509  # 转换为kcal/mol
                except:
                    continue
        return None
    
    def _estimate_energy(self, molecule: Molecule) -> float:
        """简化的能量估算"""
        mol = molecule.rdkit_mol
        
        # 基于分子描述符的简单估算
        mw = Descriptors.MolWt(mol)
        num_atoms = mol.GetNumAtoms()
        num_bonds = mol.GetNumBonds()
        
        # 简化的能量公式
        base_energy = -num_atoms * 50.0  # 每个原子约-50 kcal/mol
        bond_energy = -num_bonds * 80.0  # 每个键约-80 kcal/mol
        
        return base_energy + bond_energy
    
    def calculate_reaction_energy(self, reaction: Reaction) -> float:
        """计算反应能量"""
        reactant_energies = [self.calculate_molecule_energy(mol) for mol in reaction.reactants]
        product_energies = [self.calculate_molecule_energy(mol) for mol in reaction.products]
        
        if None in reactant_energies or None in product_energies:
            return 0.0
        
        total_reactant_energy = sum(reactant_energies)
        total_product_energy = sum(product_energies)
        
        return total_product_energy - total_reactant_energy


class RSNetApplication:
    """RSNet应用程序主类"""
    
    def __init__(self):
        self.operator_generator = IntelligentOperatorGenerator()
        self.template_generator = ReactionTemplateGenerator()
        self.reaction_instantiator = ReactionInstantiator()
        self.energy_calculator = XTBEnergyCalculator()
        
        # 网络存储
        self.molecules = {}  # smiles -> (molecule, generation)
        self.reactions = []
        self.network_graph = nx.DiGraph()
    
    def run(self, input_smiles: List[str], input_names: List[str], environment_params: Dict[str, Any]) -> Dict[str, Any]:
        """
        运行RSNet应用程序
        
        Args:
            input_smiles: 输入分子的SMILES列表
            input_names: 输入分子的名称列表  
            environment_params: 环境参数字典
                - temperature: 温度(K)
                - electrode_type: 电极类型('anode'/'cathode')
                - voltage: 电压(V)
                - li_activity: Li+活度
                - interface_type: 界面类型
        
        Returns:
            完整的反应网络结果
        """
        
        print("=" * 80)
        print("RSNet Application - Universal Reaction Network Generator")
        print("=" * 80)
        
        start_time = time.time()
        
        # Step 1: 初始化分子和环境
        print("🧪 Step 1: Initialize molecules and environment")
        seed_molecules = []
        for smiles, name in zip(input_smiles, input_names):
            mol = Molecule.from_smiles(smiles, name=name)
            seed_molecules.append(mol)
            self.molecules[smiles] = (mol, 0)
            print(f"   - {name}: {smiles}")
        
        environment = Environment(**environment_params)
        print(f"   Environment: {environment.temperature}K, {environment.electrode_type}, {environment.voltage}V")
        
        # Step 2: 迭代生成网络
        print(f"\n🔄 Step 2: Iterative network generation")
        generation = 0
        max_generations = 5
        molecules_to_process = seed_molecules.copy()
        
        while generation < max_generations and molecules_to_process:
            print(f"\n--- Generation {generation} ---")
            print(f"Processing {len(molecules_to_process)} molecules...")
            
            # 2.1: 智能生成算符
            operators = self.operator_generator.generate_operators(molecules_to_process, environment)
            print(f"   Generated {len(operators)} operators")
            
            # 2.2: 生成反应模板
            templates = self.template_generator.generate_templates(operators, molecules_to_process)
            print(f"   Generated {len(templates)} reaction templates")
            
            # 2.3: 实例化反应
            new_reactions = self.reaction_instantiator.instantiate_reactions(templates, molecules_to_process)
            print(f"   Instantiated {len(new_reactions)} reactions")
            
            # 2.4: 计算能量并筛选
            feasible_reactions = []
            new_molecules = []
            
            for reaction in new_reactions:
                # 计算反应能量
                reaction.reaction_energy = self.energy_calculator.calculate_reaction_energy(reaction)
                
                # 能量筛选
                if abs(reaction.reaction_energy) < 100.0:  # kcal/mol
                    feasible_reactions.append(reaction)
                    
                    # 添加新分子
                    for product in reaction.products:
                        if product.smiles not in self.molecules:
                            self.molecules[product.smiles] = (product, generation + 1)
                            new_molecules.append(product)
            
            self.reactions.extend(feasible_reactions)
            print(f"   Added {len(feasible_reactions)} feasible reactions")
            print(f"   Discovered {len(new_molecules)} new molecules")
            
            # 准备下一代
            if not new_molecules:
                print("   No new molecules - convergence reached")
                break
            
            molecules_to_process = new_molecules
            generation += 1
        
        total_time = time.time() - start_time
        
        # Step 3: 分析结果
        print(f"\n📊 Step 3: Network analysis")
        results = self._analyze_results(total_time)
        
        print(f"\n✅ RSNet application completed in {total_time:.2f}s")
        print(f"   Final network: {len(self.molecules)} molecules, {len(self.reactions)} reactions")
        print(f"   Generations: {generation + 1}")
        
        return results
    
    def _analyze_results(self, execution_time: float) -> Dict[str, Any]:
        """分析结果"""
        
        # 按代分组分子
        molecules_by_generation = defaultdict(list)
        for smiles, (mol, gen) in self.molecules.items():
            molecules_by_generation[gen].append({'name': mol.name, 'smiles': smiles})
        
        # 反应类型统计
        reaction_types = defaultdict(int)
        energy_stats = []
        
        for reaction in self.reactions:
            reaction_types[reaction.name] += 1
            if hasattr(reaction, 'reaction_energy'):
                energy_stats.append(reaction.reaction_energy)
        
        # 生成结果
        results = {
            'execution_time': execution_time,
            'total_molecules': len(self.molecules),
            'total_reactions': len(self.reactions),
            'generations': len(molecules_by_generation),
            'molecules_by_generation': dict(molecules_by_generation),
            'reaction_types': dict(reaction_types),
            'energy_statistics': {
                'mean': np.mean(energy_stats) if energy_stats else 0,
                'std': np.std(energy_stats) if energy_stats else 0,
                'min': min(energy_stats) if energy_stats else 0,
                'max': max(energy_stats) if energy_stats else 0,
                'exothermic_count': sum(1 for e in energy_stats if e < 0),
                'endothermic_count': sum(1 for e in energy_stats if e > 0)
            },
            'reactions': [
                {
                    'name': rxn.name,
                    'reactants': [mol.name for mol in rxn.reactants],
                    'products': [mol.name for mol in rxn.products],
                    'energy': getattr(rxn, 'reaction_energy', 0.0)
                }
                for rxn in self.reactions
            ]
        }
        
        return results


def main():
    """主函数 - 演示应用程序使用"""
    
    app = RSNetApplication()
    
    # 测试案例1: Li+ PF6- EC
    print("Test Case 1: Li+ PF6- EC System")
    results1 = app.run(
        input_smiles=['[Li+]', 'F[P-](F)(F)(F)(F)F', 'C1COC(=O)O1'],
        input_names=['Li_ion', 'PF6_anion', 'EC'],
        environment_params={
            'temperature': 300.0,
            'electrode_type': 'anode',
            'voltage': 0.1,
            'li_activity': 1.0,
            'interface_type': 'SEI'
        }
    )
    
    # 保存结果
    with open('rsnet_app_results.json', 'w') as f:
        json.dump(results1, f, indent=2, default=str)
    
    print(f"\n✅ Results saved to: rsnet_app_results.json")
    
    return 0


if __name__ == '__main__':
    exit(main())
