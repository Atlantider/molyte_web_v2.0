"""
原子映射生成器

基于 .lt 文件生成 atom_mapping.json
"""

import json
from pathlib import Path
from typing import Dict, List, Optional
import logging

from .lt_parser import LTParser

logger = logging.getLogger(__name__)


class AtomMappingGenerator:
    """生成原子映射文件"""
    
    def __init__(self, initial_salts_path: Path):
        """
        Args:
            initial_salts_path: 阳离子/阴离子初始结构文件路径
        """
        self.parser = LTParser()
        self.initial_salts_path = initial_salts_path
    
    def generate_atom_mapping(
        self,
        job_data: Dict,
        work_path: Path
    ) -> Dict:
        """
        生成原子映射文件
        
        Args:
            job_data: 任务数据，格式：
                {
                    "name": "MD-20251119-0001-xxx",
                    "cations": [{"name": "Li", "number": 50}],
                    "anions": [{"name": "PF6", "number": 50}],
                    "solvents": [
                        {"name": "TEP", "smiles": "...", "number": 100},
                        {"name": "CO2", "smiles": "...", "number": 10}
                    ]
                }
            work_path: 工作目录路径
        
        Returns:
            {
              "molecules": [
                {
                  "molecule_id": 1,
                  "molecule_name": "Li",
                  "molecule_type": "cation",
                  "atoms": [
                    {
                      "atom_id": 1,
                      "atom_index": 1,
                      "atom_name": "Li",
                      "label": "Li_Li",
                      "element": "Li",
                      "type": "simple_ion"
                    }
                  ]
                },
                ...
              ]
            }
        """
        mapping = {"molecules": []}
        
        atom_id_counter = 1
        molecule_id_counter = 1
        
        # 处理阳离子
        for cation in job_data.get("cations", []):
            name = cation["name"]
            number = cation["number"]
            
            # 获取分子信息
            mol_info = self._get_molecule_info(name, "cation", work_path)
            
            # 为每个分子实例创建映射
            for i in range(number):
                mol_entry = self._create_molecule_entry(
                    molecule_id_counter,
                    name,
                    "cation",
                    mol_info,
                    atom_id_counter
                )
                
                mapping["molecules"].append(mol_entry)
                atom_id_counter += len(mol_info["atoms"])
                molecule_id_counter += 1
        
        # 处理阴离子
        for anion in job_data.get("anions", []):
            name = anion["name"]
            number = anion["number"]
            
            mol_info = self._get_molecule_info(name, "anion", work_path)
            
            for i in range(number):
                mol_entry = self._create_molecule_entry(
                    molecule_id_counter,
                    name,
                    "anion",
                    mol_info,
                    atom_id_counter
                )
                
                mapping["molecules"].append(mol_entry)
                atom_id_counter += len(mol_info["atoms"])
                molecule_id_counter += 1
        
        # 处理溶剂
        for solvent in job_data.get("solvents", []):
            name = solvent["name"]
            number = solvent["number"]
            
            mol_info = self._get_molecule_info(name, "solvent", work_path)
            
            for i in range(number):
                mol_entry = self._create_molecule_entry(
                    molecule_id_counter,
                    name,
                    "solvent",
                    mol_info,
                    atom_id_counter
                )
                
                mapping["molecules"].append(mol_entry)
                atom_id_counter += len(mol_info["atoms"])
                molecule_id_counter += 1
        
        # 从 data 文件中修正原子映射
        data_file = work_path / f"{work_path.name}.data"
        if data_file.exists():
            logger.info(f"Correcting atom mapping from data file: {data_file}")
            mapping = self._correct_mapping_from_data_file(mapping, data_file)
        else:
            logger.warning(f"Data file not found: {data_file}, atom mapping may be incorrect")

        # 保存到文件
        mapping_file = work_path / "atom_mapping.json"
        with open(mapping_file, 'w') as f:
            json.dump(mapping, f, indent=2)

        logger.info(f"Generated atom mapping: {len(mapping['molecules'])} molecules, "
                   f"{atom_id_counter - 1} atoms")

        return mapping

    def _get_molecule_info(
        self,
        name: str,
        mol_type: str,
        work_path: Path
    ) -> Dict:
        """
        获取分子信息（从 .lt 文件读取）

        Args:
            name: 分子名称
            mol_type: 分子类型（cation, anion, solvent）
            work_path: 工作目录

        Returns:
            分子信息字典
        """
        # 尝试从工作目录查找 LT 文件（优先）
        lt_file = work_path / f"{name}.lt"
        if lt_file.exists():
            try:
                logger.info(f"Parsing LT file for {name}: {lt_file}")
                return self.parser.parse_lt_file(lt_file, name)
            except Exception as e:
                logger.error(f"Failed to parse LT file {lt_file}: {e}")

        # 尝试从初始盐目录查找 LT 文件（阳离子/阴离子）
        lt_file = self.initial_salts_path / f"{name}.lt"
        if lt_file.exists():
            try:
                logger.info(f"Parsing LT file for {name}: {lt_file}")
                return self.parser.parse_lt_file(lt_file, name)
            except Exception as e:
                logger.error(f"Failed to parse LT file {lt_file}: {e}")

        # 如果找不到 LT 文件，记录警告并创建简单结构
        logger.warning(f"No LT file found for {name}, creating simple structure")
        return {
            "molecule_name": name,
            "atoms": [
                {
                    "index": 1,
                    "name": name,
                    "type": "unknown",
                    "element": name,
                    "charge": 0.0,
                    "mass": 1.0
                }
            ],
            "bonds": []
        }

    def _create_molecule_entry(
        self,
        molecule_id: int,
        molecule_name: str,
        molecule_type: str,
        mol_info: Dict,
        atom_id_start: int
    ) -> Dict:
        """
        创建单个分子的映射条目

        Args:
            molecule_id: 分子 ID
            molecule_name: 分子名称
            molecule_type: 分子类型
            mol_info: 分子信息
            atom_id_start: 起始原子 ID

        Returns:
            分子映射条目
        """
        mol_entry = {
            "molecule_id": molecule_id,
            "molecule_name": molecule_name,
            "molecule_type": molecule_type,
            "atoms": []
        }

        atom_id = atom_id_start
        for atom in mol_info["atoms"]:
            # 生成标签：优先使用环境信息
            environment = atom.get("environment", "")
            element = atom["element"]

            if environment and environment != element:
                # 有化学环境信息，使用 "分子名_元素(环境)" 格式
                label = f"{molecule_name}_{element}({environment})"
            else:
                # 没有环境信息，使用原子名称
                label = f"{molecule_name}_{atom['name']}"

            mol_entry["atoms"].append({
                "atom_id": atom_id,
                "atom_index": atom["index"],
                "atom_name": atom["name"],
                "label": label,
                "element": atom["element"],
                "environment": environment,
                "type": atom.get("type", "unknown"),
                "charge": atom.get("charge", 0.0),
                "mass": atom.get("mass", 1.0)
            })
            atom_id += 1

        return mol_entry

    def _correct_mapping_from_data_file(self, mapping: Dict, data_file: Path) -> Dict:
        """
        从 LAMMPS data 文件中读取实际的原子信息，修正 atom_mapping

        问题：Moltemplate 生成的 data 文件中，原子的顺序可能与 .lt 文件中的定义顺序不同
        解决：从 data 文件的 Atoms 部分读取实际的 atom_id, mol_id, atom_type，
             然后从 Masses 部分读取 atom_type 对应的元素

        Args:
            mapping: 初始的 atom_mapping（基于 .lt 文件生成）
            data_file: LAMMPS data 文件路径

        Returns:
            修正后的 atom_mapping
        """
        try:
            with open(data_file, 'r', encoding='utf-8') as f:
                content = f.read()

            # 1. 解析 Masses 部分，建立 type -> element 的映射
            type_to_element = {}
            masses_section = content.split('Masses')[1].split('Atoms')[0] if 'Masses' in content else ""
            for line in masses_section.strip().split('\n'):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        atom_type = int(parts[0])
                        mass = round(float(parts[1]))
                        element = self._mass_to_element(mass)
                        type_to_element[atom_type] = element
                    except ValueError:
                        continue

            logger.info(f"Parsed {len(type_to_element)} atom types from Masses section")

            # 2. 解析 Atoms 部分，建立 atom_id -> (mol_id, atom_type, element, charge) 的映射
            atom_data = {}  # {atom_id: {'mol_id': int, 'atom_type': int, 'element': str, 'charge': float}}
            atoms_section = content.split('Atoms')[1]
            if 'Bonds' in atoms_section:
                atoms_section = atoms_section.split('Bonds')[0]
            elif 'Velocities' in atoms_section:
                atoms_section = atoms_section.split('Velocities')[0]

            for line in atoms_section.strip().split('\n'):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        atom_id = int(parts[0])
                        mol_id = int(parts[1])
                        atom_type = int(parts[2])
                        charge = float(parts[3]) if len(parts) > 3 else 0.0
                        element = type_to_element.get(atom_type, 'X')
                        atom_data[atom_id] = {
                            'mol_id': mol_id,
                            'atom_type': atom_type,
                            'element': element,
                            'charge': charge
                        }
                    except ValueError:
                        continue

            logger.info(f"Parsed {len(atom_data)} atoms from Atoms section")

            # 3. 按照 mol_id 分组原子（优化性能）
            mol_id_to_atoms_data = {}
            for atom_id, data in atom_data.items():
                mol_id = data['mol_id']
                if mol_id not in mol_id_to_atoms_data:
                    mol_id_to_atoms_data[mol_id] = []
                mol_id_to_atoms_data[mol_id].append((atom_id, data))

            # 对每个分子的原子按 atom_id 排序
            for mol_id in mol_id_to_atoms_data:
                mol_id_to_atoms_data[mol_id].sort(key=lambda x: x[0])

            logger.info(f"Grouped atoms by molecule: {len(mol_id_to_atoms_data)} molecules")

            # 4. 修正 mapping 中的元素信息
            for mol in mapping['molecules']:
                mol_id = mol['molecule_id']

                # 找到该分子的所有原子（从 data 文件中）
                mol_atoms_in_data = mol_id_to_atoms_data.get(mol_id, [])

                if not mol_atoms_in_data:
                    logger.warning(f"No atoms found in data file for molecule {mol_id} ({mol['molecule_name']})")
                    continue

                # 更新 mapping 中的原子信息
                if len(mol_atoms_in_data) != len(mol['atoms']):
                    logger.warning(f"Atom count mismatch for molecule {mol_id}: "
                                 f"mapping has {len(mol['atoms'])}, data file has {len(mol_atoms_in_data)}")
                    continue

                # 使用电荷和元素匹配原子
                # 为 mapping 中的每个原子找到 data 文件中对应的原子
                matched_atoms = self._match_atoms_by_charge_and_element(
                    mol['atoms'], mol_atoms_in_data
                )

                # 更新 mapping
                for i, (atom_id, atom_info) in enumerate(matched_atoms):
                    if i < len(mol['atoms']):
                        # 更新元素信息（使用 data 文件中的实际元素）
                        actual_element = atom_info['element']
                        mol['atoms'][i]['element'] = actual_element

                        # 更新 atom_id（确保与 data 文件一致）
                        mol['atoms'][i]['atom_id'] = atom_id

                        # 更新 label
                        atom_name = mol['atoms'][i]['atom_name']
                        mol_name = mol['molecule_name']
                        mol['atoms'][i]['label'] = f"{mol_name}_{atom_name}"

            logger.info("Corrected atom mapping from data file")
            return mapping

        except Exception as e:
            logger.error(f"Failed to correct atom mapping from data file: {e}")
            import traceback
            traceback.print_exc()
            return mapping

    def _match_atoms_by_charge_and_element(
        self,
        mapping_atoms: List[Dict],
        data_atoms: List[tuple]
    ) -> List[tuple]:
        """
        使用电荷和元素匹配原子

        Args:
            mapping_atoms: mapping 中的原子列表（按 atom_index 排序）
            data_atoms: data 文件中的原子列表 [(atom_id, atom_info), ...]（按 atom_id 排序）

        Returns:
            匹配后的原子列表 [(atom_id, atom_info), ...]（按 mapping 中的顺序）
        """
        matched = []
        used_indices = set()

        for map_atom in mapping_atoms:
            expected_charge = map_atom.get('charge', 0.0)
            expected_element = map_atom['element']

            # 在 data_atoms 中查找匹配的原子
            best_match = None
            best_match_idx = None
            min_charge_diff = float('inf')

            for idx, (atom_id, atom_info) in enumerate(data_atoms):
                if idx in used_indices:
                    continue

                # 检查元素是否匹配
                if atom_info['element'] != expected_element:
                    continue

                # 检查电荷是否接近
                charge_diff = abs(atom_info['charge'] - expected_charge)
                if charge_diff < min_charge_diff:
                    min_charge_diff = charge_diff
                    best_match = (atom_id, atom_info)
                    best_match_idx = idx

            if best_match:
                matched.append(best_match)
                used_indices.add(best_match_idx)
            else:
                # 如果没有找到匹配的原子，使用第一个未使用的原子
                for idx, (atom_id, atom_info) in enumerate(data_atoms):
                    if idx not in used_indices:
                        matched.append((atom_id, atom_info))
                        used_indices.add(idx)
                        break

        return matched

    def _mass_to_element(self, mass: int) -> str:
        """根据原子质量返回元素符号"""
        mass_to_element_map = {
            1: 'H', 4: 'He', 7: 'Li', 9: 'Be', 11: 'B', 12: 'C',
            14: 'N', 16: 'O', 19: 'F', 20: 'Ne', 23: 'Na', 24: 'Mg',
            27: 'Al', 28: 'Si', 31: 'P', 32: 'S', 35: 'Cl', 39: 'K',
            40: 'Ca', 45: 'Sc', 48: 'Ti', 51: 'V', 52: 'Cr', 55: 'Mn',
            56: 'Fe', 59: 'Co', 64: 'Cu', 65: 'Zn', 70: 'Ga', 73: 'Ge',
            75: 'As', 79: 'Se'
        }
        return mass_to_element_map.get(mass, 'X')

