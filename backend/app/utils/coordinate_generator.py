"""
渐进式3D坐标生成模块

提供多层次的坐标生成策略，避免直接使用随机坐标。

策略优先级:
1. 标准距离几何算法 (最可靠)
2. ETKDG v3 改进算法
3. UFF 力场优化
4. XTB 半经验优化 (如果可用)
5. 验证的随机坐标 (最后手段，需严格检查)
"""
import logging
from pathlib import Path
from typing import Optional, Dict, Any, Tuple
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class CoordinateGenerationResult:
    """坐标生成结果"""
    xyz_content: str
    source: str  # "standard_dg", "etkdg", "uff", "xtb", "random", "user"
    quality: str  # "excellent", "good", "acceptable", "poor"
    warnings: list[str]
    min_distance: float
    max_distance: float
    energy: Optional[float] = None


class CoordinateGenerator:
    """3D坐标生成器"""
    
    def __init__(self, enable_xtb: bool = True):
        """
        初始化坐标生成器
        
        Args:
            enable_xtb: 是否启用XTB优化
        """
        self.enable_xtb = enable_xtb
        self._check_dependencies()
    
    def _check_dependencies(self):
        """检查依赖库"""
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            self.rdkit_available = True
        except ImportError:
            logger.error("RDKit not available")
            self.rdkit_available = False
        
        if self.enable_xtb:
            try:
                from app.services.xtb_optimizer import XTBOptimizer
                self.xtb_available = True
            except ImportError:
                logger.debug("XTB optimizer not available")
                self.xtb_available = False
        else:
            self.xtb_available = False
    
    def generate_coordinates(
        self,
        smiles: str,
        molecule_name: Optional[str] = None,
        charge: int = 0,
        multiplicity: int = 1
    ) -> CoordinateGenerationResult:
        """
        使用渐进式策略生成3D坐标
        
        Args:
            smiles: SMILES字符串
            molecule_name: 分子名称
            charge: 电荷
            multiplicity: 自旋多重度
            
        Returns:
            CoordinateGenerationResult
            
        Raises:
            ValueError: 所有方法都失败时
        """
        if not self.rdkit_available:
            raise RuntimeError("RDKit is required for coordinate generation")
        
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        mol = Chem.AddHs(mol)
        name = molecule_name or smiles
        
        # 策略1: 标准距离几何
        result = self._try_standard_dg(mol, name)
        if result:
            logger.info(f"✓ Standard DG succeeded for {name}")
            return result
        
        # 策略2: ETKDG v3
        result = self._try_etkdg(mol, name)
        if result:
            logger.info(f"✓ ETKDG succeeded for {name}")
            return result
        
        # 策略3: UFF力场
        result = self._try_uff(mol, name)
        if result:
            logger.info(f"✓ UFF succeeded for {name}")
            return result
        
        # 策略4: XTB优化
        if self.xtb_available:
            result = self._try_xtb(mol, name, charge, multiplicity)
            if result:
                logger.warning(f"⚠ XTB succeeded (after crude coords) for {name}")
                return result
        
        # 策略5: 验证的随机坐标 (最后手段)
        result = self._try_validated_random(mol, name)
        if result:
            logger.error(f"⚠⚠ Using RANDOM coordinates for {name} - Structure may be unreliable!")
            return result
        
        raise ValueError(
            f"Failed to generate valid 3D coordinates for {smiles}\n"
            f"Tried: Standard DG, ETKDG, UFF, XTB, Random (10 attempts)\n"
            f"Suggestion: Provide XYZ coordinates manually or check SMILES validity"
        )
    
    def _try_standard_dg(self, mol, name: str) -> Optional[CoordinateGenerationResult]:
        """策略1: 标准距离几何算法"""
        try:
            from rdkit.Chem import AllChem
            
            result = AllChem.EmbedMolecule(mol, randomSeed=42)
            if result != 0:
                return None
            
            # MMFF优化
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            
            xyz = self._mol_to_xyz(mol, name)
            check = self._validate_structure(xyz)
            
            if check['is_valid']:
                return CoordinateGenerationResult(
                    xyz_content=xyz,
                    source="standard_dg",
                    quality="excellent",
                    warnings=[],
                    min_distance=check['min_distance'],
                    max_distance=check['max_distance']
                )
        except Exception as e:
            logger.debug(f"Standard DG failed: {e}")
        
        return None
    
    def _try_etkdg(self, mol, name: str) -> Optional[CoordinateGenerationResult]:
        """策略2: ETKDG v3算法"""
        try:
            from rdkit.Chem import AllChem
            
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            params.useRandomCoords = False  # 明确禁用随机
            
            result = AllChem.EmbedMolecule(mol, params)
            if result != 0:
                return None
            
            # 更严格的MMFF优化
            AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
            
            xyz = self._mol_to_xyz(mol, name)
            check = self._validate_structure(xyz)
            
            if check['is_valid']:
                return CoordinateGenerationResult(
                    xyz_content=xyz,
                    source="etkdg",
                    quality="excellent",
                    warnings=[],
                    min_distance=check['min_distance'],
                    max_distance=check['max_distance']
                )
        except Exception as e:
            logger.debug(f"ETKDG failed: {e}")
        
        return None
    
    def _try_uff(self, mol, name: str) -> Optional[CoordinateGenerationResult]:
        """策略3: UFF力场优化"""
        try:
            from rdkit.Chem import AllChem
            
            # 先用标准方法生成初始坐标
            result = AllChem.EmbedMolecule(mol, randomSeed=42)
            if result != 0:
                return None
            
            # 使用UFF力场优化
            AllChem.UFFOptimizeMolecule(mol, maxIters=1000)
            
            xyz = self._mol_to_xyz(mol, name)
            check = self._validate_structure(xyz)
            
            if check['is_valid']:
                return CoordinateGenerationResult(
                    xyz_content=xyz,
                    source="uff",
                    quality="good",
                    warnings=["Used UFF force field"],
                    min_distance=check['min_distance'],
                    max_distance=check['max_distance']
                )
        except Exception as e:
            logger.debug(f"UFF failed: {e}")
        
        return None
    
    def _try_xtb(
        self,
        mol,
        name: str,
        charge: int,
        multiplicity: int
    ) -> Optional[CoordinateGenerationResult]:
        """策略4: XTB半经验优化"""
        try:
            from rdkit.Chem import AllChem
            from app.services.xtb_optimizer import XTBOptimizer
            import tempfile
            
            # 先获取粗略坐标(即使是随机的)
            result = AllChem.EmbedMolecule(mol, randomSeed=42)
            if result == -1:
                # 如果标准失败，尝试随机
                result = AllChem.EmbedMolecule(
                    mol,
                    useRandomCoords=True,
                    randomSeed=42,
                    maxAttempts=100
                )
                if result == -1:
                    return None
            
            crude_xyz = self._mol_to_xyz(mol, name)
            
            # 使用XTB优化
            with tempfile.TemporaryDirectory() as tmpdir:
                optimizer = XTBOptimizer(
                    charge=charge,
                    multiplicity=multiplicity
                )
                
                xtb_result = optimizer.optimize_xyz(
                    crude_xyz,
                    Path(tmpdir),
                    convergence='normal'
                )
                
                if xtb_result and xtb_result.get('converged'):
                    xyz = xtb_result['optimized_xyz']
                    check = self._validate_structure(xyz)
                    
                    if check['is_valid']:
                        return CoordinateGenerationResult(
                            xyz_content=xyz,
                            source="xtb",
                            quality="good",
                            warnings=["Used XTB pre-optimization from crude coords"],
                            min_distance=check['min_distance'],
                            max_distance=check['max_distance'],
                            energy=xtb_result.get('energy')
                        )
        except Exception as e:
            logger.debug(f"XTB failed: {e}")
        
        return None
    
    def _try_validated_random(self, mol, name: str) -> Optional[CoordinateGenerationResult]:
        """策略5: 验证的随机坐标 (最后手段)"""
        from rdkit.Chem import AllChem
        
        warnings = []
        
        # 尝试10次随机生成
        for attempt in range(10):
            try:
                result = AllChem.EmbedMolecule(
                    mol,
                    useRandomCoords=True,
                    randomSeed=42 + attempt,  # 不同种子
                    maxAttempts=100
                )
                
                if result != 0:
                    continue
                
                # 强力优化
                try:
                    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
                except:
                    AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
                
                xyz = self._mol_to_xyz(mol, name)
                check = self._validate_structure(xyz)
                
                # 严格的质量要求
                if check['is_valid'] and check['min_distance'] > 0.8:
                    warnings.append(f"RANDOM coordinates used (attempt {attempt + 1}/10)")
                    warnings.append(f"Min atom distance: {check['min_distance']:.2f} Å")
                    warnings.append("Structure may be unreliable - manual verification recommended")
                    
                    return CoordinateGenerationResult(
                        xyz_content=xyz,
                        source="random",
                        quality="poor",
                        warnings=warnings,
                        min_distance=check['min_distance'],
                        max_distance=check['max_distance']
                    )
            except Exception as e:
                logger.debug(f"Random attempt {attempt + 1} failed: {e}")
                continue
        
        return None
    
    def _mol_to_xyz(self, mol, name: str) -> str:
        """将RDKit分子对象转换为XYZ格式"""
        conf = mol.GetConformer()
        num_atoms = mol.GetNumAtoms()
        
        xyz_lines = [str(num_atoms), name or "molecule"]
        
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            symbol = atom.GetSymbol()
            xyz_lines.append(f"{symbol} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")
        
        return '\n'.join(xyz_lines)
    
    def _validate_structure(self, xyz_content: str) -> Dict[str, Any]:
        """验证结构质量"""
        try:
            # 使用qc_safety的检查函数
            from deployment.qc_safety import check_structure
            result = check_structure(xyz_content)
            
            return {
                'is_valid': result.is_valid,
                'min_distance': result.min_distance,
                'max_distance': result.max_distance,
                'warnings': result.warnings,
                'errors': result.errors
            }
        except Exception as e:
            logger.warning(f"Structure validation failed: {e}")
            # 简单验证
            return self._simple_validation(xyz_content)
    
    def _simple_validation(self, xyz_content: str) -> Dict[str, Any]:
        """简单的结构验证"""
        import math
        
        lines = xyz_content.strip().split('\n')
        atoms = []
        
        for line in lines[2:]:  # 跳过前两行
            parts = line.split()
            if len(parts) >= 4:
                try:
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                    atoms.append((x, y, z))
                except:
                    continue
        
        if len(atoms) < 2:
            return {'is_valid': False, 'min_distance': 0, 'max_distance': 0}
        
        # 计算最小距离
        min_dist = float('inf')
        max_dist = 0
        
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                dx = atoms[i][0] - atoms[j][0]
                dy = atoms[i][1] - atoms[j][1]
                dz = atoms[i][2] - atoms[j][2]
                dist = math.sqrt(dx*dx + dy*dy + dz*dz)
                min_dist = min(min_dist, dist)
                max_dist = max(max_dist, dist)
        
        # 原子间距至少0.5 Å
        is_valid = min_dist > 0.5
        
        return {
            'is_valid': is_valid,
            'min_distance': min_dist,
            'max_distance': max_dist,
            'warnings': [] if is_valid else [f"Min distance too small: {min_dist:.2f} Å"],
            'errors': []
        }


# 便利函数
def generate_3d_coordinates(
    smiles: str,
    molecule_name: Optional[str] = None,
    charge: int = 0,
    multiplicity: int = 1,
    enable_xtb: bool = True
) -> CoordinateGenerationResult:
    """
    生成3D坐标的便利函数
    
    Args:
        smiles: SMILES字符串
        molecule_name: 分子名称
        charge: 电荷
        multiplicity: 自旋多重度
        enable_xtb: 是否启用XTB
        
    Returns:
        CoordinateGenerationResult
        
    Example:
        >>> result = generate_3d_coordinates("CCO", "ethanol")
        >>> print(result.xyz_content)
        >>> print(f"Source: {result.source}, Quality: {result.quality}")
    """
    generator = CoordinateGenerator(enable_xtb=enable_xtb)
    return generator.generate_coordinates(
        smiles=smiles,
        molecule_name=molecule_name,
        charge=charge,
        multiplicity=multiplicity
    )
