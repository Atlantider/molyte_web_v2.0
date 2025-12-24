#!/usr/bin/env python3
"""调试 PF6 的电子数计算"""

from rdkit import Chem

# PF6 的 SMILES
pf6_smiles = "F[P-](F)(F)(F)(F)F"

mol = Chem.MolFromSmiles(pf6_smiles)
if mol:
    print(f"PF6 SMILES: {pf6_smiles}")
    print(f"原子数: {mol.GetNumAtoms()}")
    
    # 不添加氢原子
    print("\n不添加氢原子:")
    total_electrons = 0
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        symbol = atom.GetSymbol()
        print(f"  {symbol}: 原子序数={atomic_num}")
        total_electrons += atomic_num
    print(f"总电子数: {total_electrons}")
    
    # 添加氢原子
    mol_with_h = Chem.AddHs(mol)
    print(f"\n添加氢原子后:")
    print(f"原子数: {mol_with_h.GetNumAtoms()}")
    total_electrons_with_h = 0
    for atom in mol_with_h.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        symbol = atom.GetSymbol()
        print(f"  {symbol}: 原子序数={atomic_num}")
        total_electrons_with_h += atomic_num
    print(f"总电子数: {total_electrons_with_h}")
    
    # 计算形式电荷
    formal_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    print(f"\n形式电荷: {formal_charge}")
    
    # 计算实际电子数
    actual_electrons = total_electrons_with_h - formal_charge
    print(f"实际电子数 (总电子数 - 形式电荷): {actual_electrons}")
    
    # K+ 的电子数
    print(f"\nK+ 的电子数: 19")
    print(f"K+ + 3×PF6-:")
    print(f"  总电子数: 19 + 3×{total_electrons_with_h} = {19 + 3*total_electrons_with_h}")
    print(f"  总电荷: 1 + 3×(-1) = -2")
    print(f"  实际电子数: {19 + 3*total_electrons_with_h} - (-2) = {19 + 3*total_electrons_with_h + 2}")
    print(f"  自旋多重度: {'偶数 → 1' if (19 + 3*total_electrons_with_h + 2) % 2 == 0 else '奇数 → 2'}")
else:
    print("无法解析 SMILES")

