-- 修改 qc_results 表的 smiles 字段为可空
-- 原因：从 cluster 中提取的 dimer/ligand 任务使用 XYZ 坐标，没有 SMILES

ALTER TABLE qc_results ALTER COLUMN smiles DROP NOT NULL;

-- 添加注释
COMMENT ON COLUMN qc_results.smiles IS '分子 SMILES 字符串（对于从 cluster 提取的任务可为空）';

