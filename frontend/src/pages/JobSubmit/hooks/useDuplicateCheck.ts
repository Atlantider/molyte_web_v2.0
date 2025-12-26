/**
 * QC重复计算检查Hook
 * 
 * 管理QC计算的重复检查逻辑
 */
import { useState, useCallback } from 'react';
import type { ElectrolyteSystem, MDJob } from '../../types';
import type { DuplicateCheckResponse } from '../../api/qc';
import { checkDuplicateCalculations } from '../../api/qc';

interface QCMoleculeParams {
    functional: string;
    basis_set: string;
    solvent_model: string;
}

interface UseDuplicateCheckOptions {
    job: MDJob | null;
    electrolyte: ElectrolyteSystem | null;
    moleculeParams: Record<string, QCMoleculeParams>;
}

export function useDuplicateCheck({ job, electrolyte, moleculeParams }: UseDuplicateCheckOptions) {
    const [checking, setChecking] = useState(false);
    const [result, setResult] = useState<DuplicateCheckResponse | null>(null);

    const checkDuplicates = useCallback(async () => {
        if (!job?.config?.qc_enabled || !electrolyte) {
            return null;
        }

        const molecules: any[] = [];
        const config = job.config;
        const solventModel = config.qc_solvent_model || 'pcm';
        const solventName = config.qc_solvent_name || 'water';
        const functional = config.qc_functional || 'B3LYP';
        const basisSet = config.qc_basis_set || '6-31G(d)';

        // 收集溶剂分子
        electrolyte.solvents?.forEach((s: any) => {
            const customParams = moleculeParams[`solvent_${electrolyte.solvents?.indexOf(s)}`];
            molecules.push({
                smiles: s.smiles,
                molecule_name: s.name,
                functional: customParams?.functional || functional,
                basis_set: customParams?.basis_set || basisSet,
                solvent_model: customParams?.solvent_model || solventModel,
                solvent_name: solventModel !== 'gas' ? solventName : undefined,
                charge: 0,
                spin_multiplicity: 1,
            });
        });

        // 收集阳离子
        electrolyte.cations?.forEach((c: any) => {
            const customParams = moleculeParams[`cation_${electrolyte.cations?.indexOf(c)}`];
            molecules.push({
                smiles: c.smiles,
                molecule_name: c.name,
                functional: customParams?.functional || functional,
                basis_set: customParams?.basis_set || basisSet,
                solvent_model: customParams?.solvent_model || solventModel,
                solvent_name: solventModel !== 'gas' ? solventName : undefined,
                charge: 1,
                spin_multiplicity: 1,
            });
        });

        // 收集阴离子
        electrolyte.anions?.forEach((a: any) => {
            const customParams = moleculeParams[`anion_${electrolyte.anions?.indexOf(a)}`];
            molecules.push({
                smiles: a.smiles,
                molecule_name: a.name,
                functional: customParams?.functional || functional,
                basis_set: customParams?.basis_set || basisSet,
                solvent_model: customParams?.solvent_model || solventModel,
                solvent_name: solventModel !== 'gas' ? solventName : undefined,
                charge: -1,
                spin_multiplicity: 1,
            });
        });

        if (molecules.length === 0) {
            return null;
        }

        try {
            setChecking(true);
            const checkResult = await checkDuplicateCalculations(molecules);
            setResult(checkResult);
            return checkResult;
        } catch (error) {
            console.error('检查重复计算失败:', error);
            return null;
        } finally {
            setChecking(false);
        }
    }, [job, electrolyte, moleculeParams]);

    return {
        checking,
        result,
        checkDuplicates,
        clearResult: () => setResult(null),
    };
}
