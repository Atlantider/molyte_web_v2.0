/**
 * QC重复计算检查提示组件
 */
import React from 'react';
import { Alert } from 'antd';
import { CheckCircleOutlined } from '@ant-design/icons';
import type { DuplicateCheckResponse } from '../../api/qc';

interface DuplicateCheckAlertProps {
    result: DuplicateCheckResponse | null;
}

export function DuplicateCheckAlert({ result }: DuplicateCheckAlertProps) {
    if (!result || result.existing_count === 0) {
        return null;
    }

    return (
        <Alert
            type="success"
            showIcon
            icon={<CheckCircleOutlined />}
            style={{ marginBottom: 16 }}
            message={
                <span>
                    检测到 <strong>{result.existing_count}</strong> 个分子已有计算结果，
                    将直接复用，无需重复计算！
                </span>
            }
            description={
                result.new_count > 0 ? (
                    <span>另外 {result.new_count} 个分子将执行新计算。</span>
                ) : (
                    <span>所有QC计算都将复用已有结果，节省计算时间和资源！</span>
                )
            }
        />
    );
}
