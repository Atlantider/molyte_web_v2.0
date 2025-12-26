/**
 * 阴离子选择器 - 带有自动生成功能
 */
import React, { useState, useEffect } from 'react';
import {
  Select,
  Button,
  Space,
  Tooltip,
  Spin,
  message,
} from 'antd';
import { PlusOutlined, ReloadOutlined } from '@ant-design/icons';
import type { IonInfo } from '../types';
import { getAnionLibrary } from '../api/forcefield';
import AnionGenerationModal from './AnionGenerationModal';

interface AnionSelectorWithGenerationProps {
  availableAnions: IonInfo[];
  selectedAnions: IonInfo[];
  onAddAnion: (ionName: string) => void;
  onRefresh?: () => void;
}

export const AnionSelectorWithGeneration: React.FC<AnionSelectorWithGenerationProps> = ({
  availableAnions,
  selectedAnions,
  onAddAnion,
  onRefresh,
}) => {
  const [generationModalVisible, setGenerationModalVisible] = useState(false);
  const [refreshing, setRefreshing] = useState(false);

  const handleRefreshAnions = async () => {
    try {
      setRefreshing(true);
      // 重新加载阴离子库
      const library = await getAnionLibrary();
      message.success(`已刷新阴离子库，共 ${library.total} 个阴离子`);
      onRefresh?.();
    } catch (error) {
      console.error('刷新阴离子库失败:', error);
      message.error('刷新失败，请重试');
    } finally {
      setRefreshing(false);
    }
  };

  const handleGenerationSuccess = () => {
    // 生成成功后刷新阴离子列表
    handleRefreshAnions();
  };

  return (
    <>
      <Space direction="vertical" style={{ width: '100%' }} size="small">
        <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
          <Select
            placeholder="选择阴离子添加"
            style={{ flex: 1 }}
            onChange={onAddAnion}
            value={undefined}
            showSearch
            size="small"
            filterOption={(input, option) =>
              (option?.label as string)?.toLowerCase().includes(input.toLowerCase())
            }
          >
            {availableAnions
              .filter(ion => !selectedAnions.find(a => a.name === ion.name))
              .map(ion => (
                <Select.Option key={ion.name} value={ion.name}>
                  {ion.name} ({ion.charge})
                </Select.Option>
              ))}
          </Select>

          <Tooltip title="生成新阴离子势函数">
            <Button
              type="primary"
              icon={<PlusOutlined />}
              size="small"
              onClick={() => setGenerationModalVisible(true)}
            >
              新阴离子
            </Button>
          </Tooltip>

          <Tooltip title="刷新阴离子库">
            <Button
              icon={<ReloadOutlined />}
              size="small"
              loading={refreshing}
              onClick={handleRefreshAnions}
            />
          </Tooltip>
        </div>
      </Space>

      <AnionGenerationModal
        visible={generationModalVisible}
        onClose={() => setGenerationModalVisible(false)}
        onSuccess={handleGenerationSuccess}
      />
    </>
  );
};

export default AnionSelectorWithGeneration;

