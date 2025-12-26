/**
 * Anion Force Field Auto-Generation Modal
 */
import React, { useState, useEffect } from 'react';
import {
  Modal,
  Form,
  Input,
  Select,
  Button,
  message,
  Spin,
  Alert,
  Progress,
  Space,
  Divider,
  Tag,
  Collapse,
} from 'antd';
import { PlusOutlined, CheckCircleOutlined, CloseCircleOutlined } from '@ant-design/icons';
import {
  submitAnionGeneration,
  getAnionGenerationStatus,
  AnionGenerationRequest,
  AnionGenerationStatusResponse,
} from '../api/forcefield';

interface AnionGenerationModalProps {
  visible: boolean;
  onClose: () => void;
  onSuccess?: () => void;
}

export const AnionGenerationModal: React.FC<AnionGenerationModalProps> = ({
  visible,
  onClose,
  onSuccess,
}) => {
  const [form] = Form.useForm();
  const [loading, setLoading] = useState(false);
  const [jobId, setJobId] = useState<string | null>(null);
  const [jobStatus, setJobStatus] = useState<AnionGenerationStatusResponse | null>(null);
  const [statusLoading, setStatusLoading] = useState(false);
  const [pollInterval, setPollInterval] = useState<ReturnType<typeof setTimeout> | null>(null);
  const [savedFormData, setSavedFormData] = useState<any>(null);

  // Poll job status
  useEffect(() => {
    if (!jobId) return;

    const pollStatus = async () => {
      try {
        setStatusLoading(true);
        const status = await getAnionGenerationStatus(jobId);
        setJobStatus(status);

        // Stop polling if job is completed or failed
        if (status.status === 'success' || status.status === 'failed') {
          if (pollInterval) {
            clearInterval(pollInterval);
            setPollInterval(null);
          }
          setStatusLoading(false);

          if (status.status === 'success') {
            message.success('阴离子力场生成成功！');
            setTimeout(() => {
              onSuccess?.();
              handleClose();
            }, 1500);
          }
        }
      } catch (error) {
        console.error('Failed to get job status:', error);
      }
    };

    // Initial poll
    pollStatus();

    // Set up interval for subsequent polls
    if (!pollInterval && jobStatus?.status !== 'success' && jobStatus?.status !== 'failed') {
      const interval = setInterval(pollStatus, 3000); // Poll every 3 seconds
      setPollInterval(interval);
    }

    return () => {
      if (pollInterval) {
        clearInterval(pollInterval);
      }
    };
  }, [jobId, jobStatus?.status]);

  const handleSave = () => {
    const values = form.getFieldsValue();
    setSavedFormData(values);
    message.success('表单数据已保存');
  };

  const handleLoad = () => {
    if (savedFormData) {
      form.setFieldsValue(savedFormData);
      message.success('表单数据已加载');
    } else {
      message.warning('没有保存的数据');
    }
  };

  const handleSubmit = async (values: any) => {
    try {
      setLoading(true);

      const request: AnionGenerationRequest = {
        anion_name: values.anion_name,
        display_name: values.display_name,
        charge: values.charge || -1,
        identifier_type: values.identifier_type,
        identifier_value: values.identifier_value,
      };

      const response = await submitAnionGeneration(request);
      setJobId(response.job_id);
      message.info('任务已提交，正在处理中...');
    } catch (error: any) {
      const errorMsg = error.response?.data?.detail || '提交失败，请重试';
      message.error(errorMsg);
    } finally {
      setLoading(false);
    }
  };

  const handleClose = () => {
    if (pollInterval) {
      clearInterval(pollInterval);
      setPollInterval(null);
    }
    setJobId(null);
    setJobStatus(null);
    form.resetFields();
    onClose();
  };

  const getStatusColor = (status: string) => {
    switch (status) {
      case 'pending':
        return 'processing';
      case 'running':
        return 'processing';
      case 'success':
        return 'success';
      case 'failed':
        return 'error';
      default:
        return 'default';
    }
  };

  const getStatusText = (status: string) => {
    switch (status) {
      case 'pending':
        return '等待处理';
      case 'running':
        return '处理中';
      case 'success':
        return '成功';
      case 'failed':
        return '失败';
      default:
        return status;
    }
  };

  return (
    <Modal
      title="生成新阴离子势函数"
      open={visible}
      onCancel={handleClose}
      footer={null}
      width={600}
      destroyOnClose
    >
      {!jobId ? (
        // Form for submitting new job
        <Form
          form={form}
          layout="vertical"
          onFinish={handleSubmit}
          disabled={loading}
        >
          <Form.Item
            label="阴离子短名"
            name="anion_name"
            rules={[
              { required: true, message: '请输入阴离子短名' },
              { min: 1, max: 50, message: '长度应在 1-50 字符之间' },
            ]}
            tooltip="如 FSI, Cl, BF4 等，用作文件名和数据库 key"
          >
            <Input placeholder="例如: FSI" />
          </Form.Item>

          <Form.Item
            label="显示名称"
            name="display_name"
            tooltip="阴离子的完整名称，可选"
          >
            <Input placeholder="例如: bis(fluorosulfonyl)imide" />
          </Form.Item>

          <Form.Item
            label="电荷"
            name="charge"
            initialValue={-1}
          >
            <Input type="number" placeholder="-1" />
          </Form.Item>

          <Form.Item
            label="身份类型"
            name="identifier_type"
            rules={[{ required: true, message: '请选择身份类型' }]}
          >
            <Select placeholder="选择输入方式">
              <Select.Option value="smiles">SMILES</Select.Option>
              <Select.Option value="cas">CAS 编号</Select.Option>
            </Select>
          </Form.Item>

          <Form.Item
            label="身份值"
            name="identifier_value"
            rules={[{ required: true, message: '请输入身份值' }]}
            tooltip="SMILES 字符串或 CAS 编号"
          >
            <Input.TextArea
              placeholder="例如: N(S(=O)(=O)F)S(=O)(=O)F (SMILES) 或 372-64-0 (CAS)"
              rows={3}
            />
          </Form.Item>

          <Form.Item>
            <Space direction="vertical" style={{ width: '100%' }}>
              <Space style={{ width: '100%' }}>
                <Button onClick={handleSave} style={{ flex: 1 }}>
                  保存表单
                </Button>
                <Button onClick={handleLoad} style={{ flex: 1 }}>
                  加载表单
                </Button>
              </Space>
              <Button type="primary" htmlType="submit" loading={loading} block>
                <PlusOutlined /> 提交生成任务
              </Button>
            </Space>
          </Form.Item>
        </Form>
      ) : (
        // Job status display
        <Spin spinning={statusLoading}>
          <Space direction="vertical" style={{ width: '100%' }} size="large">
            <div>
              <div style={{ marginBottom: '8px' }}>
                <strong>任务 ID:</strong> {jobId}
              </div>
              <div style={{ marginBottom: '8px' }}>
                <strong>状态:</strong>{' '}
                <Tag color={getStatusColor(jobStatus?.status || 'pending')}>
                  {getStatusText(jobStatus?.status || 'pending')}
                </Tag>
              </div>
              {jobStatus?.anion_key && (
                <div style={{ marginBottom: '8px' }}>
                  <strong>阴离子:</strong> {jobStatus.anion_key}
                </div>
              )}
            </div>

            {jobStatus?.status === 'running' && (
              <Progress percent={50} status="active" />
            )}

            {jobStatus?.message && (
              <Alert
                message={jobStatus.message}
                type={jobStatus.status === 'failed' ? 'error' : 'info'}
                showIcon
              />
            )}

            {jobStatus?.status === 'success' && jobStatus?.files && (
              <Collapse
                items={[
                  {
                    key: '1',
                    label: '生成的文件',
                    children: (
                      <Space direction="vertical" style={{ width: '100%' }}>
                        <div>
                          <strong>.lt 文件:</strong>
                          <br />
                          <code>{jobStatus.files.lt_path}</code>
                        </div>
                        <div>
                          <strong>.pdb 文件:</strong>
                          <br />
                          <code>{jobStatus.files.pdb_path}</code>
                        </div>
                      </Space>
                    ),
                  },
                ]}
              />
            )}

            <Divider />

            <Space style={{ width: '100%', justifyContent: 'flex-end' }}>
              {jobStatus?.status === 'success' || jobStatus?.status === 'failed' ? (
                <Button type="primary" onClick={handleClose}>
                  关闭
                </Button>
              ) : (
                <Button onClick={handleClose}>
                  返回
                </Button>
              )}
            </Space>
          </Space>
        </Spin>
      )}
    </Modal>
  );
};

export default AnionGenerationModal;

