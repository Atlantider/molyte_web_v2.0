import React, { useState, useEffect } from 'react';
import { Input, Button, Form, message } from 'antd';
import { MobileOutlined } from '@ant-design/icons';
import { sendSmsCode } from '../api/sms';

interface PhoneVerificationInputProps {
  onPhoneChange?: (phone: string) => void;
  onCodeChange?: (code: string) => void;
  purpose?: 'register' | 'bind' | 'reset';
  disabled?: boolean;
}

/**
 * 手机验证码输入组件
 * 现代 SaaS 风格，带倒计时功能
 */
const PhoneVerificationInput: React.FC<PhoneVerificationInputProps> = ({
  onPhoneChange,
  onCodeChange,
  purpose = 'register',
  disabled = false,
}) => {
  const [phone, setPhone] = useState('');
  const [code, setCode] = useState('');
  const [countdown, setCountdown] = useState(0);
  const [sending, setSending] = useState(false);

  // 倒计时逻辑
  useEffect(() => {
    if (countdown > 0) {
      const timer = setTimeout(() => setCountdown(countdown - 1), 1000);
      return () => clearTimeout(timer);
    }
  }, [countdown]);

  // 手机号验证
  const isValidPhone = (phone: string) => {
    return /^1[3-9]\d{9}$/.test(phone);
  };

  // 发送验证码
  const handleSendCode = async () => {
    if (!isValidPhone(phone)) {
      message.error('请输入有效的手机号');
      return;
    }

    setSending(true);
    try {
      const response = await sendSmsCode({ phone, purpose });
      message.success(response.message || '验证码已发送');
      setCountdown(60); // 60秒倒计时

      // 开发环境下显示验证码
      if (response.code && import.meta.env.DEV) {
        message.info(`测试验证码: ${response.code}`, 10);
      }
    } catch (error: any) {
      message.error(error.response?.data?.detail || '发送失败，请稍后重试');
    } finally {
      setSending(false);
    }
  };

  // 手机号变化
  const handlePhoneChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const value = e.target.value;
    setPhone(value);
    onPhoneChange?.(value);
  };

  // 验证码变化
  const handleCodeChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const value = e.target.value;
    setCode(value);
    onCodeChange?.(value);
  };

  return (
    <div className="phone-verification-input">
      <Form.Item
        label="手机号"
        rules={[
          { required: true, message: '请输入手机号' },
          { pattern: /^1[3-9]\d{9}$/, message: '请输入有效的手机号' },
        ]}
      >
        <Input
          prefix={<MobileOutlined />}
          placeholder="请输入手机号"
          value={phone}
          onChange={handlePhoneChange}
          disabled={disabled}
          size="large"
        />
      </Form.Item>

      <Form.Item
        label="验证码"
        rules={[{ required: true, message: '请输入验证码' }]}
      >
        <Input.Group compact style={{ display: 'flex' }}>
          <Input
            placeholder="请输入验证码"
            value={code}
            onChange={handleCodeChange}
            disabled={disabled}
            size="large"
            style={{ flex: 1 }}
            maxLength={6}
          />
          <Button
            type="primary"
            size="large"
            onClick={handleSendCode}
            disabled={disabled || countdown > 0 || !isValidPhone(phone)}
            loading={sending}
            style={{ marginLeft: 8 }}
          >
            {countdown > 0 ? `${countdown}秒后重试` : '获取验证码'}
          </Button>
        </Input.Group>
      </Form.Item>
    </div>
  );
};

export default PhoneVerificationInput;

