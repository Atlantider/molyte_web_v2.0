/**
 * 主应用组件
 */
import { useEffect } from 'react';
import { RouterProvider } from 'react-router-dom';
import { ConfigProvider } from 'antd';
import zhCN from 'antd/locale/zh_CN';
import router from './router';
import { useThemeStore, getAntdTheme } from './stores/themeStore';
import { ErrorBoundary } from './components/ErrorBoundary';

function App() {
  const { mode } = useThemeStore();
  const themeConfig = getAntdTheme(mode);

  // 更新 body 背景色和主题属性 - Robust Sync
  useEffect(() => {
    // 强制同步 data-theme 属性，确保 DOM 与 Store 状态一致
    const targetMode = mode === 'dark' ? 'dark' : 'light';

    // 1. 设置 html 属性 (index.css 依赖此属性)
    document.documentElement.setAttribute('data-theme', targetMode);

    // 2. 设置 body 属性 (部分全局样式可能依赖此)
    document.body.setAttribute('data-theme', targetMode);
    document.body.style.backgroundColor = mode === 'dark' ? '#141414' : '#F5F7FA';
    document.body.style.color = mode === 'dark' ? '#E8E8E8' : '#333333'; // 强制设置 body 文字颜色
    document.body.style.transition = 'background-color 0.3s, color 0.3s';

    // Debug output
    console.log(`[Theme] Mode switched to: ${mode} (Attribute set to: ${targetMode})`);
  }, [mode]);

  return (
    <ErrorBoundary>
      <ConfigProvider locale={zhCN} theme={themeConfig}>
        <RouterProvider router={router} />
      </ConfigProvider>
    </ErrorBoundary>
  );
}

export default App;

