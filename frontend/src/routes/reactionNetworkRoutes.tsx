/**
 * Frontend Route Configuration Update
 * 添加反应网络模块路由
 * 
 * 将以下路由添加到你的主路由配置文件中（通常是 App.tsx 或 routes.tsx）
 */

import ReactionNetworkJobs from './pages/ReactionNetworkJobs';
import ReactionNetworkCreate from './pages/ReactionNetworkCreate';
import ReactionNetworkDetail from './pages/ReactionNetworkDetail';

// 在路由配置中添加以下路由：
const reactionNetworkRoutes = [
    {
        path: '/reaction-network',
        element: <ReactionNetworkJobs />,
        meta: {
            title: '反应网络',
            icon: 'ExperimentOutlined',
            requiresAuth: true
        }
    },
    {
        path: '/reaction-network/create',
        element: <ReactionNetworkCreate />,
        meta: {
            title: '创建反应网络任务',
            requiresAuth: true
        }
    },
    {
        path: '/reaction-network/:id',
        element: <ReactionNetworkDetail />,
        meta: {
            title: '反应网络详情',
            requiresAuth: true
        }
    }
];

// 导航菜单配置（添加到主菜单）
const menuItem = {
    key: 'reaction-network',
    icon: <ExperimentOutlined />,
    label: '反应网络',
    path: '/reaction-network'
};

export { reactionNetworkRoutes, menuItem };
