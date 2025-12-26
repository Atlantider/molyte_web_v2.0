/**
 * 统一设计系统
 * Modern Tech-Style Design System
 */

// ==================== 渐变系统 ====================
export const gradients = {
    // 主色调渐变
    primary: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)',
    primaryHover: 'linear-gradient(135deg, #5568d3 0%, #63408a 100%)',

    // 成功状态
    success: 'linear-gradient(135deg, #11998e 0%, #38ef7d 100%)',
    successHover: 'linear-gradient(135deg, #0d7a72 0%, #2bc765 100%)',

    // 警告状态
    warning: 'linear-gradient(135deg, #f093fb 0%, #f5576c 100%)',
    warningHover: 'linear-gradient(135deg, #d77ee0 0%, #db4456 100%)',

    // 信息状态
    info: 'linear-gradient(135deg, #4facfe 0%, #00f2fe 100%)',
    infoHover: 'linear-gradient(135deg, #3d8dd6 0%, #00cad4 100%)',

    // 深色渐变
    dark: 'linear-gradient(135deg, #2d3561 0%, #1a1f3a 100%)',
    darkHover: 'linear-gradient(135deg, #252a4f 0%, #151829 100%)',

    // 特殊效果渐变
    cosmic: 'linear-gradient(135deg, #fa709a 0%, #fee140 100%)',
    ocean: 'linear-gradient(135deg, #a8edea 0%, #fed6e3 100%)',
    sunset: 'linear-gradient(135deg, #ffecd2 0%, #fcb69f 100%)',
    aurora: 'linear-gradient(135deg, #a8edea 0%, #fed6e3 100%)',
};

// ==================== Glass Morphism ====================
export const glassMorphism = {
    light: {
        background: 'rgba(255, 255, 255, 0.85)',
        backdropFilter: 'blur(12px)',
        WebkitBackdropFilter: 'blur(12px)',
        border: '1px solid rgba(255, 255, 255, 0.25)',
        boxShadow: '0 8px 32px 0 rgba(31, 38, 135, 0.15)',
    },
    lightStrong: {
        background: 'rgba(255, 255, 255, 0.95)',
        backdropFilter: 'blur(16px)',
        WebkitBackdropFilter: 'blur(16px)',
        border: '1px solid rgba(255, 255, 255, 0.3)',
        boxShadow: '0 8px 32px 0 rgba(31, 38, 135, 0.2)',
    },
    dark: {
        background: 'rgba(30, 30, 30, 0.85)',
        backdropFilter: 'blur(12px)',
        WebkitBackdropFilter: 'blur(12px)',
        border: '1px solid rgba(255, 255, 255, 0.1)',
        boxShadow: '0 8px 32px 0 rgba(0, 0, 0, 0.3)',
    },
    darkStrong: {
        background: 'rgba(20, 20, 20, 0.95)',
        backdropFilter: 'blur(16px)',
        WebkitBackdropFilter: 'blur(16px)',
        border: '1px solid rgba(255, 255, 255, 0.15)',
        boxShadow: '0 8px 32px 0 rgba(0, 0, 0, 0.4)',
    },
};

// ==================== 阴影系统 ====================
export const shadows = {
    sm: '0 1px 2px 0 rgba(0, 0, 0, 0.05)',
    md: '0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)',
    lg: '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)',
    xl: '0 20px 25px -5px rgba(0, 0, 0, 0.1), 0 10px 10px -5px rgba(0, 0, 0, 0.04)',
    '2xl': '0 25px 50px -12px rgba(0, 0, 0, 0.25)',
    inner: 'inset 0 2px 4px 0 rgba(0, 0, 0, 0.06)',

    // Dark mode shadows
    darkSm: '0 1px 2px 0 rgba(0, 0, 0, 0.3)',
    darkMd: '0 4px 6px -1px rgba(0, 0, 0, 0.4), 0 2px 4px -1px rgba(0, 0, 0, 0.3)',
    darkLg: '0 10px 15px -3px rgba(0, 0, 0, 0.5), 0 4px 6px -2px rgba(0, 0, 0, 0.4)',
    darkXl: '0 20px 25px -5px rgba(0, 0, 0, 0.6), 0 10px 10px -5px rgba(0, 0, 0, 0.5)',

    // Glow effects
    glowPrimary: '0 0 20px rgba(102, 126, 234, 0.5)',
    glowSuccess: '0 0 20px rgba(56, 239, 125, 0.5)',
    glowWarning: '0 0 20px rgba(245, 87, 108, 0.5)',
    glowInfo: '0 0 20px rgba(79, 172, 254, 0.5)',
};

// ==================== 间距系统 ====================
export const spacing = {
    xs: '4px',
    sm: '8px',
    md: '12px',
    base: '16px',
    lg: '24px',
    xl: '32px',
    '2xl': '48px',
    '3xl': '64px',
    '4xl': '96px',
};

// ==================== 圆角系统 ====================
export const borderRadius = {
    none: '0',
    sm: '6px',
    md: '12px',
    lg: '16px',
    xl: '24px',
    '2xl': '32px',
    full: '9999px',
};

// ==================== 动画曲线 ====================
export const easings = {
    // 标准缓动
    standard: 'cubic-bezier(0.4, 0.0, 0.2, 1)',
    // 加速
    accelerate: 'cubic-bezier(0.4, 0.0, 1, 1)',
    // 减速
    decelerate: 'cubic-bezier(0.0, 0.0, 0.2, 1)',
    // 弹性
    spring: 'cubic-bezier(0.68, -0.55, 0.265, 1.55)',
    // 平滑
    smooth: 'cubic-bezier(0.25, 0.1, 0.25, 1)',
};

// ==================== 动画持续时间 ====================
export const durations = {
    fast: '150ms',
    normal: '300ms',
    slow: '500ms',
    slower: '800ms',
};

// ==================== 通用动画样式 ====================
export const animations = {
    // Hover上移
    hoverLift: {
        transition: `all ${durations.normal} ${easings.standard}`,
        '&:hover': {
            transform: 'translateY(-4px)',
            boxShadow: shadows.xl,
        },
    },

    // Hover缩放
    hoverScale: {
        transition: `transform ${durations.normal} ${easings.spring}`,
        '&:hover': {
            transform: 'scale(1.05)',
        },
    },

    // 淡入
    fadeIn: {
        animation: 'fadeIn 0.5s ease-in',
        '@keyframes fadeIn': {
            from: { opacity: 0 },
            to: { opacity: 1 },
        },
    },

    // 从下方滑入
    slideUp: {
        animation: 'slideUp 0.4s ease-out',
        '@keyframes slideUp': {
            from: {
                opacity: 0,
                transform: 'translateY(20px)',
            },
            to: {
                opacity: 1,
                transform: 'translateY(0)',
            },
        },
    },

    // 脉冲效果
    pulse: {
        animation: 'pulse 2s cubic-bezier(0.4, 0, 0.6, 1) infinite',
        '@keyframes pulse': {
            '0%, 100%': { opacity: 1 },
            '50%': { opacity: 0.5 },
        },
    },
};

// ==================== 卡片样式生成器 ====================
export const createCardStyle = (isDark: boolean, variant: 'glass' | 'solid' | 'gradient' = 'glass') => {
    const baseStyle = {
        borderRadius: borderRadius.lg,
        transition: `all ${durations.normal} ${easings.standard}`,
    };

    switch (variant) {
        case 'glass':
            return {
                ...baseStyle,
                ...(isDark ? glassMorphism.dark : glassMorphism.light),
            };

        case 'solid':
            return {
                ...baseStyle,
                background: isDark ? '#1a1a1a' : '#ffffff',
                border: `1px solid ${isDark ? 'rgba(255, 255, 255, 0.1)' : 'rgba(0, 0, 0, 0.06)'}`,
                boxShadow: isDark ? shadows.darkMd : shadows.md,
            };

        case 'gradient':
            return {
                ...baseStyle,
                background: isDark ? gradients.dark : gradients.primary,
                border: 'none',
                boxShadow: isDark ? shadows.darkLg : shadows.lg,
                color: '#ffffff',
            };

        default:
            return baseStyle;
    }
};

// ==================== 按钮样式生成器 ====================
export const createButtonStyle = (
    isDark: boolean,
    variant: 'primary' | 'success' | 'warning' | 'ghost' = 'primary'
) => {
    const baseStyle = {
        borderRadius: borderRadius.md,
        padding: `${spacing.md} ${spacing.lg}`,
        fontWeight: 600,
        transition: `all ${durations.normal} ${easings.standard}`,
        cursor: 'pointer',
        border: 'none',
    };

    switch (variant) {
        case 'primary':
            return {
                ...baseStyle,
                background: gradients.primary,
                color: '#ffffff',
                boxShadow: shadows.md,
                '&:hover': {
                    background: gradients.primaryHover,
                    boxShadow: shadows.lg,
                    transform: 'translateY(-2px)',
                },
            };

        case 'success':
            return {
                ...baseStyle,
                background: gradients.success,
                color: '#ffffff',
                boxShadow: shadows.md,
                '&:hover': {
                    background: gradients.successHover,
                    boxShadow: shadows.lg,
                    transform: 'translateY(-2px)',
                },
            };

        case 'warning':
            return {
                ...baseStyle,
                background: gradients.warning,
                color: '#ffffff',
                boxShadow: shadows.md,
                '&:hover': {
                    background: gradients.warningHover,
                    boxShadow: shadows.lg,
                    transform: 'translateY(-2px)',
                },
            };

        case 'ghost':
            return {
                ...baseStyle,
                background: 'transparent',
                color: isDark ? '#ffffff' : '#000000',
                border: `1px solid ${isDark ? 'rgba(255, 255, 255, 0.2)' : 'rgba(0, 0, 0, 0.1)'}`,
                '&:hover': {
                    background: isDark ? 'rgba(255, 255, 255, 0.1)' : 'rgba(0, 0, 0, 0.05)',
                    borderColor: isDark ? 'rgba(255, 255, 255, 0.3)' : 'rgba(0, 0, 0, 0.2)',
                },
            };

        default:
            return baseStyle;
    }
};

// ==================== 输入框样式生成器 ====================
export const createInputStyle = (isDark: boolean) => ({
    borderRadius: borderRadius.md,
    padding: `${spacing.md} ${spacing.base}`,
    border: `1px solid ${isDark ? 'rgba(255, 255, 255, 0.15)' : 'rgba(0, 0, 0, 0.1)'}`,
    background: isDark ? 'rgba(255, 255, 255, 0.05)' : 'rgba(255, 255, 255, 0.8)',
    transition: `all ${durations.normal} ${easings.standard}`,
    '&:focus': {
        outline: 'none',
        borderColor: '#667eea',
        boxShadow: '0 0 0 3px rgba(102, 126, 234, 0.1)',
        background: isDark ? 'rgba(255, 255, 255, 0.08)' : '#ffffff',
    },
    '&:hover': {
        borderColor: isDark ? 'rgba(255, 255, 255, 0.25)' : 'rgba(0, 0, 0, 0.2)'
    },
});
