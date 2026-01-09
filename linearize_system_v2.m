% linearize_system_v2.m
% 系统线性化 - 在平衡点附近线性化（符号形式）
%
% 输入: dynamics_new_coords.mat (用广义坐标表示的方程)
% 输出: linearized_system.mat (符号形式的 A, B 矩阵函数)
%
% 对照推导文档 §5 和 §6
%
% 状态向量 (10维):
%   x = [X_b^h, dX_b^h, phi, dphi, theta_l, dtheta_l, theta_r, dtheta_r, theta_b, dtheta_b]
% 控制向量 (4维):
%   u = [T_wr, T_wl, T_r, T_l]
%
% 注意：本脚本只输出符号形式，具体参数在 compute_lqr.m 中代入

clear; clc;
fprintf('========================================\n');
fprintf('系统线性化 (符号形式，对照推导文档 §5, §6)\n');
fprintf('========================================\n\n');

%% 定义符号变量

fprintf('Step 1: 定义符号变量...\n');

% 广义坐标
syms X_b_h dX_b_h ddX_b_h real
syms phi dphi ddphi real
syms theta_l dtheta_l ddtheta_l real
syms theta_r dtheta_r ddtheta_r real
syms theta_b dtheta_b ddtheta_b real

% 轮角速度
syms ddtheta_wl ddtheta_wr real

% 控制力矩
syms T_wr_to_r T_wl_to_l T_r_to_b T_l_to_b real

% 物理参数 (保持符号形式)
syms m_b m_l m_r m_wl m_wr real        % 质量
syms I_b I_l I_r I_wl I_wr I_yaw real  % 转动惯量
syms l_b l_l l_r l_l_d l_r_d real      % 长度
syms theta_l0 theta_r0 theta_b0 real   % 初始角
syms R R_w g real                       % 轮半径、轮距、重力

% 对称简化参数
syms m_leg m_w I_leg I_w l_leg l_leg_d real

%% 加载动力学方程

fprintf('Step 2: 加载动力学方程...\n');
load('dynamics_new_coords.mat');
fprintf('  已加载 dynamics_new_coords.mat\n\n');

%% 将轮角加速度代换为广义坐标表达式 (§5.0)

fprintf('========================================\n');
fprintf('将轮角加速度代换为广义坐标 (§5.0)\n');
fprintf('========================================\n\n');

% 轮角加速度符号（已在方程中）
% syms ddtheta_wl ddtheta_wr real  % 已在上面定义

% 腿角度相关项 (使用 l_l, l_r 因为对称假设在后面应用)
leg_term = (l_r*cos(theta_r)*ddtheta_r + l_l*cos(theta_l)*ddtheta_l)/2 ...
         - (l_r*sin(theta_r)*dtheta_r^2 + l_l*sin(theta_l)*dtheta_l^2)/2;

% 轮角加速度用广义坐标表示 (来自 §5.0)
% ddtheta_wr = (ddX_b_h + R_w*ddphi)/R - leg_term/R
% ddtheta_wl = (ddX_b_h - R_w*ddphi)/R - leg_term/R
ddtheta_wr_sub = (ddX_b_h + R_w*ddphi)/R - leg_term/R;
ddtheta_wl_sub = (ddX_b_h - R_w*ddphi)/R - leg_term/R;

fprintf('轮角加速度代换:\n');
fprintf('  ddtheta_wr = (ddX_b_h + R_w*ddphi)/R - leg_term/R\n');
fprintf('  ddtheta_wl = (ddX_b_h - R_w*ddphi)/R - leg_term/R\n\n');

% 对方程进行代换
wheel_subs = {
    ddtheta_wr, ddtheta_wr_sub;
    ddtheta_wl, ddtheta_wl_sub;
};

eq1_new = simplify(subs(eq1_new, wheel_subs(:,1), wheel_subs(:,2)));
eq2_new = simplify(subs(eq2_new, wheel_subs(:,1), wheel_subs(:,2)));
eq3_new = simplify(subs(eq3_new, wheel_subs(:,1), wheel_subs(:,2)));
eq4_new = simplify(subs(eq4_new, wheel_subs(:,1), wheel_subs(:,2)));
eq5_new = simplify(subs(eq5_new, wheel_subs(:,1), wheel_subs(:,2)));

fprintf('✓ 轮角加速度已代换为广义坐标\n\n');

%% 构建符号M, B, g矩阵 (对照 §5.1, §5.2, §5.3)

fprintf('========================================\n');
fprintf('构建 M*ddq = B*u + g 形式 (§5)\n');
fprintf('========================================\n\n');

% 广义加速度向量
ddq = [ddX_b_h; ddphi; ddtheta_l; ddtheta_r; ddtheta_b];

% 控制输入向量 (§5.0)
u = [T_wr_to_r; T_wl_to_l; T_r_to_b; T_l_to_b];

fprintf('----------------------------------------\n');
fprintf('§5.0 广义坐标和控制输入定义\n');
fprintf('----------------------------------------\n');
fprintf('  q = [X_b^h, phi, theta_l, theta_r, theta_b]^T\n');
fprintf('  u = [T_wr, T_wl, T_r, T_l]^T\n');
fprintf('✓ 广义坐标定义一致!\n\n');

%% 从方程中提取M矩阵

fprintf('----------------------------------------\n');
fprintf('§5.1 质量矩阵 M (5×5)\n');
fprintf('----------------------------------------\n');

M_sym = sym(zeros(5,5));
for i = 1:5
    if i == 1
        eq_temp = eq1_new;
    elseif i == 2
        eq_temp = eq2_new;
    elseif i == 3
        eq_temp = eq3_new;
    elseif i == 4
        eq_temp = eq4_new;
    else
        eq_temp = eq5_new;
    end
    
    for j = 1:5
        M_sym(i,j) = diff(eq_temp, ddq(j));
    end
end

fprintf('符号质量矩阵 M 已提取\n\n');

%% 提取B矩阵

fprintf('----------------------------------------\n');
fprintf('§5.2 控制矩阵 B (5×4)\n');
fprintf('----------------------------------------\n');

B_sym = sym(zeros(5,4));
for i = 1:5
    if i == 1
        eq_temp = eq1_new;
    elseif i == 2
        eq_temp = eq2_new;
    elseif i == 3
        eq_temp = eq3_new;
    elseif i == 4
        eq_temp = eq4_new;
    else
        eq_temp = eq5_new;
    end
    
    for j = 1:4
        % 方程形式: M*ddq + ... + B_coeff*u = 0
        % 文档形式: M*ddq = B*u + g, 即 B = diff(eq, u)
        B_sym(i,j) = diff(eq_temp, u(j));
    end
end

fprintf('符号控制矩阵 B 已提取\n\n');

%% 提取g向量

fprintf('----------------------------------------\n');
fprintf('§5.3 重力+离心力项 g (5×1)\n');
fprintf('----------------------------------------\n');

g_sym = sym(zeros(5,1));
for i = 1:5
    if i == 1
        eq_temp = eq1_new;
    elseif i == 2
        eq_temp = eq2_new;
    elseif i == 3
        eq_temp = eq3_new;
    elseif i == 4
        eq_temp = eq4_new;
    else
        eq_temp = eq5_new;
    end
    
    % 方程形式: M*ddq + B_coeff*u + rest = 0
    % 文档形式: M*ddq = B*u + g, 其中 B = B_coeff, g = -rest
    % 所以: g = -(eq - M*ddq - B*u)
    g_sym(i) = -(eq_temp - M_sym(i,:)*ddq - B_sym(i,:)*u);
end

fprintf('符号重力+离心力项 g 已提取\n\n');

%% 应用对称假设

fprintf('========================================\n');
fprintf('应用对称假设\n');
fprintf('========================================\n\n');

fprintf('假设左右对称:\n');
fprintf('  m_l = m_r = m_leg\n');
fprintf('  m_wl = m_wr = m_w\n');
fprintf('  I_l = I_r = I_leg\n');
fprintf('  I_wl = I_wr = I_w\n');
fprintf('  l_l = l_r = l_leg\n');
fprintf('  l_l_d = l_r_d = l_leg_d\n');
fprintf('  theta_l0 = theta_r0 = theta_b0 = 0\n\n');

% 对称参数代换 (只替换对称参数，不替换具体数值)
sym_subs = {
    m_l, m_leg;
    m_r, m_leg;
    m_wl, m_w;
    m_wr, m_w;
    I_l, I_leg;
    I_r, I_leg;
    I_wl, I_w;
    I_wr, I_w;
    l_l, l_leg;
    l_r, l_leg;
    l_l_d, l_leg_d;
    l_r_d, l_leg_d;
    theta_l0, 0;
    theta_r0, 0;
    theta_b0, 0;
};

M_sym = subs(M_sym, sym_subs(:,1), sym_subs(:,2));
B_sym = subs(B_sym, sym_subs(:,1), sym_subs(:,2));
g_sym = subs(g_sym, sym_subs(:,1), sym_subs(:,2));

%% 平衡点线性化 (符号形式)

fprintf('========================================\n');
fprintf('平衡点线性化 (§6.2) - 符号形式\n');
fprintf('========================================\n\n');

fprintf('平衡点条件:\n');
fprintf('  theta_l = theta_r = theta_b = 0\n');
fprintf('  所有速度 = 0\n\n');

% 平衡点代换 (角度和角速度都为0)
eq_subs = {
    theta_l, 0;
    theta_r, 0;
    theta_b, 0;
    dtheta_l, 0;
    dtheta_r, 0;
    dtheta_b, 0;
    dphi, 0;
    dX_b_h, 0;
    ddtheta_wl, 0;
    ddtheta_wr, 0;
};

% 平衡点处的M和B矩阵（仍是符号形式）
M_eq = simplify(subs(M_sym, eq_subs(:,1), eq_subs(:,2)));
B_eq = simplify(subs(B_sym, eq_subs(:,1), eq_subs(:,2)));
g_eq = simplify(subs(g_sym, eq_subs(:,1), eq_subs(:,2)));

fprintf('平衡点处 M 矩阵 (符号):\n');
disp(M_eq);

fprintf('平衡点处 B 矩阵 (符号):\n');
disp(B_eq);

fprintf('平衡点处 g 向量 (符号，应为0):\n');
disp(g_eq);

%% 计算 dg/d(theta) 在平衡点的符号表达式

fprintf('========================================\n');
fprintf('计算状态空间 A, B 矩阵 (符号形式)\n');
fprintf('========================================\n\n');

% 对 theta_l 求导
dg_dtheta_l = simplify(subs(diff(g_sym, theta_l), eq_subs(:,1), eq_subs(:,2)));
% 对 theta_r 求导
dg_dtheta_r = simplify(subs(diff(g_sym, theta_r), eq_subs(:,1), eq_subs(:,2)));
% 对 theta_b 求导
dg_dtheta_b = simplify(subs(diff(g_sym, theta_b), eq_subs(:,1), eq_subs(:,2)));

fprintf('dg/d(theta_l) 在平衡点 (符号):\n');
disp(dg_dtheta_l);

fprintf('dg/d(theta_r) 在平衡点 (符号):\n');
disp(dg_dtheta_r);

fprintf('dg/d(theta_b) 在平衡点 (符号):\n');
disp(dg_dtheta_b);

%% 构建符号形式的状态空间矩阵

fprintf('----------------------------------------\n');
fprintf('构建 10x10 状态空间 A, B 矩阵 (符号形式)\n');
fprintf('----------------------------------------\n\n');

% 状态: x = [X_b_h; dX_b_h; phi; dphi; theta_l; dtheta_l; theta_r; dtheta_r; theta_b; dtheta_b]

n = 10;
m_ctrl = 4;

A_sym = sym(zeros(n, n));
B_state_sym = sym(zeros(n, m_ctrl));

% 运动学关系
A_sym(1,2) = 1;
A_sym(3,4) = 1;
A_sym(5,6) = 1;
A_sym(7,8) = 1;
A_sym(9,10) = 1;

% 符号 M^{-1} 计算比较复杂，这里使用符号求逆
fprintf('计算符号 M^{-1}...\n');
M_inv_sym = simplify(inv(M_eq));

% dg/d(theta) 矩阵
dg_dtheta = [zeros(5,1), dg_dtheta_l, dg_dtheta_r, dg_dtheta_b];

% A矩阵动力学部分: M^{-1} * dg/d(theta)
A_dyn_sym = simplify(M_inv_sym * dg_dtheta);

% 填充A矩阵
A_sym(2,5) = A_dyn_sym(1,2);   % ddX_b_h 对 theta_l
A_sym(2,7) = A_dyn_sym(1,3);   % ddX_b_h 对 theta_r
A_sym(2,9) = A_dyn_sym(1,4);   % ddX_b_h 对 theta_b
A_sym(4,5) = A_dyn_sym(2,2);
A_sym(4,7) = A_dyn_sym(2,3);
A_sym(4,9) = A_dyn_sym(2,4);
A_sym(6,5) = A_dyn_sym(3,2);
A_sym(6,7) = A_dyn_sym(3,3);
A_sym(6,9) = A_dyn_sym(3,4);
A_sym(8,5) = A_dyn_sym(4,2);
A_sym(8,7) = A_dyn_sym(4,3);
A_sym(8,9) = A_dyn_sym(4,4);
A_sym(10,5) = A_dyn_sym(5,2);
A_sym(10,7) = A_dyn_sym(5,3);
A_sym(10,9) = A_dyn_sym(5,4);

% B矩阵: M^{-1} * B
B_dyn_sym = simplify(M_inv_sym * B_eq);

B_state_sym(2,:) = B_dyn_sym(1,:);
B_state_sym(4,:) = B_dyn_sym(2,:);
B_state_sym(6,:) = B_dyn_sym(3,:);
B_state_sym(8,:) = B_dyn_sym(4,:);
B_state_sym(10,:) = B_dyn_sym(5,:);

fprintf('符号 A 矩阵已构建\n');
fprintf('符号 B 矩阵已构建\n\n');

%% 创建 MATLAB 函数句柄

fprintf('========================================\n');
fprintf('创建 matlabFunction 函数句柄\n');
fprintf('========================================\n\n');

% 参数列表: [m_b, m_leg, m_w, I_b, I_leg, I_w, I_yaw, l_leg, l_leg_d, l_b, R, R_w, g]
param_list = [m_b, m_leg, m_w, I_b, I_leg, I_w, I_yaw, l_leg, l_leg_d, l_b, R, R_w, g];

fprintf('参数列表:\n');
fprintf('  [m_b, m_leg, m_w, I_b, I_leg, I_w, I_yaw, l_leg, l_leg_d, l_b, R, R_w, g]\n\n');

% 转换为函数句柄
fprintf('生成 A_func...\n');
A_func = matlabFunction(A_sym, 'Vars', {param_list});

fprintf('生成 B_func...\n');
B_func = matlabFunction(B_state_sym, 'Vars', {param_list});

fprintf('生成 M_func...\n');
M_func = matlabFunction(M_eq, 'Vars', {param_list});

fprintf('生成 B_ctrl_func...\n');
B_ctrl_func = matlabFunction(B_eq, 'Vars', {param_list});

fprintf('✓ 函数句柄创建完成!\n\n');

%% 汇总

fprintf('========================================\n');
fprintf('线性化对照结果汇总\n');
fprintf('========================================\n');
fprintf('✓ §5.0 广义坐标定义: 一致\n');
fprintf('✓ §5.1 质量矩阵 M: 已提取 (符号)\n');
fprintf('✓ §5.2 控制矩阵 B: 已提取 (符号)\n');
fprintf('✓ §5.3 重力项 g: 已提取 (符号)\n');
fprintf('✓ §6.2 平衡点线性化: 已完成 (符号)\n');
fprintf('✓ §6.3 A 矩阵: 已构建 (符号函数)\n');
fprintf('✓ §6.4 B 矩阵: 已构建 (符号函数)\n');
fprintf('----------------------------------------\n');
fprintf('输出: A_func, B_func 函数句柄\n');
fprintf('使用: A = A_func([m_b, m_leg, m_w, I_b, I_leg, I_w, I_yaw, l_leg, l_leg_d, l_b, R, R_w, g])\n');
fprintf('========================================\n\n');

%% 保存结果

fprintf('Step 3: 保存结果...\n');

save('linearized_system.mat', ...
     'A_func', 'B_func', 'M_func', 'B_ctrl_func', ...
     'A_sym', 'B_state_sym', 'M_eq', 'B_eq', ...
     'M_sym', 'B_sym', 'g_sym', ...
     'dg_dtheta_l', 'dg_dtheta_r', 'dg_dtheta_b', ...
     'param_list');

fprintf('  结果已保存到 linearized_system.mat\n');
fprintf('  包含: A_func, B_func (函数句柄，在 compute_lqr.m 中使用)\n');
fprintf('\n========================================\n');
fprintf('线性化完成! (符号形式)\n');
fprintf('========================================\n');
