% linearize_system_v2.m
% 系统线性化 - 在平衡点附近线性化（符号形式）
%
% 输入: dynamics_new_coords.mat (用广义坐标表示的方程)
% 输出: linearized_system.mat (符号形式的 A, B 矩阵函数)
%
% 对照推导文档 §4.3, §5, §6
%
% ========== 状态向量 (10维, 对应推导 §4.3) ==========
%   X = [X_b^h, V_b^h, phi, dphi, theta_l, dtheta_l, theta_r, dtheta_r, theta_b, dtheta_b]^T
%   其中:
%     X_b^h    : 机体水平位置 (m)
%     V_b^h    : 机体水平速度 (m/s)
%     phi      : 偏航角 Yaw (rad)
%     dphi     : 偏航角速度 (rad/s)
%     theta_l  : 左腿与Z轴负方向夹角 (rad)
%     dtheta_l : 左腿角速度 (rad/s)
%     theta_r  : 右腿与Z轴负方向夹角 (rad)
%     dtheta_r : 右腿角速度 (rad/s)
%     theta_b  : 机体俯仰角 Pitch (rad)
%     dtheta_b : 机体俯仰角速度 (rad/s)
%
% ========== 控制向量 (4维, 对应推导 §4.3) ==========
%   U = [T_{r→b}, T_{l→b}, T_{wr→r}, T_{wl→l}]^T
%   其中:
%     T_{r→b}  : 右腿对机体的扭矩（右髋关节电机）(Nm)
%     T_{l→b}  : 左腿对机体的扭矩（左髋关节电机）(Nm)
%     T_{wr→r} : 右轮对右腿的扭矩（右轮电机）(Nm)
%     T_{wl→l} : 左轮对左腿的扭矩（左轮电机）(Nm)
%
% 注意：本脚本只输出符号形式，具体参数在 compute_lqr.m 中代入

clear; clc;
fprintf('========================================\n');
fprintf('系统线性化 (对照推导文档 §4.3, §5, §6)\n');
fprintf('========================================\n\n');

%% ========================================
%  Step 1: 定义符号变量
%  ========================================

fprintf('Step 1: 定义符号变量...\n\n');

% ========== 广义坐标及其导数 ==========
syms X_b_h dX_b_h ddX_b_h real     % 机体水平位置/速度/加速度
syms phi dphi ddphi real           % yaw角
syms theta_l dtheta_l ddtheta_l real
syms theta_r dtheta_r ddtheta_r real
syms theta_b dtheta_b ddtheta_b real

% ========== 轮角加速度 ==========
syms ddtheta_wl ddtheta_wr real

% ========== 控制力矩 (对应 §1.4 和 §4.3) ==========
% T_{r→b}: 右腿对机体的扭矩（右髋关节电机）
% T_{l→b}: 左腿对机体的扭矩（左髋关节电机）
% T_{wr→r}: 右轮对右腿的扭矩（右轮电机）
% T_{wl→l}: 左轮对左腿的扭矩（左轮电机）
syms T_r_to_b T_l_to_b T_wr_to_r T_wl_to_l real

% ========== 物理参数 (保持符号形式) ==========
syms m_b m_l m_r m_wl m_wr real        % 质量
syms I_b I_l I_r I_wl I_wr I_yaw real  % 转动惯量
syms l_b l_l l_r l_l_d l_r_d real      % 长度
syms theta_l0 theta_r0 theta_b0 real   % 初始角偏置
syms R R_w g real                       % 轮半径、半轮距、重力

%% ========================================
%  Step 2: 加载动力学方程
%  ========================================

fprintf('Step 2: 加载动力学方程...\n');
load('dynamics_new_coords.mat');
fprintf('  已加载 dynamics_new_coords.mat (推导 §3 代入运动学后的方程)\n\n');

%% ========================================
%  §5.0 将轮角加速度代换为广义坐标
%  ========================================

fprintf('========================================\n');
fprintf('§5.0 将轮角加速度代换为广义坐标\n');
fprintf('========================================\n\n');

% 腿角度相关项 (leg_term)
% leg_term = (l_r*cos(theta_r)*ddtheta_r + l_l*cos(theta_l)*ddtheta_l)/2
%          - (l_r*sin(theta_r)*dtheta_r^2 + l_l*sin(theta_l)*dtheta_l^2)/2
leg_term = (l_r*cos(theta_r)*ddtheta_r + l_l*cos(theta_l)*ddtheta_l)/2 ...
         - (l_r*sin(theta_r)*dtheta_r^2 + l_l*sin(theta_l)*dtheta_l^2)/2;

% 轮角加速度用广义坐标表示
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

%% ========================================
%  §5 构建 M*ddq = B*u + g 形式
%  ========================================

fprintf('========================================\n');
fprintf('§5 构建 M*ddq = B*u + g 形式\n');
fprintf('========================================\n\n');

% ========== 广义加速度向量 (5维) ==========
% ddq = [ddX_b_h, ddphi, ddtheta_l, ddtheta_r, ddtheta_b]^T
ddq = [ddX_b_h; ddphi; ddtheta_l; ddtheta_r; ddtheta_b];

% ========== 控制输入向量 (4维, 对应 §4.3) ==========
% U = [T_{r→b}, T_{l→b}, T_{wr→r}, T_{wl→l}]^T
u = [T_r_to_b; T_l_to_b; T_wr_to_r; T_wl_to_l];

fprintf('----------------------------------------\n');
fprintf('§4.3 状态向量和控制输入定义\n');
fprintf('----------------------------------------\n');
fprintf('状态向量 X (10维):\n');
fprintf('  X = [X_b^h, V_b^h, phi, dphi, theta_l, dtheta_l, theta_r, dtheta_r, theta_b, dtheta_b]^T\n\n');
fprintf('控制输入向量 U (4维):\n');
fprintf('  U = [T_{r→b}, T_{l→b}, T_{wr→r}, T_{wl→l}]^T\n');
fprintf('  其中:\n');
fprintf('    T_{r→b}  : 右腿对机体的扭矩（右髋关节电机）\n');
fprintf('    T_{l→b}  : 左腿对机体的扭矩（左髋关节电机）\n');
fprintf('    T_{wr→r} : 右轮对右腿的扭矩（右轮电机）\n');
fprintf('    T_{wl→l} : 左轮对左腿的扭矩（左轮电机）\n');
fprintf('✓ 与推导文档 §4.3 一致!\n\n');

%% ========================================
%  §5.1 提取质量矩阵 M (5×5)
%  ========================================

fprintf('----------------------------------------\n');
fprintf('§5.1 质量矩阵 M (5×5)\n');
fprintf('----------------------------------------\n');

M_sym = sym(zeros(5,5));
eqs = {eq1_new, eq2_new, eq3_new, eq4_new, eq5_new};

for i = 1:5
    for j = 1:5
        M_sym(i,j) = diff(eqs{i}, ddq(j));
    end
end

fprintf('符号质量矩阵 M 已提取\n\n');

%% ========================================
%  §5.2 提取控制矩阵 B (5×4)
%  ========================================

fprintf('----------------------------------------\n');
fprintf('§5.2 控制矩阵 B (5×4)\n');
fprintf('----------------------------------------\n');

% 方程形式: eq = M*ddq + B_raw*u + rest = 0
% 标准形式: M*ddq = B*u + g
% 所以: B = -B_raw (需要取负号!)

B_raw = sym(zeros(5,4));
for i = 1:5
    for j = 1:4
        B_raw(i,j) = diff(eqs{i}, u(j));
    end
end

% 转换为标准形式: B = -B_raw
B_sym = -B_raw;

fprintf('符号控制矩阵 B 已提取 (已转换为标准形式)\n\n');

%% ========================================
%  §5.3 提取重力+离心力项 g (5×1)
%  ========================================

fprintf('----------------------------------------\n');
fprintf('§5.3 重力+离心力项 g (5×1)\n');
fprintf('----------------------------------------\n');

g_sym = sym(zeros(5,1));
for i = 1:5
    % 方程形式: M*ddq + B_raw*u + rest = 0
    % 标准形式: M*ddq = B*u + g
    % 所以: g = -rest = -(eq - M*ddq - B_raw*u)
    % 注意这里用 B_raw 而不是 B_sym
    g_sym(i) = -(eqs{i} - M_sym(i,:)*ddq - B_raw(i,:)*u);
end

fprintf('符号重力+离心力项 g 已提取\n\n');

fprintf('========================================\n');
fprintf('使用左右腿参数不同的建模\n');
fprintf('========================================\n\n');
fprintf('左腿参数: m_l, I_l, l_l, l_l_d\n');
fprintf('右腿参数: m_r, I_r, l_r, l_r_d\n');
fprintf('左轮参数: m_wl, I_wl\n');
fprintf('右轮参数: m_wr, I_wr\n\n');

%% ========================================
%  §6 平衡点线性化
%  ========================================

fprintf('========================================\n');
fprintf('§6 平衡点线性化\n');
fprintf('========================================\n\n');

fprintf('平衡点条件:\n');
fprintf('  theta_l = theta_r = theta_b = 0\n');
fprintf('  所有速度 = 0\n\n');

% 平衡点代换
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

% 平衡点处的 M, B, g 矩阵
M_eq = simplify(subs(M_sym, eq_subs(:,1), eq_subs(:,2)));
B_eq = simplify(subs(B_sym, eq_subs(:,1), eq_subs(:,2)));
g_eq = simplify(subs(g_sym, eq_subs(:,1), eq_subs(:,2)));

fprintf('平衡点处 M 矩阵:\n');
disp(M_eq);

fprintf('平衡点处 B 矩阵:\n');
disp(B_eq);

fprintf('平衡点处 g 向量 (应为0):\n');
disp(g_eq);

%% ========================================
%  计算 dg/d(theta) 在平衡点
%  ========================================

fprintf('========================================\n');
fprintf('计算状态空间 A, B 矩阵\n');
fprintf('========================================\n\n');

% 对角度求导
dg_dtheta_l = simplify(subs(diff(g_sym, theta_l), eq_subs(:,1), eq_subs(:,2)));
dg_dtheta_r = simplify(subs(diff(g_sym, theta_r), eq_subs(:,1), eq_subs(:,2)));
dg_dtheta_b = simplify(subs(diff(g_sym, theta_b), eq_subs(:,1), eq_subs(:,2)));

fprintf('dg/d(theta_l) 在平衡点:\n');
disp(dg_dtheta_l);

fprintf('dg/d(theta_r) 在平衡点:\n');
disp(dg_dtheta_r);

fprintf('dg/d(theta_b) 在平衡点:\n');
disp(dg_dtheta_b);

%% ========================================
%  构建 10x10 状态空间矩阵
%  ========================================

fprintf('----------------------------------------\n');
fprintf('构建 10x10 状态空间 A, B 矩阵\n');
fprintf('----------------------------------------\n\n');

% 状态向量 (对应 §4.3):
% X = [X_b^h, V_b^h, phi, dphi, theta_l, dtheta_l, theta_r, dtheta_r, theta_b, dtheta_b]^T
%      [  1  ,   2  ,  3 ,  4  ,    5   ,    6    ,    7   ,    8    ,    9   ,   10    ]

n = 10;
m_ctrl = 4;

A_sym = sym(zeros(n, n));
B_state_sym = sym(zeros(n, m_ctrl));

% ========== 运动学关系 (位置-速度) ==========
A_sym(1,2) = 1;   % dX_b^h/dt = V_b^h
A_sym(3,4) = 1;   % dphi/dt = dphi
A_sym(5,6) = 1;   % dtheta_l/dt = dtheta_l
A_sym(7,8) = 1;   % dtheta_r/dt = dtheta_r
A_sym(9,10) = 1;  % dtheta_b/dt = dtheta_b

% ========== 计算 M^{-1} ==========
fprintf('计算符号 M^{-1}...\n');
M_inv_sym = simplify(inv(M_eq));

% ========== dg/d(theta) 矩阵 ==========
% 列顺序: [无, theta_l, theta_r, theta_b]
% 对应状态: [-, 5, 7, 9]
dg_dtheta = [zeros(5,1), dg_dtheta_l, dg_dtheta_r, dg_dtheta_b];

% ========== A矩阵动力学部分: M^{-1} * dg/d(theta) ==========
A_dyn_sym = simplify(M_inv_sym * dg_dtheta);

% 填充A矩阵 (加速度行)
% 行2: ddX_b^h 对应方程1
% 行4: ddphi 对应方程2
% 行6: ddtheta_l 对应方程3
% 行8: ddtheta_r 对应方程4
% 行10: ddtheta_b 对应方程5
%
% 列5: theta_l, 列7: theta_r, 列9: theta_b

A_sym(2,5) = A_dyn_sym(1,2);   % ddX_b_h 对 theta_l
A_sym(2,7) = A_dyn_sym(1,3);   % ddX_b_h 对 theta_r
A_sym(2,9) = A_dyn_sym(1,4);   % ddX_b_h 对 theta_b

A_sym(4,5) = A_dyn_sym(2,2);   % ddphi 对 theta_l
A_sym(4,7) = A_dyn_sym(2,3);   % ddphi 对 theta_r
A_sym(4,9) = A_dyn_sym(2,4);   % ddphi 对 theta_b

A_sym(6,5) = A_dyn_sym(3,2);   % ddtheta_l 对 theta_l
A_sym(6,7) = A_dyn_sym(3,3);   % ddtheta_l 对 theta_r
A_sym(6,9) = A_dyn_sym(3,4);   % ddtheta_l 对 theta_b

A_sym(8,5) = A_dyn_sym(4,2);   % ddtheta_r 对 theta_l
A_sym(8,7) = A_dyn_sym(4,3);   % ddtheta_r 对 theta_r
A_sym(8,9) = A_dyn_sym(4,4);   % ddtheta_r 对 theta_b

A_sym(10,5) = A_dyn_sym(5,2);  % ddtheta_b 对 theta_l
A_sym(10,7) = A_dyn_sym(5,3);  % ddtheta_b 对 theta_r
A_sym(10,9) = A_dyn_sym(5,4);  % ddtheta_b 对 theta_b

% ========== B矩阵: M^{-1} * B ==========
B_dyn_sym = simplify(M_inv_sym * B_eq);

B_state_sym(2,:) = B_dyn_sym(1,:);   % ddX_b_h
B_state_sym(4,:) = B_dyn_sym(2,:);   % ddphi
B_state_sym(6,:) = B_dyn_sym(3,:);   % ddtheta_l
B_state_sym(8,:) = B_dyn_sym(4,:);   % ddtheta_r
B_state_sym(10,:) = B_dyn_sym(5,:);  % ddtheta_b

fprintf('符号 A 矩阵已构建 (10×10)\n');
fprintf('符号 B 矩阵已构建 (10×4)\n\n');

%% ========================================
%  创建 MATLAB 函数句柄
%  ========================================

fprintf('========================================\n');
fprintf('创建 matlabFunction 函数句柄\n');
fprintf('========================================\n\n');

% 参数列表 (左右腿参数分开)
param_list = [m_b, m_l, m_r, m_wl, m_wr, I_b, I_l, I_r, I_wl, I_wr, I_yaw, l_l, l_r, l_l_d, l_r_d, l_b, R, R_w, g, theta_l0, theta_r0, theta_b0];

fprintf('参数列表 (左右腿参数分开):\n');
fprintf('  [m_b, m_l, m_r, m_wl, m_wr, I_b, I_l, I_r, I_wl, I_wr, I_yaw, l_l, l_r, l_l_d, l_r_d, l_b, R, R_w, g, theta_l0, theta_r0, theta_b0]\n\n');

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

%% ========================================
%  汇总
%  ========================================

fprintf('========================================\n');
fprintf('线性化结果汇总\n');
fprintf('========================================\n');
fprintf('状态向量 X (10维, §4.3):\n');
fprintf('  [X_b^h, V_b^h, phi, dphi, theta_l, dtheta_l, theta_r, dtheta_r, theta_b, dtheta_b]^T\n\n');
fprintf('控制向量 U (4维, §4.3):\n');
fprintf('  [T_{r→b}, T_{l→b}, T_{wr→r}, T_{wl→l}]^T\n\n');
fprintf('----------------------------------------\n');
fprintf('✓ §5.0 广义坐标定义: 一致\n');
fprintf('✓ §5.1 质量矩阵 M: 已提取\n');
fprintf('✓ §5.2 控制矩阵 B: 已提取\n');
fprintf('✓ §5.3 重力项 g: 已提取\n');
fprintf('✓ §6 平衡点线性化: 已完成\n');
fprintf('----------------------------------------\n');
fprintf('输出: A_func, B_func 函数句柄\n');
fprintf('========================================\n\n');

%% ========================================
%  保存结果
%  ========================================

fprintf('Step 3: 保存结果...\n');

save('linearized_system.mat', ...
     'A_func', 'B_func', 'M_func', 'B_ctrl_func', ...
     'A_sym', 'B_state_sym', 'M_eq', 'B_eq', ...
     'M_sym', 'B_sym', 'g_sym', ...
     'dg_dtheta_l', 'dg_dtheta_r', 'dg_dtheta_b', ...
     'param_list', ...
     'u', 'ddq');

fprintf('  结果已保存到 linearized_system.mat\n');
fprintf('\n========================================\n');
fprintf('线性化完成!\n');
fprintf('========================================\n');
