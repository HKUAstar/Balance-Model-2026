% compute_lqr_v2.m
% 轮腿机器人 LQR 控制器计算 - 整合版
%
% 功能: 加载动力学方程 → 代入数值参数 → 求解平衡点 → 线性化 → 计算 LQR
%
% 输入: dynamics_new_coords.mat (用广义坐标表示的方程)
% 输出: lqr_results.mat, lqr_fitting_results.mat
%
% ========== 状态向量 (10维, 对应推导 §4.3) ==========
%   X = [X_b^h, V_b^h, phi, dphi, theta_l, dtheta_l, theta_r, dtheta_r, theta_b, dtheta_b]^T
%
% ========== 控制向量 (4维, 对应推导 §4.3) ==========
%   U = [T_{r→b}, T_{l→b}, T_{wr→r}, T_{wl→l}]^T

clear; clc;
tic;
fprintf('========================================\n');
fprintf('轮腿机器人 LQR 控制器计算 (整合版)\n');
fprintf('========================================\n\n');

%% ========================================
%  Part 1: 定义物理参数 (先于符号计算)
%  ========================================

fprintf('========================================\n');
fprintf('Part 1: 定义物理参数\n');
fprintf('========================================\n\n');

% ==================== 物理常数 ====================
g_val = -9.81;              % 重力加速度 (m/s^2)

% ==================== 几何参数 ====================
R_val = 0.055;              % 轮子半径 (m)
R_w_val = (0.455+0.445)/4;  % 轮距/2 (m)

% ==================== 机体参数 ====================
m_b_val = 7.846;            % 机体质量 (kg)
I_b_val = 0.150;            % 机体俯仰转动惯量 (kg·m²)
l_b_val = 0.055;            % 机体质心到俯仰轴距离 (m)
I_yaw_val = 0.465;          % 整体yaw轴转动惯量 (kg·m²)
theta_b0_val = 0;           % 质心偏移角度 (rad)

% ==================== 轮子参数 ====================
m_wl_val = 0.19;            % 左轮质量 (kg)
m_wr_val = 0.19;            % 右轮质量 (kg)
I_wl_val = 0.000207897;     % 左轮转动惯量 (kg·m²)
I_wr_val = 0.000207897;     % 右轮转动惯量 (kg·m²)

% ==================== 腿部参数 (默认腿长 0.20m) ====================
l_l_val = 0.20;             % 左腿长度 (m)
l_r_val = 0.20;             % 右腿长度 (m)
m_l_val = 1.62;             % 左腿质量 (kg)
m_r_val = 1.62;             % 右腿质量 (kg)
I_l_val = 0.0339;           % 左腿转动惯量 (kg·m²)
I_r_val = 0.0339;           % 右腿转动惯量 (kg·m²)
l_l_d_val = 0.084685385397954;  % 左腿质心到轮轴距离 (m)
l_r_d_val = 0.084685385397954;  % 右腿质心到轮轴距离 (m)
theta_l0_val = 0.8152113973036387;  % 左腿偏移角度 (rad)
theta_r0_val = 0.8152113973036387;  % 右腿偏移角度 (rad)

fprintf('✓ 物理参数设置完成\n\n');

%% ========================================
%  Part 2: 定义符号变量并加载动力学
%  ========================================

fprintf('========================================\n');
fprintf('Part 2: 加载动力学方程\n');
fprintf('========================================\n\n');

% ========== 广义坐标及其导数 ==========
syms X_b_h dX_b_h ddX_b_h real
syms phi dphi ddphi real
syms theta_l dtheta_l ddtheta_l real
syms theta_r dtheta_r ddtheta_r real
syms theta_b dtheta_b ddtheta_b real

% ========== 轮角加速度 ==========
syms ddtheta_wl ddtheta_wr real

% ========== 控制力矩 ==========
syms T_r_to_b T_l_to_b T_wr_to_r T_wl_to_l real

% ========== 物理参数符号 ==========
syms m_b m_l m_r m_wl m_wr real
syms I_b I_l I_r I_wl I_wr I_yaw real
syms l_b l_l l_r l_l_d l_r_d real
syms theta_l0 theta_r0 theta_b0 real
syms R R_w g real

% 加载动力学方程
load('dynamics_new_coords.mat');
fprintf('✓ 已加载 dynamics_new_coords.mat\n\n');

%% ========================================
%  Part 3: 符号处理 - 轮角加速度代换
%  ========================================

fprintf('========================================\n');
fprintf('Part 3: 轮角加速度代换\n');
fprintf('========================================\n\n');

leg_term = (l_r*cos(theta_r)*ddtheta_r + l_l*cos(theta_l)*ddtheta_l)/2 ...
         - (l_r*sin(theta_r)*dtheta_r^2 + l_l*sin(theta_l)*dtheta_l^2)/2;

ddtheta_wr_sub = (ddX_b_h + R_w*ddphi)/R - leg_term/R;
ddtheta_wl_sub = (ddX_b_h - R_w*ddphi)/R - leg_term/R;

wheel_subs = {ddtheta_wr, ddtheta_wr_sub; ddtheta_wl, ddtheta_wl_sub};

eq1_new = simplify(subs(eq1_new, wheel_subs(:,1), wheel_subs(:,2)));
eq2_new = simplify(subs(eq2_new, wheel_subs(:,1), wheel_subs(:,2)));
eq3_new = simplify(subs(eq3_new, wheel_subs(:,1), wheel_subs(:,2)));
eq4_new = simplify(subs(eq4_new, wheel_subs(:,1), wheel_subs(:,2)));
eq5_new = simplify(subs(eq5_new, wheel_subs(:,1), wheel_subs(:,2)));

fprintf('✓ 轮角加速度已代换为广义坐标\n\n');

%% ========================================
%  Part 4: 提取 M, B, g 矩阵 (符号形式)
%  ========================================

fprintf('========================================\n');
fprintf('Part 4: 提取 M, B, g 矩阵\n');
fprintf('========================================\n\n');

ddq = [ddX_b_h; ddphi; ddtheta_l; ddtheta_r; ddtheta_b];
u = [T_r_to_b; T_l_to_b; T_wr_to_r; T_wl_to_l];
eqs = {eq1_new, eq2_new, eq3_new, eq4_new, eq5_new};

% 提取 M 矩阵
M_sym = sym(zeros(5,5));
for i = 1:5
    for j = 1:5
        M_sym(i,j) = diff(eqs{i}, ddq(j));
    end
end

% 提取 B 矩阵
B_raw = sym(zeros(5,4));
for i = 1:5
    for j = 1:4
        B_raw(i,j) = diff(eqs{i}, u(j));
    end
end
B_sym = -B_raw;

% 提取 g 向量
g_sym = sym(zeros(5,1));
for i = 1:5
    g_sym(i) = -(eqs{i} - M_sym(i,:)*ddq - B_raw(i,:)*u);
end

fprintf('✓ M, B, g 符号矩阵已提取\n\n');

%% ========================================
%  Part 5: 代入物理参数数值
%  ========================================

fprintf('========================================\n');
fprintf('Part 5: 代入物理参数数值\n');
fprintf('========================================\n\n');

% 物理参数代换表
param_subs = {
    m_b, m_b_val;
    m_l, m_l_val;
    m_r, m_r_val;
    m_wl, m_wl_val;
    m_wr, m_wr_val;
    I_b, I_b_val;
    I_l, I_l_val;
    I_r, I_r_val;
    I_wl, I_wl_val;
    I_wr, I_wr_val;
    I_yaw, I_yaw_val;
    l_l, l_l_val;
    l_r, l_r_val;
    l_l_d, l_l_d_val;
    l_r_d, l_r_d_val;
    l_b, l_b_val;
    R, R_val;
    R_w, R_w_val;
    g, g_val;
    theta_l0, theta_l0_val;
    theta_r0, theta_r0_val;
    theta_b0, theta_b0_val;
};

% 代入物理参数
M_param = simplify(subs(M_sym, param_subs(:,1), param_subs(:,2)));
B_param = simplify(subs(B_sym, param_subs(:,1), param_subs(:,2)));
g_param = simplify(subs(g_sym, param_subs(:,1), param_subs(:,2)));

fprintf('✓ 物理参数已代入 (M, B, g 现在只含状态变量)\n\n');

%% ========================================
%  Part 6: 求解平衡点 theta_l*, theta_r*, theta_b*
%  ========================================

fprintf('========================================\n');
fprintf('Part 6: 求解平衡点\n');
fprintf('========================================\n\n');

fprintf('平衡点条件:\n');
fprintf('  所有速度 = 0, 加速度 = 0, 控制输入 = 0\n');
fprintf('  求解 theta_l*, theta_r*, theta_b* 使得 g = 0\n\n');

% 代入平衡点条件 (速度=0, 加速度=0, 控制=0, phi=0, X=0)，保留 theta_l, theta_r, theta_b
eq_subs_partial = {
    dtheta_l, 0;
    dtheta_r, 0;
    dtheta_b, 0;
    ddtheta_l, 0;
    ddtheta_r, 0;
    ddtheta_b, 0;
    ddX_b_h, 0;
    ddphi, 0;
    phi, 0;
    dphi, 0;
    X_b_h, 0;
    dX_b_h, 0;
    T_r_to_b, 0;
    T_l_to_b, 0;
    T_wr_to_r, 0;
    T_wl_to_l, 0;
};

g_at_eq = simplify(subs(g_param, eq_subs_partial(:,1), eq_subs_partial(:,2)));

fprintf('g 向量 (含 theta_l, theta_r, theta_b):\n');
disp(g_at_eq);

% 符号求解 g = 0
fprintf('正在符号求解 g = 0 ...\n');

% g 向量中只有 g(2), g(3), g(4) 非零，且各自只含一个 theta
% g(2) 只含 theta_b → 解 g(2)=0 得 theta_b
% g(3) 只含 theta_r → 解 g(3)=0 得 theta_r
% g(4) 只含 theta_l → 解 g(4)=0 得 theta_l

% 解 theta_b (从 g(2)=0)
if g_at_eq(2) == 0
    theta_b_star = 0;
    fprintf('  g(2) = 0 恒成立, theta_b* = 0\n');
else
    sol_b = solve(g_at_eq(2) == 0, theta_b, 'Real', true);
    % 选择接近0的解
    sol_b_vals = double(sol_b);
    [~, idx_b] = min(abs(sol_b_vals));
    theta_b_star = sol_b_vals(idx_b);
    fprintf('  解 g(2)=0: theta_b* = %.10f rad (%.6f deg)\n', theta_b_star, rad2deg(theta_b_star));
end

% 解 theta_r (从 g(3)=0)
if g_at_eq(3) == 0
    theta_r_star = 0;
    fprintf('  g(3) = 0 恒成立, theta_r* = 0\n');
else
    sol_r = solve(g_at_eq(3) == 0, theta_r, 'Real', true);
    sol_r_vals = double(sol_r);
    [~, idx_r] = min(abs(sol_r_vals));
    theta_r_star = sol_r_vals(idx_r);
    fprintf('  解 g(3)=0: theta_r* = %.10f rad (%.6f deg)\n', theta_r_star, rad2deg(theta_r_star));
end

% 解 theta_l (从 g(4)=0)
if g_at_eq(4) == 0
    theta_l_star = 0;
    fprintf('  g(4) = 0 恒成立, theta_l* = 0\n');
else
    sol_l = solve(g_at_eq(4) == 0, theta_l, 'Real', true);
    sol_l_vals = double(sol_l);
    [~, idx_l] = min(abs(sol_l_vals));
    theta_l_star = sol_l_vals(idx_l);
    fprintf('  解 g(4)=0: theta_l* = %.10f rad (%.6f deg)\n', theta_l_star, rad2deg(theta_l_star));
end

% 验证
theta_star = [theta_l_star; theta_r_star; theta_b_star];
g_func = matlabFunction(g_at_eq, 'Vars', {[theta_l; theta_r; theta_b]});
fval = g_func(theta_star);

fprintf('\n平衡点求解结果:\n');
fprintf('  theta_l* = %.10f rad (%.6f deg)\n', theta_l_star, rad2deg(theta_l_star));
fprintf('  theta_r* = %.10f rad (%.6f deg)\n', theta_r_star, rad2deg(theta_r_star));
fprintf('  theta_b* = %.10f rad (%.6f deg)\n', theta_b_star, rad2deg(theta_b_star));
fprintf('  残差 ||g||: %.2e\n', norm(fval));
fprintf('  ✓ 符号求解成功!\n\n');

%% ========================================
%  Part 7: 在平衡点处线性化
%  ========================================

fprintf('========================================\n');
fprintf('Part 7: 在平衡点处线性化\n');
fprintf('========================================\n\n');

% 完整的平衡点代换
eq_subs_full = {
    theta_l, theta_l_star;
    theta_r, theta_r_star;
    theta_b, theta_b_star;
    dtheta_l, 0;
    dtheta_r, 0;
    dtheta_b, 0;
    phi, 0;
    dphi, 0;
    X_b_h, 0;
    dX_b_h, 0;
};

% 平衡点处的 M, B, g 矩阵 (数值)
M_eq = double(subs(M_param, eq_subs_full(:,1), eq_subs_full(:,2)));
B_eq = double(subs(B_param, eq_subs_full(:,1), eq_subs_full(:,2)));
g_eq = double(subs(g_param, eq_subs_full(:,1), eq_subs_full(:,2)));

fprintf('平衡点处 M 矩阵:\n');
disp(M_eq);

fprintf('平衡点处 B 矩阵:\n');
disp(B_eq);

fprintf('平衡点处 g 向量 (验证: 应为0):\n');
disp(g_eq);

% 计算 dg/d(theta) 在平衡点
dg_dtheta_l = double(subs(diff(g_param, theta_l), eq_subs_full(:,1), eq_subs_full(:,2)));
dg_dtheta_r = double(subs(diff(g_param, theta_r), eq_subs_full(:,1), eq_subs_full(:,2)));
dg_dtheta_b = double(subs(diff(g_param, theta_b), eq_subs_full(:,1), eq_subs_full(:,2)));

fprintf('dg/d(theta_l) 在平衡点:\n');
disp(dg_dtheta_l);

fprintf('dg/d(theta_r) 在平衡点:\n');
disp(dg_dtheta_r);

fprintf('dg/d(theta_b) 在平衡点:\n');
disp(dg_dtheta_b);

%% ========================================
%  Part 8: 构建 10x10 状态空间矩阵
%  ========================================

fprintf('========================================\n');
fprintf('Part 8: 构建状态空间 A, B 矩阵\n');
fprintf('========================================\n\n');

n = 10;
m_ctrl = 4;

A_num = zeros(n, n);
B_num = zeros(n, m_ctrl);

% 运动学关系 (位置-速度)
A_num(1,2) = 1;   % dX_b^h/dt = V_b^h
A_num(3,4) = 1;   % dphi/dt = dphi
A_num(5,6) = 1;   % dtheta_l/dt = dtheta_l
A_num(7,8) = 1;   % dtheta_r/dt = dtheta_r
A_num(9,10) = 1;  % dtheta_b/dt = dtheta_b

% 计算 M^{-1}
M_inv = inv(M_eq);

% dg/d(theta) 矩阵
dg_dtheta = [zeros(5,1), dg_dtheta_l, dg_dtheta_r, dg_dtheta_b];

% A矩阵动力学部分: M^{-1} * dg/d(theta)
A_dyn = M_inv * dg_dtheta;

% 填充A矩阵 (加速度行)
A_num(2,5) = A_dyn(1,2);   A_num(2,7) = A_dyn(1,3);   A_num(2,9) = A_dyn(1,4);
A_num(4,5) = A_dyn(2,2);   A_num(4,7) = A_dyn(2,3);   A_num(4,9) = A_dyn(2,4);
A_num(6,5) = A_dyn(3,2);   A_num(6,7) = A_dyn(3,3);   A_num(6,9) = A_dyn(3,4);
A_num(8,5) = A_dyn(4,2);   A_num(8,7) = A_dyn(4,3);   A_num(8,9) = A_dyn(4,4);
A_num(10,5) = A_dyn(5,2);  A_num(10,7) = A_dyn(5,3);  A_num(10,9) = A_dyn(5,4);

% B矩阵: M^{-1} * B
B_dyn = M_inv * B_eq;

B_num(2,:) = B_dyn(1,:);
B_num(4,:) = B_dyn(2,:);
B_num(6,:) = B_dyn(3,:);
B_num(8,:) = B_dyn(4,:);
B_num(10,:) = B_dyn(5,:);

fprintf('数值 A 矩阵 (10×10):\n');
disp(A_num);

fprintf('数值 B 矩阵 (10×4):\n');
disp(B_num);

%% ========================================
%  Part 9: 检查可控性
%  ========================================

fprintf('========================================\n');
fprintf('Part 9: 检查系统可控性\n');
fprintf('========================================\n\n');

Co = ctrb(A_num, B_num);
rank_Co = rank(Co);
fprintf('可控性矩阵秩: %d (系统维度: 10)\n', rank_Co);

if rank_Co < 10
    fprintf('  ⚠ 系统不完全可控 (秩=%d < 10)\n', rank_Co);
    fprintf('  物理原因: X_b^h 和 phi 是积分器状态\n');
    fprintf('  解决方案: 这是正常的! LQR仍可计算\n\n');
else
    fprintf('  ✓ 系统完全可控\n\n');
end

%% ========================================
%  Part 10: 设置 LQR 权重
%  ========================================

fprintf('========================================\n');
fprintf('Part 10: 设置 LQR 权重\n');
fprintf('========================================\n\n');

% Q矩阵: 状态权重
lqr_Q = diag([100, 1, 4000, 1, 1000, 10, 1000, 10, 40000, 1]);

% R矩阵: 控制输入权重
lqr_R = diag([1, 1, 10, 10]);

fprintf('Q矩阵 (状态权重):\n');
disp(lqr_Q);

fprintf('R矩阵 (控制权重):\n');
disp(lqr_R);

%% ========================================
%  Part 11: 计算 LQR 增益
%  ========================================

fprintf('========================================\n');
fprintf('Part 11: 计算 LQR 增益\n');
fprintf('========================================\n\n');

try
    [K, S, e] = lqr(A_num, B_num, lqr_Q, lqr_R);
    
    fprintf('✓ LQR增益矩阵 K (4×10):\n');
    disp(K);
    
    fprintf('闭环特征值:\n');
    disp(e);
    
    if all(real(e) < 1e-6)
        fprintf('✓ 闭环系统稳定!\n\n');
    else
        warning('闭环系统不稳定!');
    end
catch ME
    fprintf('✗ LQR计算失败: %s\n', ME.message);
    K = [];
end

%% ========================================
%  Part 12: 格式化输出
%  ========================================

fprintf('========================================\n');
fprintf('Part 12: 格式化输出 (C代码)\n');
fprintf('========================================\n\n');

if ~isempty(K)
    fprintf('// LQR增益矩阵 K[4][10]\n');
    fprintf('// 控制律: u = -K * X\n');
    fprintf('// 平衡点: theta_l*=%.6f, theta_r*=%.6f, theta_b*=%.6f (rad)\n\n', ...
            theta_l_star, theta_r_star, theta_b_star);
    
    fprintf('float K[4][10] = {\n');
    control_names = {'T_r_to_b', 'T_l_to_b', 'T_wr_to_r', 'T_wl_to_l'};
    for i = 1:4
        fprintf('    {');
        for j = 1:10
            if j < 10
                fprintf('%11.6ff, ', K(i,j));
            else
                fprintf('%11.6ff', K(i,j));
            end
        end
        if i < 4
            fprintf('},  // %s\n', control_names{i});
        else
            fprintf('}   // %s\n', control_names{i});
        end
    end
    fprintf('};\n\n');
    
    % 输出平衡点偏移量
    fprintf('// ═══════════════════════════════════════════════════════════════════════\n');
    fprintf('// 平衡点偏移量 (Equilibrium Point Offsets)\n');
    fprintf('// ═══════════════════════════════════════════════════════════════════════\n');
    fprintf('// 控制律: u = -K * (X - X_eq)\n');
    fprintf('// 其中 X_eq 是平衡点状态向量\n');
    fprintf('// 状态向量中只有 theta_l, theta_r, theta_b 有非零平衡点\n');
    fprintf('// 对应索引: theta_l -> X[4], theta_r -> X[6], theta_b -> X[8]\n');
    fprintf('// ═══════════════════════════════════════════════════════════════════════\n\n');
    
    fprintf('// 平衡点角度 (rad)\n');
    fprintf('float theta_l_eq = %.10ff;  // 左腿平衡点角度\n', theta_l_star);
    fprintf('float theta_r_eq = %.10ff;  // 右腿平衡点角度\n', theta_r_star);
    fprintf('float theta_b_eq = %.10ff;  // 机体俯仰平衡点角度\n\n', theta_b_star);
    
    fprintf('// 平衡点角度 (deg) - 仅供参考\n');
    fprintf('// theta_l_eq = %.6f deg\n', rad2deg(theta_l_star));
    fprintf('// theta_r_eq = %.6f deg\n', rad2deg(theta_r_star));
    fprintf('// theta_b_eq = %.6f deg\n\n', rad2deg(theta_b_star));
    
    fprintf('// 完整平衡点状态向量 X_eq[10]\n');
    fprintf('float X_eq[10] = {0.0f, 0.0f, 0.0f, 0.0f, %.10ff, 0.0f, %.10ff, 0.0f, %.10ff, 0.0f};\n', ...
            theta_l_star, theta_r_star, theta_b_star);
    fprintf('//                X_b^h  V_b^h  phi   dphi   theta_l  dtheta_l  theta_r  dtheta_r  theta_b  dtheta_b\n\n');
end

%% ========================================
%  Part 13: 腿长拟合功能
%  ========================================

fprintf('========================================\n');
fprintf('Part 13: 腿长拟合功能\n');
fprintf('========================================\n\n');

% 腿长参数查找表
Leg_data = [
    0.13, 0.13680373277071062, 0.06482330830804611, 0.01298843908;
    0.14, 0.14177326863693313, 0.06550465403312959, 0.01309498073;
    0.15, 0.14696847451069225, 0.06629277864141764, 0.01320078273;
    0.16, 0.15233205965915383, 0.06716142047336403, 0.01330738064;
    0.17, 0.15785793169809365, 0.06810379284592012, 0.013415874730;
    0.18, 0.16353393531619057, 0.06910968094268703, 0.01352706998;
    0.19, 0.16932790703247946, 0.07018361703417686, 0.0136415662890;
    0.20, 0.17523986903670066, 0.07130926798109767, 0.01375981819;
    0.21, 0.18125438063671734, 0.07249103737704407, 0.0138821752200;
    0.22, 0.18735804573062775, 0.07371456640312009, 0.01400890984;
    0.23, 0.1935440737919919,  0.07498472177717272, 0.01414023708;
    0.24, 0.19981164355462372, 0.07629346564418214, 0.01427632855;
    0.25, 0.20613804355334314, 0.07763950669601141, 0.01441732272;
    0.26, 0.2125254867068889,  0.07902583438344704, 0.0145633324100;
    0.27, 0.21897323717751446, 0.08044798692322885, 0.0147144504;
    0.28, 0.22546613138118995, 0.08189613177678175, 0.01487075369;
    0.29, 0.23200655335571885, 0.08337530089900726, 0.0150323068;
    0.30, 0.23859429875837357, 0.08488368158839482, 0.01519916419;
    0.31, 0.2452153170175142,  0.08641615416112892, 0.01537137233;
    0.32, 0.2518818852160671,  0.0879709275840604,  0.01554897121;
    0.33, 0.2585709469371994,  0.08955408756723503, 0.01573199564;
    0.34, 0.2652931949749182,  0.09115963635293857, 0.01592047625;
    0.35, 0.2720462993683244,  0.09278032657842933, 0.01611444043;
    0.36, 0.27883182906547804, 0.09442027801272353, 0.0163139131;
    0.37, 0.2856372146622355,  0.09608027060744573, 0.01651891744;
    0.38, 0.29246232389831,    0.0977599657323999,  0.0167294757;
    0.39, 0.29930908322334626, 0.09946520647945191, 0.01694561016;
    0.40, 0.30618331143287353, 0.10117420718740523, 0.01716734433;
];

enable_fitting = true;

if enable_fitting
    fprintf('正在计算不同腿长下的K矩阵和平衡点偏移量...\n\n');
    
    num_legs = size(Leg_data, 1);
    sample_size_2d = num_legs^2;
    K_sample_2d = zeros(sample_size_2d, 44);  % [l_l, l_r, K(40个)]
    offset_sample_2d = zeros(sample_size_2d, 5);  % [l_l, l_r, theta_l*, theta_r*, theta_b*]
    
    tic_fit = tic;
    idx = 0;
    
    for i = 1:num_legs
        for j = 1:num_legs
            idx = idx + 1;
            
            % 左腿参数
            l_l_fit = Leg_data(i, 1);
            theta_l0_fit = Leg_data(i, 2);
            l_l_d_fit = Leg_data(i, 3);
            I_l_fit = Leg_data(i, 4);
            
            % 右腿参数
            l_r_fit = Leg_data(j, 1);
            theta_r0_fit = Leg_data(j, 2);
            l_r_d_fit = Leg_data(j, 3);
            I_r_fit = Leg_data(j, 4);
            
            % 构建参数代换表
            param_subs_fit = {
                m_b, m_b_val; m_l, m_l_val; m_r, m_r_val;
                m_wl, m_wl_val; m_wr, m_wr_val;
                I_b, I_b_val; I_l, I_l_fit; I_r, I_r_fit;
                I_wl, I_wl_val; I_wr, I_wr_val; I_yaw, I_yaw_val;
                l_l, l_l_fit; l_r, l_r_fit;
                l_l_d, l_l_d_fit; l_r_d, l_r_d_fit;
                l_b, l_b_val; R, R_val; R_w, R_w_val; g, g_val;
                theta_l0, theta_l0_fit; theta_r0, theta_r0_fit; theta_b0, theta_b0_val;
            };
            
            % 代入物理参数
            M_fit = subs(M_sym, param_subs_fit(:,1), param_subs_fit(:,2));
            B_fit = subs(B_sym, param_subs_fit(:,1), param_subs_fit(:,2));
            g_fit = subs(g_sym, param_subs_fit(:,1), param_subs_fit(:,2));
            
            % 求解平衡点 (符号求解，每个方程独立)
            g_fit_eq = subs(g_fit, eq_subs_partial(:,1), eq_subs_partial(:,2));
            
            try
                % 符号求解各 theta
                % theta_b 从 g(2)=0
                if g_fit_eq(2) == 0
                    theta_b_fit = 0;
                else
                    sol_b = solve(g_fit_eq(2) == 0, theta_b, 'Real', true);
                    sol_b_vals = double(sol_b);
                    [~, idx_b] = min(abs(sol_b_vals));
                    theta_b_fit = sol_b_vals(idx_b);
                end
                
                % theta_r 从 g(3)=0
                if g_fit_eq(3) == 0
                    theta_r_fit_val = 0;
                else
                    sol_r = solve(g_fit_eq(3) == 0, theta_r, 'Real', true);
                    sol_r_vals = double(sol_r);
                    [~, idx_r] = min(abs(sol_r_vals));
                    theta_r_fit_val = sol_r_vals(idx_r);
                end
                
                % theta_l 从 g(4)=0
                if g_fit_eq(4) == 0
                    theta_l_fit_val = 0;
                else
                    sol_l = solve(g_fit_eq(4) == 0, theta_l, 'Real', true);
                    sol_l_vals = double(sol_l);
                    [~, idx_l] = min(abs(sol_l_vals));
                    theta_l_fit_val = sol_l_vals(idx_l);
                end
                
                theta_star_fit = [theta_l_fit_val; theta_r_fit_val; theta_b_fit];
                
                % 完整平衡点代换 (包括加速度和控制输入)
                eq_subs_fit = {
                    theta_l, theta_star_fit(1); theta_r, theta_star_fit(2); theta_b, theta_star_fit(3);
                    dtheta_l, 0; dtheta_r, 0; dtheta_b, 0;
                    ddtheta_l, 0; ddtheta_r, 0; ddtheta_b, 0;
                    ddX_b_h, 0; ddphi, 0;
                    phi, 0; dphi, 0; X_b_h, 0; dX_b_h, 0;
                    T_r_to_b, 0; T_l_to_b, 0; T_wr_to_r, 0; T_wl_to_l, 0;
                };
                
                M_eq_fit = double(subs(M_fit, eq_subs_fit(:,1), eq_subs_fit(:,2)));
                B_eq_fit = double(subs(B_fit, eq_subs_fit(:,1), eq_subs_fit(:,2)));
                
                dg_dtheta_l_fit = double(subs(diff(g_fit, theta_l), eq_subs_fit(:,1), eq_subs_fit(:,2)));
                dg_dtheta_r_fit = double(subs(diff(g_fit, theta_r), eq_subs_fit(:,1), eq_subs_fit(:,2)));
                dg_dtheta_b_fit = double(subs(diff(g_fit, theta_b), eq_subs_fit(:,1), eq_subs_fit(:,2)));
                
                M_inv_fit = inv(M_eq_fit);
                dg_dtheta_fit = [zeros(5,1), dg_dtheta_l_fit, dg_dtheta_r_fit, dg_dtheta_b_fit];
                A_dyn_fit = M_inv_fit * dg_dtheta_fit;
                
                A_fit_num = zeros(10, 10);
                B_fit_num = zeros(10, 4);
                
                A_fit_num(1,2) = 1; A_fit_num(3,4) = 1; A_fit_num(5,6) = 1;
                A_fit_num(7,8) = 1; A_fit_num(9,10) = 1;
                
                A_fit_num(2,5) = A_dyn_fit(1,2); A_fit_num(2,7) = A_dyn_fit(1,3); A_fit_num(2,9) = A_dyn_fit(1,4);
                A_fit_num(4,5) = A_dyn_fit(2,2); A_fit_num(4,7) = A_dyn_fit(2,3); A_fit_num(4,9) = A_dyn_fit(2,4);
                A_fit_num(6,5) = A_dyn_fit(3,2); A_fit_num(6,7) = A_dyn_fit(3,3); A_fit_num(6,9) = A_dyn_fit(3,4);
                A_fit_num(8,5) = A_dyn_fit(4,2); A_fit_num(8,7) = A_dyn_fit(4,3); A_fit_num(8,9) = A_dyn_fit(4,4);
                A_fit_num(10,5) = A_dyn_fit(5,2); A_fit_num(10,7) = A_dyn_fit(5,3); A_fit_num(10,9) = A_dyn_fit(5,4);
                
                B_dyn_fit = M_inv_fit * B_eq_fit;
                B_fit_num(2,:) = B_dyn_fit(1,:); B_fit_num(4,:) = B_dyn_fit(2,:);
                B_fit_num(6,:) = B_dyn_fit(3,:); B_fit_num(8,:) = B_dyn_fit(4,:);
                B_fit_num(10,:) = B_dyn_fit(5,:);
                
                %因为matlab是按照列优先展开，这里转置后再展开可以切换为行优先
                K_fit = lqr(A_fit_num, B_fit_num, lqr_Q, lqr_R);
                K_fit_trans = K_fit';
                
                K_sample_2d(idx, 1) = l_l_fit;
                K_sample_2d(idx, 2) = l_r_fit;
                K_sample_2d(idx, 3:42) = K_fit_trans(:)';
                
                % 存储平衡点偏移量
                offset_sample_2d(idx, 1) = l_l_fit;
                offset_sample_2d(idx, 2) = l_r_fit;
                offset_sample_2d(idx, 3) = theta_star_fit(1);  % theta_l*
                offset_sample_2d(idx, 4) = theta_star_fit(2);  % theta_r*
                offset_sample_2d(idx, 5) = theta_star_fit(3);  % theta_b*
            catch ME
                warning('计算失败 l_l=%.2f, l_r=%.2f: %s', l_l_fit, l_r_fit, ME.message);
            end
            
            if mod(idx, 49) == 0
                fprintf('  进度: %d/%d (%.1f秒)\n', idx, sample_size_2d, toc(tic_fit));
            end
        end
    end
    
    fprintf('  ✓ %d 个样本计算完成! 耗时: %.2f秒\n', sample_size_2d, toc(tic_fit));
    
    % 二维多项式拟合 - K矩阵
    fprintf('\n正在进行二维多项式拟合 (K矩阵)...\n');
    
    K_Fit_Coefficients = zeros(40, 6);
    l_l_samples = K_sample_2d(:, 1);
    l_r_samples = K_sample_2d(:, 2);
    
    for n_fit = 1:40
        K_values = K_sample_2d(:, n_fit+2);
        try
            K_Surface_Fit = fit([l_l_samples, l_r_samples], K_values, 'poly22');
            K_Fit_Coefficients(n_fit, :) = coeffvalues(K_Surface_Fit);
        catch
            warning('K拟合失败: 元素 %d', n_fit);
        end
    end
    
    fprintf('  ✓ K矩阵拟合完成\n');
    
    % 二维多项式拟合 - 平衡点偏移量
    fprintf('正在进行二维多项式拟合 (平衡点偏移量)...\n');
    
    Offset_Fit_Coefficients = zeros(3, 6);  % [theta_l*, theta_r*, theta_b*] 各6个系数
    offset_names = {'theta_l_eq', 'theta_r_eq', 'theta_b_eq'};
    
    for n_off = 1:3
        offset_values = offset_sample_2d(:, n_off+2);
        try
            Offset_Surface_Fit = fit([l_l_samples, l_r_samples], offset_values, 'poly22');
            Offset_Fit_Coefficients(n_off, :) = coeffvalues(Offset_Surface_Fit);
        catch
            warning('偏移量拟合失败: %s', offset_names{n_off});
        end
    end
    
    fprintf('  ✓ 偏移量拟合完成\n\n');
    
    % 输出K矩阵拟合系数
    fprintf('// ═══════════════════════════════════════════════════════════════════════\n');
    fprintf('// K矩阵拟合系数 K_Fit_Coefficients[40][6]\n');
    fprintf('// ═══════════════════════════════════════════════════════════════════════\n');
    fprintf('// 拟合多项式: K_ij(l_l, l_r) = p00 + p10*l_l + p01*l_r + p20*l_l^2 + p11*l_l*l_r + p02*l_r^2\n');
    fprintf('// 系数顺序: [p00, p10, p01, p20, p11, p02]\n');
    fprintf('// ═══════════════════════════════════════════════════════════════════════\n\n');
    
    fprintf('float K_Fit_Coefficients[40][6] = {\n');
    for n_fit = 1:40
        row = ceil(n_fit/10) - 1;
        col = mod(n_fit-1, 10);
        fprintf('    {');
        for c = 1:6
            if c < 6
                fprintf('%12.6gf, ', K_Fit_Coefficients(n_fit,c));
            else
                fprintf('%12.6gf', K_Fit_Coefficients(n_fit,c));
            end
        end
        if n_fit < 40
            fprintf('},  // K[%d][%d]\n', row, col);
        else
            fprintf('}   // K[%d][%d]\n', row, col);
        end
    end
    fprintf('};\n\n');
    
    % 输出平衡点偏移量拟合系数
    fprintf('// ═══════════════════════════════════════════════════════════════════════\n');
    fprintf('// 平衡点偏移量拟合系数 Offset_Fit_Coefficients[3][6]\n');
    fprintf('// ═══════════════════════════════════════════════════════════════════════\n');
    fprintf('// 拟合多项式: theta_eq(l_l, l_r) = p00 + p10*l_l + p01*l_r + p20*l_l^2 + p11*l_l*l_r + p02*l_r^2\n');
    fprintf('// 系数顺序: [p00, p10, p01, p20, p11, p02]\n');
    fprintf('// ═══════════════════════════════════════════════════════════════════════\n\n');
    
    fprintf('float Offset_Fit_Coefficients[3][6] = {\n');
    for n_off = 1:3
        fprintf('    {');
        for c = 1:6
            if c < 6
                fprintf('%12.6gf, ', Offset_Fit_Coefficients(n_off,c));
            else
                fprintf('%12.6gf', Offset_Fit_Coefficients(n_off,c));
            end
        end
        if n_off < 3
            fprintf('},  // %s\n', offset_names{n_off});
        else
            fprintf('}   // %s\n', offset_names{n_off});
        end
    end
    fprintf('};\n\n');
    
    % 输出使用说明
    fprintf('// ═══════════════════════════════════════════════════════════════════════\n');
    fprintf('// 使用方法 (C代码示例)\n');
    fprintf('// ═══════════════════════════════════════════════════════════════════════\n');
    fprintf('// // 1. 计算平衡点偏移量\n');
    fprintf('// float theta_l_eq = Offset_Fit_Coefficients[0][0] + Offset_Fit_Coefficients[0][1]*l_l + ...;\n');
    fprintf('// float theta_r_eq = Offset_Fit_Coefficients[1][0] + Offset_Fit_Coefficients[1][1]*l_l + ...;\n');
    fprintf('// float theta_b_eq = Offset_Fit_Coefficients[2][0] + Offset_Fit_Coefficients[2][1]*l_l + ...;\n');
    fprintf('// \n');
    fprintf('// // 2. 计算状态误差 (相对于平衡点)\n');
    fprintf('// float X_err[10] = {X[0], X[1], X[2], X[3], X[4]-theta_l_eq, X[5], X[6]-theta_r_eq, X[7], X[8]-theta_b_eq, X[9]};\n');
    fprintf('// \n');
    fprintf('// // 3. 计算控制输出\n');
    fprintf('// u[i] = -sum(K[i][j] * X_err[j]);\n');
    fprintf('// ═══════════════════════════════════════════════════════════════════════\n\n');
    
    save('lqr_fitting_results.mat', 'K_sample_2d', 'K_Fit_Coefficients', ...
         'offset_sample_2d', 'Offset_Fit_Coefficients', 'Leg_data');
    fprintf('拟合结果已保存到 lqr_fitting_results.mat\n');
    fprintf('  包含: K矩阵拟合系数, 平衡点偏移量拟合系数\n');
end

%% ========================================
%  Part 14: 保存结果
%  ========================================

fprintf('\n========================================\n');
fprintf('Part 14: 保存结果\n');
fprintf('========================================\n\n');

save('lqr_results.mat', 'A_num', 'B_num', 'K', 'lqr_Q', 'lqr_R', 'e', ...
     'theta_l_star', 'theta_r_star', 'theta_b_star', 'M_eq', 'B_eq', 'g_eq');

fprintf('✓ LQR结果已保存到 lqr_results.mat\n');
fprintf('  包含: A_num, B_num, K, 平衡点 theta_*_star\n');

elapsed_time = toc;
fprintf('\n========================================\n');
fprintf('计算完成! 总耗时: %.2f秒\n', elapsed_time);
fprintf('========================================\n');
