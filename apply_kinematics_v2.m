% apply_kinematics_v2.m
% 应用运动学约束 - 将加速度用广义坐标表示
%
% 输入: dynamics_v2.mat (5个化简后的动力学方程，对应推导 §3)
% 输出: dynamics_new_coords.mat (用广义坐标表示的方程)
%
% 对照推导文档:
%   §4.1 运动学约束
%   §4.2 代入运动学约束后的加速度表达式
%   §5.0 广义坐标定义

clear; clc;
fprintf('========================================\n');
fprintf('应用运动学约束 (对照推导 §4)\n');
fprintf('========================================\n\n');

%% 加载化简后的方程

fprintf('Step 1: 加载化简后的方程...\n');
load('dynamics_v2.mat');
fprintf('  已加载 dynamics_v2.mat (推导 §3 的5个方程)\n\n');

%% 定义符号变量

fprintf('Step 2: 定义符号变量...\n\n');

% ========== 物理参数 ==========
syms g real                                    % 重力加速度
syms R R_w real                                % 轮半径, 半轮距
syms m_b m_l m_r m_wl m_wr real               % 质量
syms I_b I_l I_r I_wl I_wr I_yaw real         % 转动惯量
syms l_b l_l l_r l_l_d l_r_d real             % 长度参数
syms theta_l0 theta_r0 theta_b0 real          % 初始角偏置

% ========== 控制输入 (推导 §1.4) ==========
syms T_wr_to_r T_wl_to_l T_r_to_b T_l_to_b real

% ========== 广义坐标及其导数 (推导 §5.0) ==========
% q = [X_b^h, phi, theta_l, theta_r, theta_b]
syms X_b_h dX_b_h ddX_b_h real     % 机体水平位置
syms phi dphi ddphi real           % yaw角
syms theta_l dtheta_l ddtheta_l real
syms theta_r dtheta_r ddtheta_r real
syms theta_b dtheta_b ddtheta_b real

% ========== 轮角速度 ==========
syms theta_wl dtheta_wl ddtheta_wl real
syms theta_wr dtheta_wr ddtheta_wr real

% ========== 加速度变量 (原方程中出现的) ==========
syms a_b_h a_b_v real              % 机体加速度
syms a_l_h a_l_v a_r_h a_r_v real  % 腿转轴加速度
syms a_wl_h a_wl_v a_wr_h a_wr_v real  % 轮加速度

%% ========================================
%  §4.1 运动学约束
%  ========================================

fprintf('========================================\n');
fprintf('§4.1 运动学约束\n');
fprintf('========================================\n\n');

% ----- 纯滚动约束 -----
fprintf('【纯滚动约束】\n');
fprintf('  a_wr^h = R * ddtheta_wr\n');
fprintf('  a_wl^h = R * ddtheta_wl\n\n');

a_wr_h_expr = R * ddtheta_wr;
a_wl_h_expr = R * ddtheta_wl;

% ----- 轮不离地约束 -----
fprintf('【轮不离地约束】\n');
fprintf('  a_wr^v = 0\n');
fprintf('  a_wl^v = 0\n\n');

a_wr_v_expr = sym(0);
a_wl_v_expr = sym(0);

% ----- 腿转轴加速度 (由轮得到) -----
fprintf('【腿转轴加速度】\n');
fprintf('  a_r^h = a_wr^h + l_r*cos(theta_r)*ddtheta_r - l_r*sin(theta_r)*dtheta_r^2\n');
fprintf('  a_r^v = a_wr^v - l_r*sin(theta_r)*ddtheta_r - l_r*cos(theta_r)*dtheta_r^2\n');
fprintf('  a_l^h = a_wl^h + l_l*cos(theta_l)*ddtheta_l - l_l*sin(theta_l)*dtheta_l^2\n');
fprintf('  a_l^v = a_wl^v - l_l*sin(theta_l)*ddtheta_l - l_l*cos(theta_l)*dtheta_l^2\n\n');

a_r_h_expr = a_wr_h_expr + l_r*cos(theta_r)*ddtheta_r - l_r*sin(theta_r)*dtheta_r^2;
a_r_v_expr = a_wr_v_expr - l_r*sin(theta_r)*ddtheta_r - l_r*cos(theta_r)*dtheta_r^2;
a_l_h_expr = a_wl_h_expr + l_l*cos(theta_l)*ddtheta_l - l_l*sin(theta_l)*dtheta_l^2;
a_l_v_expr = a_wl_v_expr - l_l*sin(theta_l)*ddtheta_l - l_l*cos(theta_l)*dtheta_l^2;

% ----- 机体加速度 (左右腿转轴平均) -----
fprintf('【机体加速度】\n');
fprintf('  a_b^h = (a_r^h + a_l^h) / 2\n');
fprintf('  a_b^v = (a_r^v + a_l^v) / 2\n\n');

a_b_h_expr = (a_r_h_expr + a_l_h_expr) / 2;
a_b_v_expr = (a_r_v_expr + a_l_v_expr) / 2;

% ----- Yaw角加速度 -----
fprintf('【Yaw角加速度】\n');
fprintf('  ddphi = (a_wr^h - a_wl^h) / (2*R_w) = R*(ddtheta_wr - ddtheta_wl) / (2*R_w)\n\n');

ddphi_expr = (a_wr_h_expr - a_wl_h_expr) / (2*R_w);

%% ========================================
%  §4.2 代入运动学约束后的加速度表达式
%  ========================================

fprintf('========================================\n');
fprintf('§4.2 展开后的加速度表达式\n');
fprintf('========================================\n\n');

% 展开 a_b_h 表达式
a_b_h_expanded = simplify(a_b_h_expr);
fprintf('a_b^h (展开) = '); disp(a_b_h_expanded);

% 展开 a_b_v 表达式
a_b_v_expanded = simplify(a_b_v_expr);
fprintf('a_b^v (展开) = '); disp(a_b_v_expanded);

%% ========================================
%  §5.0 轮角加速度用广义坐标表示
%  ========================================

fprintf('========================================\n');
fprintf('§5.0 轮角加速度用广义坐标表示\n');
fprintf('========================================\n\n');

% 从 a_b_h = ddX_b_h 和 ddphi 的关系反解 ddtheta_wr, ddtheta_wl
%
% 由 §4.2:
%   ddX_b_h = R*(ddtheta_wr + ddtheta_wl)/2 + (腿角加速度项) - (腿角速度平方项)
%   ddphi = R*(ddtheta_wr - ddtheta_wl)/(2*R_w)
%
% 令:
%   leg_term = (l_r*cos(theta_r)*ddtheta_r + l_l*cos(theta_l)*ddtheta_l)/2
%            - (l_r*sin(theta_r)*dtheta_r^2 + l_l*sin(theta_l)*dtheta_l^2)/2
%
% 则:
%   ddX_b_h = R*(ddtheta_wr + ddtheta_wl)/2 + leg_term
%   ddphi = R*(ddtheta_wr - ddtheta_wl)/(2*R_w)
%
% 解得:
%   ddtheta_wr = (ddX_b_h + R_w*ddphi)/R - leg_term/R
%   ddtheta_wl = (ddX_b_h - R_w*ddphi)/R - leg_term/R

leg_accel_term = (l_r*cos(theta_r)*ddtheta_r + l_l*cos(theta_l)*ddtheta_l)/2;
leg_vel_sq_term = (l_r*sin(theta_r)*dtheta_r^2 + l_l*sin(theta_l)*dtheta_l^2)/2;
leg_term = leg_accel_term - leg_vel_sq_term;

ddtheta_wr_expr = (ddX_b_h + R_w*ddphi)/R - leg_term/R;
ddtheta_wl_expr = (ddX_b_h - R_w*ddphi)/R - leg_term/R;

fprintf('ddtheta_wr = (ddX_b_h + R_w*ddphi)/R - leg_term/R\n');
fprintf('ddtheta_wl = (ddX_b_h - R_w*ddphi)/R - leg_term/R\n');
fprintf('其中 leg_term = (l_r*cos(theta_r)*ddtheta_r + l_l*cos(theta_l)*ddtheta_l)/2\n');
fprintf('             - (l_r*sin(theta_r)*dtheta_r^2 + l_l*sin(theta_l)*dtheta_l^2)/2\n\n');

%% ========================================
%  Step 3: 构建代换规则
%  ========================================

fprintf('========================================\n');
fprintf('Step 3: 构建代换规则\n');
fprintf('========================================\n\n');

% 将原方程中的加速度变量替换为用广义坐标表示的表达式
% 注意: a_b_h 不直接替换为 ddX_b_h，而是保持为展开形式以便后续提取系数

kinematics_subs = {
    % 轮加速度 (纯滚动)
    a_wr_h, R*ddtheta_wr;
    a_wl_h, R*ddtheta_wl;
    a_wr_v, sym(0);
    a_wl_v, sym(0);
    
    % 腿转轴加速度
    a_r_h, R*ddtheta_wr + l_r*cos(theta_r)*ddtheta_r - l_r*sin(theta_r)*dtheta_r^2;
    a_r_v, -l_r*sin(theta_r)*ddtheta_r - l_r*cos(theta_r)*dtheta_r^2;
    a_l_h, R*ddtheta_wl + l_l*cos(theta_l)*ddtheta_l - l_l*sin(theta_l)*dtheta_l^2;
    a_l_v, -l_l*sin(theta_l)*ddtheta_l - l_l*cos(theta_l)*dtheta_l^2;
    
    % 机体加速度 (左右腿转轴平均)
    a_b_h, (R*(ddtheta_wr + ddtheta_wl)/2 ...
          + (l_r*cos(theta_r)*ddtheta_r + l_l*cos(theta_l)*ddtheta_l)/2 ...
          - (l_r*sin(theta_r)*dtheta_r^2 + l_l*sin(theta_l)*dtheta_l^2)/2);
    a_b_v, (-(l_r*sin(theta_r)*ddtheta_r + l_l*sin(theta_l)*ddtheta_l)/2 ...
          - (l_r*cos(theta_r)*dtheta_r^2 + l_l*cos(theta_l)*dtheta_l^2)/2);
};

fprintf('代换规则:\n');
for i = 1:size(kinematics_subs, 1)
    fprintf('  %s -> (表达式)\n', char(kinematics_subs{i,1}));
end
fprintf('\n');

%% ========================================
%  Step 4: 对5个方程进行代换
%  ========================================

fprintf('========================================\n');
fprintf('Step 4: 代换动力学方程 (推导 §3)\n');
fprintf('========================================\n\n');

% eq1: §3.1 整体水平动量方程
% eq2: §3.2 机体转动方程
% eq3: §3.3 右腿转动方程
% eq4: §3.4 左腿转动方程
% eq5: §3.5 Yaw转动方程

eq1_new = subs(eq1, kinematics_subs(:,1), kinematics_subs(:,2));
eq2_new = subs(eq2, kinematics_subs(:,1), kinematics_subs(:,2));
eq3_new = subs(eq3, kinematics_subs(:,1), kinematics_subs(:,2));
eq4_new = subs(eq4, kinematics_subs(:,1), kinematics_subs(:,2));
eq5_new = subs(eq5, kinematics_subs(:,1), kinematics_subs(:,2));

% 化简
fprintf('正在化简方程...\n');
eq1_new = simplify(eq1_new);
eq2_new = simplify(eq2_new);
eq3_new = simplify(eq3_new);
eq4_new = simplify(eq4_new);
eq5_new = simplify(eq5_new);

fprintf('\n代换后的方程:\n\n');

fprintf('方程1 (§3.1 整体水平动量):\n');
disp(eq1_new);

fprintf('方程2 (§3.2 机体转动):\n');
disp(eq2_new);

fprintf('方程3 (§3.3 右腿转动):\n');
disp(eq3_new);

fprintf('方程4 (§3.4 左腿转动):\n');
disp(eq4_new);

fprintf('方程5 (§3.5 Yaw转动):\n');
disp(eq5_new);

%% ========================================
%  验证运动学约束的正确性
%  ========================================

fprintf('========================================\n');
fprintf('验证运动学约束\n');
fprintf('========================================\n\n');

% 验证 a_b_h 展开式与 §4.2 一致
a_b_h_ref = R*(ddtheta_wr + ddtheta_wl)/2 ...
          + (l_r*cos(theta_r)*ddtheta_r + l_l*cos(theta_l)*ddtheta_l)/2 ...
          - (l_r*sin(theta_r)*dtheta_r^2 + l_l*sin(theta_l)*dtheta_l^2)/2;

diff_h = simplify(a_b_h_expr - a_b_h_ref);
if isequal(diff_h, sym(0))
    fprintf('✓ a_b^h 展开式与 §4.2 一致\n');
else
    fprintf('✗ a_b^h 差异: '); disp(diff_h);
end

% 验证 a_b_v 展开式
a_b_v_ref = -(l_r*sin(theta_r)*ddtheta_r + l_l*sin(theta_l)*ddtheta_l)/2 ...
          - (l_r*cos(theta_r)*dtheta_r^2 + l_l*cos(theta_l)*dtheta_l^2)/2;

diff_v = simplify(a_b_v_expr - a_b_v_ref);
if isequal(diff_v, sym(0))
    fprintf('✓ a_b^v 展开式与 §4.2 一致\n');
else
    fprintf('✗ a_b^v 差异: '); disp(diff_v);
end

% 验证 ddphi 表达式
ddphi_ref = R*(ddtheta_wr - ddtheta_wl)/(2*R_w);
diff_phi = simplify(ddphi_expr - ddphi_ref);
if isequal(diff_phi, sym(0))
    fprintf('✓ ddphi 表达式与 §4.1 一致\n');
else
    fprintf('✗ ddphi 差异: '); disp(diff_phi);
end

fprintf('\n');

%% ========================================
%  汇总
%  ========================================

fprintf('========================================\n');
fprintf('运动学约束应用结果汇总\n');
fprintf('========================================\n');
fprintf('✓ §4.1 纯滚动约束: 已应用\n');
fprintf('✓ §4.1 轮不离地约束: 已应用\n');
fprintf('✓ §4.1 腿转轴加速度: 已应用\n');
fprintf('✓ §4.1 机体加速度: 已应用\n');
fprintf('✓ §4.1 Yaw角加速度: 已应用\n');
fprintf('✓ §5.0 轮角加速度表达式: 已建立\n');
fprintf('----------------------------------------\n');
fprintf('5个动力学方程已转换为广义坐标表示\n');
fprintf('========================================\n\n');

%% ========================================
%  保存结果
%  ========================================

fprintf('Step 5: 保存结果...\n');

save('dynamics_new_coords.mat', ...
     'eq1_new', 'eq2_new', 'eq3_new', 'eq4_new', 'eq5_new', ...
     'kinematics_subs', ...
     'ddtheta_wr_expr', 'ddtheta_wl_expr', ...
     'leg_term');

fprintf('  结果已保存到 dynamics_new_coords.mat\n');
fprintf('\n========================================\n');
fprintf('运动学约束应用完成!\n');
fprintf('========================================\n');
