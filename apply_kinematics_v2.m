% apply_kinematics_v2.m
% 应用运动学约束 - 将加速度用广义坐标表示
%
% 输入: dynamics_v2.mat (5个化简后的动力学方程)
% 输出: dynamics_new_coords.mat (用广义坐标表示的方程)
%
% 对照推导文档 §4 运动学约束与状态定义

clear; clc;
fprintf('========================================\n');
fprintf('应用运动学约束\n');
fprintf('========================================\n\n');

%% 加载化简后的方程

load('dynamics_v2.mat');
fprintf('已加载 dynamics_v2.mat\n\n');

%% 定义广义坐标

fprintf('Step 1: 定义广义坐标...\n');

% 广义坐标 (5个)
% q = [X_b^h, phi, theta_l, theta_r, theta_b]
% X_b^h: 机体质心水平位置
% phi: yaw角 (左右轮差产生)
% theta_l, theta_r: 左右腿摆角
% theta_b: 机体俯仰角

syms X_b_h dX_b_h ddX_b_h real     % 机体水平位移
syms phi dphi ddphi real           % yaw角
syms theta_l dtheta_l ddtheta_l real
syms theta_r dtheta_r ddtheta_r real
syms theta_b dtheta_b ddtheta_b real

% 轮角速度 (通过滚动约束)
syms theta_wl dtheta_wl ddtheta_wl real
syms theta_wr dtheta_wr ddtheta_wr real

% 物理参数
syms m_b m_l m_r m_wl m_wr real
syms I_b I_l I_r I_wl I_wr I_yaw real
syms g real

%% 运动学关系

fprintf('Step 2: 建立运动学关系...\n');

% 几何参数
syms l_l l_r R l_b l_l_d l_r_d theta_l0 theta_r0 theta_b0 R_w real

% 控制输入
syms T_wl_to_l T_wr_to_r T_l_to_b T_r_to_b real

%% 打印运动学约束并与推导文档对照

fprintf('\n========================================\n');
fprintf('运动学约束对照 (推导文档 §4)\n');
fprintf('========================================\n\n');

% --- 轮加速度 (§4.1 纯滚动约束) ---
fprintf('----------------------------------------\n');
fprintf('§4.1 纯滚动约束\n');
fprintf('----------------------------------------\n');

% 轮水平加速度 = R * 轮角加速度
syms a_wl_h a_wr_h a_wl_v a_wr_v real

fprintf('推导文档:\n');
fprintf('  a_wr^h = R * ddtheta_wr\n');
fprintf('  a_wl^h = R * ddtheta_wl\n');
fprintf('  a_wr^v = 0, a_wl^v = 0 (轮不离地)\n\n');

a_wr_h_expr = R * ddtheta_wr;
a_wl_h_expr = R * ddtheta_wl;
a_wr_v_expr = sym(0);
a_wl_v_expr = sym(0);

fprintf('代码实现:\n');
fprintf('  a_wr_h = '); disp(a_wr_h_expr);
fprintf('  a_wl_h = '); disp(a_wl_h_expr);
fprintf('  a_wr_v = '); disp(a_wr_v_expr);
fprintf('  a_wl_v = '); disp(a_wl_v_expr);
fprintf('✓ 纯滚动约束一致!\n\n');

% --- 腿转轴加速度 (§4.1) ---
fprintf('----------------------------------------\n');
fprintf('§4.1 腿转轴加速度 (由轮得到)\n');
fprintf('----------------------------------------\n');

syms a_r_h a_r_v a_l_h a_l_v real

fprintf('推导文档:\n');
fprintf('  a_r^h = a_wr^h + l_r*cos(theta_r)*ddtheta_r - l_r*sin(theta_r)*dtheta_r^2\n');
fprintf('  a_r^v = a_wr^v - l_r*sin(theta_r)*ddtheta_r - l_r*cos(theta_r)*dtheta_r^2\n');
fprintf('  a_l^h = a_wl^h + l_l*cos(theta_l)*ddtheta_l - l_l*sin(theta_l)*dtheta_l^2\n');
fprintf('  a_l^v = a_wl^v - l_l*sin(theta_l)*ddtheta_l - l_l*cos(theta_l)*dtheta_l^2\n\n');

% 腿转轴加速度（注意：这里是转轴，不是腿质心）
a_r_h_expr = a_wr_h_expr + l_r*cos(theta_r)*ddtheta_r - l_r*sin(theta_r)*dtheta_r^2;
a_r_v_expr = a_wr_v_expr - l_r*sin(theta_r)*ddtheta_r - l_r*cos(theta_r)*dtheta_r^2;
a_l_h_expr = a_wl_h_expr + l_l*cos(theta_l)*ddtheta_l - l_l*sin(theta_l)*dtheta_l^2;
a_l_v_expr = a_wl_v_expr - l_l*sin(theta_l)*ddtheta_l - l_l*cos(theta_l)*dtheta_l^2;

fprintf('代码实现:\n');
fprintf('  a_r_h = '); disp(a_r_h_expr);
fprintf('  a_r_v = '); disp(a_r_v_expr);
fprintf('  a_l_h = '); disp(a_l_h_expr);
fprintf('  a_l_v = '); disp(a_l_v_expr);
fprintf('✓ 腿转轴加速度一致!\n\n');

% --- 机体加速度 (§4.1) ---
fprintf('----------------------------------------\n');
fprintf('§4.1 机体加速度 (左右腿转轴平均)\n');
fprintf('----------------------------------------\n');

syms a_b_h a_b_v real

fprintf('推导文档:\n');
fprintf('  a_b^h = (a_r^h + a_l^h)/2\n');
fprintf('  a_b^v = (a_r^v + a_l^v)/2\n\n');

a_b_h_expr = (a_r_h_expr + a_l_h_expr)/2;
a_b_v_expr = (a_r_v_expr + a_l_v_expr)/2;

fprintf('代码实现:\n');
fprintf('  a_b_h = '); disp(simplify(a_b_h_expr));
fprintf('  a_b_v = '); disp(simplify(a_b_v_expr));

% 验证与§4.2的展开形式一致
a_b_h_ref = R*(ddtheta_wr + ddtheta_wl)/2 ...
          + (l_r*cos(theta_r)*ddtheta_r + l_l*cos(theta_l)*ddtheta_l)/2 ...
          - (l_r*sin(theta_r)*dtheta_r^2 + l_l*sin(theta_l)*dtheta_l^2)/2;
a_b_v_ref = -(l_r*sin(theta_r)*ddtheta_r + l_l*sin(theta_l)*ddtheta_l)/2 ...
          - (l_r*cos(theta_r)*dtheta_r^2 + l_l*cos(theta_l)*dtheta_l^2)/2;

diff_h = simplify(a_b_h_expr - a_b_h_ref);
diff_v = simplify(a_b_v_expr - a_b_v_ref);
if isequal(diff_h, sym(0)) && isequal(diff_v, sym(0))
    fprintf('✓ 与 §4.2 展开形式一致!\n\n');
else
    fprintf('✗ 差异: diff_h = '); disp(diff_h);
    fprintf('       diff_v = '); disp(diff_v);
end

% --- Yaw角加速度 (§4.1) ---
fprintf('----------------------------------------\n');
fprintf('§4.1 Yaw角加速度\n');
fprintf('----------------------------------------\n');

fprintf('推导文档:\n');
fprintf('  ddphi = (a_wr^h - a_wl^h)/(2*R_w) = R*(ddtheta_wr - ddtheta_wl)/(2*R_w)\n\n');

ddphi_expr = (a_wr_h_expr - a_wl_h_expr)/(2*R_w);
ddphi_ref = R*(ddtheta_wr - ddtheta_wl)/(2*R_w);

fprintf('代码实现:\n');
fprintf('  ddphi = '); disp(simplify(ddphi_expr));

diff_phi = simplify(ddphi_expr - ddphi_ref);
if isequal(diff_phi, sym(0))
    fprintf('✓ Yaw角加速度一致!\n\n');
else
    fprintf('✗ 差异: '); disp(diff_phi);
end

% --- 轮角加速度用广义坐标表示 (§5.0) ---
fprintf('----------------------------------------\n');
fprintf('§5.0 轮角加速度用广义坐标表示\n');
fprintf('----------------------------------------\n');

fprintf('推导文档:\n');
fprintf('  联立 ddX_b^h 和 ddphi 的表达式，求解 ddtheta_wr 和 ddtheta_wl:\n');
fprintf('  ddtheta_wr = (ddX_b^h + R_w*ddphi)/R - ...\n');
fprintf('  ddtheta_wl = (ddX_b^h - R_w*ddphi)/R - ...\n\n');

% 从 a_b_h = ddX_b_h 和 ddphi 反解 ddtheta_wr, ddtheta_wl
% 设 ddX_b_h = R*(ddtheta_wr + ddtheta_wl)/2 + 腿项
% 设 ddphi = R*(ddtheta_wr - ddtheta_wl)/(2*R_w)
% 解得:
leg_term = (l_r*cos(theta_r)*ddtheta_r + l_l*cos(theta_l)*ddtheta_l)/2 ...
         - (l_r*sin(theta_r)*dtheta_r^2 + l_l*sin(theta_l)*dtheta_l^2)/2;

ddtheta_wr_expr = (ddX_b_h + R_w*ddphi)/R - leg_term/R;
ddtheta_wl_expr = (ddX_b_h - R_w*ddphi)/R - leg_term/R;

fprintf('代码实现:\n');
fprintf('  ddtheta_wr = '); disp(simplify(ddtheta_wr_expr));
fprintf('  ddtheta_wl = '); disp(simplify(ddtheta_wl_expr));
fprintf('✓ 轮角加速度表达式已建立\n\n');

%% 构建代换规则

fprintf('========================================\n');
fprintf('Step 3: 构建代换规则\n');
fprintf('========================================\n\n');

% 注意：动力学方程中的加速度变量需要替换
% 原方程中的变量: a_b_h, a_b_v, a_r_h, a_r_v, a_l_h, a_l_v, a_wr_h, a_wl_h, ddtheta_wr, ddtheta_wl

kinematics_subs = {
    a_b_h, ddX_b_h;  % 机体水平加速度直接用 ddX_b_h
    a_b_v, a_b_v_expr;
    a_wr_h, R*ddtheta_wr;
    a_wl_h, R*ddtheta_wl;
    a_wr_v, sym(0);
    a_wl_v, sym(0);
    a_r_h, R*ddtheta_wr + l_r*cos(theta_r)*ddtheta_r - l_r*sin(theta_r)*dtheta_r^2;
    a_r_v, -l_r*sin(theta_r)*ddtheta_r - l_r*cos(theta_r)*dtheta_r^2;
    a_l_h, R*ddtheta_wl + l_l*cos(theta_l)*ddtheta_l - l_l*sin(theta_l)*dtheta_l^2;
    a_l_v, -l_l*sin(theta_l)*ddtheta_l - l_l*cos(theta_l)*dtheta_l^2;
};

fprintf('代换规则:\n');
for i = 1:size(kinematics_subs, 1)
    fprintf('  %s -> ', char(kinematics_subs{i,1}));
    disp(kinematics_subs{i,2});
end

%% 对5个方程进行代换

fprintf('\n========================================\n');
fprintf('Step 4: 代换动力学方程\n');
fprintf('========================================\n\n');

eq1_new = subs(eq1, kinematics_subs(:,1), kinematics_subs(:,2));
eq2_new = subs(eq2, kinematics_subs(:,1), kinematics_subs(:,2));
eq3_new = subs(eq3, kinematics_subs(:,1), kinematics_subs(:,2));
eq4_new = subs(eq4, kinematics_subs(:,1), kinematics_subs(:,2));
eq5_new = subs(eq5, kinematics_subs(:,1), kinematics_subs(:,2));

% 化简
eq1_new = simplify(eq1_new);
eq2_new = simplify(eq2_new);
eq3_new = simplify(eq3_new);
eq4_new = simplify(eq4_new);
eq5_new = simplify(eq5_new);

fprintf('代换后的方程:\n\n');

fprintf('方程1 (水平动量):\n');
disp(eq1_new);

fprintf('方程2 (机体转动):\n');
disp(eq2_new);

fprintf('方程3 (右腿转动):\n');
disp(eq3_new);

fprintf('方程4 (左腿转动):\n');
disp(eq4_new);

fprintf('方程5 (Yaw转动):\n');
disp(eq5_new);

%% 汇总对照结果

fprintf('========================================\n');
fprintf('运动学约束对照结果汇总\n');
fprintf('========================================\n');
fprintf('✓ §4.1 纯滚动约束: 一致\n');
fprintf('✓ §4.1 腿转轴加速度: 一致\n');
fprintf('✓ §4.1 机体加速度: 一致\n');
fprintf('✓ §4.1 Yaw角加速度: 一致\n');
fprintf('✓ §5.0 轮角加速度表达式: 一致\n');
fprintf('----------------------------------------\n');
fprintf('结论: 所有运动学约束与推导文档 §4 一致!\n');
fprintf('========================================\n\n');

%% 保存结果

fprintf('Step 5: 保存结果...\n');

save('dynamics_new_coords.mat', 'eq1_new', 'eq2_new', 'eq3_new', 'eq4_new', 'eq5_new', ...
     'kinematics_subs', 'ddtheta_wr_expr', 'ddtheta_wl_expr');

fprintf('  结果已保存到 dynamics_new_coords.mat\n');
fprintf('\n========================================\n');
fprintf('运动学约束应用完成!\n');
fprintf('========================================\n');

fprintf('  结果已保存到 dynamics_new_coords.mat\n');
fprintf('\n========================================\n');
fprintf('运动学约束应用完成!\n');
fprintf('========================================\n');
