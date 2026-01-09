% compute_lqr.m
% 基于新推导的动力学模型计算LQR控制器增益矩阵
% 
% 状态向量 (10维):
%   x = [X_b^h; dX_b^h; phi; dphi; theta_l; dtheta_l; theta_r; dtheta_r; theta_b; dtheta_b]
%   其中:
%     X_b^h    : 机体水平位置 (m)
%     dX_b^h   : 机体水平速度 (m/s)
%     phi      : 偏航角 (rad)
%     dphi     : 偏航角速度 (rad/s)
%     theta_l  : 左腿与Z轴负方向夹角 (rad)
%     dtheta_l : 左腿角速度 (rad/s)
%     theta_r  : 右腿与Z轴负方向夹角 (rad)
%     dtheta_r : 右腿角速度 (rad/s)
%     theta_b  : 机体俯仰角 (rad)
%     dtheta_b : 机体俯仰角速度 (rad/s)
%
% 控制向量 (4维):
%   u = [T_wr; T_wl; T_r; T_l]
%   其中:
%     T_wr : 右轮电机力矩 (Nm)
%     T_wl : 左轮电机力矩 (Nm)
%     T_r  : 右髋关节电机力矩 (Nm)
%     T_l  : 左髋关节电机力矩 (Nm)
%
% 作者: 基于2026公式推导
% 日期: 2026/01/09

clear all; clc;
tic

%% ======================== Step 0: 定义符号变量 ========================

fprintf('========================================\n');
fprintf('轮腿机器人LQR控制器计算 (新动力学模型)\n');
fprintf('========================================\n\n');

% 物理常数
syms g real

% 机体参数
syms m_b I_b l_b theta_b0 real

% 腿部参数 (左右)
syms m_l m_r I_l I_r l_l l_r l_l_d l_r_d theta_l0 theta_r0 real

% 轮子参数 (左右)
syms m_wl m_wr I_wl I_wr real

% 几何参数
syms R R_w I_yaw real

%% ======================== Step 1: 构建线性化A、B矩阵 ========================

fprintf('Step 1: 构建线性化状态空间矩阵...\n');

% 定义公共分母 D (对称假设下)
% D = 2*I_leg*I_w + I_leg*R^2*m_b + 2*I_leg*R^2*m_leg + 2*I_leg*R^2*m_w 
%   + I_w*l_leg^2*m_b + 2*I_w*l_leg^2*m_leg + R^2*l_leg^2*m_b*m_w + 2*R^2*l_leg^2*m_leg*m_w

% 初始化 A 和 B 矩阵
A = sym(zeros(10, 10));
B = sym(zeros(10, 4));

% A矩阵: 位置-速度关系
A(1,2) = 1;  % dX/dt = V
A(3,4) = 1;  % dphi/dt = omega_phi
A(5,6) = 1;  % dtheta_l/dt = omega_l
A(7,8) = 1;  % dtheta_r/dt = omega_r
A(9,10) = 1; % dtheta_b/dt = omega_b

% 定义D (公共分母) - 使用通用左右不对称形式
D_l = 2*I_l*I_wl + I_l*R^2*m_b + 2*I_l*R^2*m_l + 2*I_l*R^2*m_wl ...
    + I_wl*l_l^2*m_b + 2*I_wl*l_l^2*m_l + R^2*l_l^2*m_b*m_wl + 2*R^2*l_l^2*m_l*m_wl;

D_r = 2*I_r*I_wr + I_r*R^2*m_b + 2*I_r*R^2*m_r + 2*I_r*R^2*m_wr ...
    + I_wr*l_r^2*m_b + 2*I_wr*l_r^2*m_r + R^2*l_r^2*m_b*m_wr + 2*R^2*l_r^2*m_r*m_wr;

% 为简化，使用对称假设下的公共D
D = 2*I_l*I_wl + I_l*R^2*m_b + 2*I_l*R^2*m_l + 2*I_l*R^2*m_wl ...
    + I_wl*l_l^2*m_b + 2*I_wl*l_l^2*m_l + R^2*l_l^2*m_b*m_wl + 2*R^2*l_l^2*m_l*m_wl;

D_yaw = I_yaw*R^2 + 2*I_wl*R_w^2;

% A矩阵: 加速度对角度的依赖 (在平衡点theta=0处线性化)
% 重力不稳定项
gravity_term_l = l_l*m_b + 2*l_l*m_l - 2*l_l_d*m_l*cos(theta_l0);
gravity_term_r = l_r*m_b + 2*l_r*m_r - 2*l_r_d*m_r*cos(theta_r0);

% A(2,5), A(2,7): 水平加速度对腿角度的耦合
A(2,5) = -g*l_l*(I_wl + R^2*m_wl)*gravity_term_l / (2*D);
A(2,7) = -g*l_r*(I_wr + R^2*m_wr)*gravity_term_r / (2*D);

% A(6,5): 左腿角加速度对自身角度 (不稳定项)
numerator_65 = 4*I_l*I_wl + 2*I_l*R^2*m_b + 4*I_l*R^2*m_l + 4*I_l*R^2*m_wl ...
             + I_wl*l_l^2*m_b + 2*I_wl*l_l^2*m_l + R^2*l_l^2*m_b*m_wl + 2*R^2*l_l^2*m_l*m_wl;
A(6,5) = -g*gravity_term_l*numerator_65 / (4*I_l*D);

% A(6,7): 左腿角加速度对右腿角度的耦合
A(6,7) = g*l_l*l_r*(I_wl + R^2*m_wl)*(m_b + 2*m_r)*gravity_term_r / (4*I_l*D);

% A(8,5): 右腿角加速度对左腿角度的耦合
A(8,5) = g*l_l*l_r*(I_wr + R^2*m_wr)*(m_b + 2*m_l)*gravity_term_l / (4*I_r*D);

% A(8,7): 右腿角加速度对自身角度 (不稳定项)
numerator_87 = 4*I_r*I_wr + 2*I_r*R^2*m_b + 4*I_r*R^2*m_r + 4*I_r*R^2*m_wr ...
             + I_wr*l_r^2*m_b + 2*I_wr*l_r^2*m_r + R^2*l_r^2*m_b*m_wr + 2*R^2*l_r^2*m_r*m_wr;
A(8,7) = -g*gravity_term_r*numerator_87 / (4*I_r*D);

% A(10,9): 机体俯仰角加速度对自身角度
A(10,9) = g*l_b*m_b*cos(theta_b0) / I_b;

% B矩阵: 控制输入 u = [T_wr; T_wl; T_r; T_l]

% B(2,:): 水平加速度对力矩的响应
B(2,1) = (I_wl*l_l - I_l*R + R^2*l_l*m_wl) / D;  % T_wr对水平加速度
B(2,2) = (I_wr*l_r - I_r*R + R^2*l_r*m_wr) / D;  % T_wl对水平加速度
B(2,3) = -l_r*(I_wr + R^2*m_wr) / D;              % T_r对水平加速度
B(2,4) = -l_l*(I_wl + R^2*m_wl) / D;              % T_l对水平加速度

% B(4,:): yaw角加速度对力矩的响应
B(4,1) = -R*R_w / D_yaw;   % T_wr对yaw
B(4,2) = R*R_w / D_yaw;    % T_wl对yaw (反向)

% B(6,:): 左腿角加速度对力矩的响应
B(6,3) = l_l*l_r*(I_wl + R^2*m_wl)*(m_b + 2*m_r) / (2*I_l*D);  % T_r
B(6,4) = -(numerator_65) / (2*I_l*D);  % T_l (直接控制)

% B(8,:): 右腿角加速度对力矩的响应  
B(8,3) = -(numerator_87) / (2*I_r*D);  % T_r (直接控制)
B(8,4) = l_l*l_r*(I_wr + R^2*m_wr)*(m_b + 2*m_l) / (2*I_r*D);  % T_l

% B(10,:): 机体俯仰角加速度对力矩的响应
B(10,3) = 1/I_b;  % T_r
B(10,4) = 1/I_b;  % T_l

% 轮电机对腿角度的影响 (较复杂的耦合项，这里简化处理)
% 完整公式见linearize_system_v2.m输出
Gamma = 2*I_l*I_yaw*R - 2*I_l*R*R_w^2*m_b - 4*I_l*R*R_w^2*m_l - 4*I_l*R*R_w^2*m_wl ...
      + I_yaw*R*l_l^2*m_b + I_yaw*R^2*l_l*m_b + 2*I_yaw*R*l_l^2*m_l + 2*I_yaw*R^2*l_l*m_l ...
      + 2*I_wl*R_w^2*l_l*m_b + 4*I_wl*R_w^2*l_l*m_l - 2*R*R_w^2*l_l^2*m_b*m_wl - 4*R*R_w^2*l_l^2*m_l*m_wl;

B(6,1) = -l_l*(I_wl + R^2*m_wl)*Gamma / (2*I_l*D_yaw*D);

% 对称处理B(8,2)
Gamma_r = 2*I_r*I_yaw*R - 2*I_r*R*R_w^2*m_b - 4*I_r*R*R_w^2*m_r - 4*I_r*R*R_w^2*m_wr ...
        + I_yaw*R*l_r^2*m_b + I_yaw*R^2*l_r*m_b + 2*I_yaw*R*l_r^2*m_r + 2*I_yaw*R^2*l_r*m_r ...
        + 2*I_wr*R_w^2*l_r*m_b + 4*I_wr*R_w^2*l_r*m_r - 2*R*R_w^2*l_r^2*m_b*m_wr - 4*R*R_w^2*l_r^2*m_r*m_wr;

B(8,2) = -l_r*(I_wr + R^2*m_wr)*Gamma_r / (2*I_r*D_yaw*D);

fprintf('  A矩阵和B矩阵符号表达式构建完成\n\n');

%% ======================== Step 2: 输入物理参数 ========================

fprintf('Step 2: 输入机器人物理参数...\n');

% 重力加速度
g_val = 9.81;  % m/s^2

% ==================== 机器人机体参数 ====================
R_val = 0.055;              % 轮子半径 (m)
R_w_val = (0.455+0.445)/4;  % 轮距/2 (m)

m_b_val = 15.564;           % 机体质量 (kg)
I_b_val = 0.212;            % 机体转动惯量 (kg·m^2)
l_b_val = 0.043;            % 机体质心到转轴距离 (m)，注意符号
theta_b0_val = 0;           % 机体质心偏置角 (rad)

I_yaw_val = 0.276;          % 整体yaw轴转动惯量 (kg·m^2)

% ==================== 轮子参数 ====================
m_wl_val = 0.18677;         % 左轮质量 (kg)
m_wr_val = 0.18677;         % 右轮质量 (kg)
I_wl_val = 207897.79e-9;    % 左轮转动惯量 (kg·m^2)
I_wr_val = 207897.79e-9;    % 右轮转动惯量 (kg·m^2)

% ==================== 腿部参数 (可变腿长) ====================
l_l_val = 0.20;             % 左腿长度 (m)
l_r_val = 0.20;             % 右腿长度 (m)
m_l_val = 1.20485;          % 左腿质量 (kg)
m_r_val = 1.20485;          % 右腿质量 (kg)
I_l_val = 0.01298843908;    % 左腿转动惯量 (kg·m^2)
I_r_val = 0.01298843908;    % 右腿转动惯量 (kg·m^2)
l_l_d_val = 0.17523986903670066;  % 左腿质心到轮轴距离 (m)
l_r_d_val = 0.17523986903670066;  % 右腿质心到轮轴距离 (m)
theta_l0_val = 0;           % 左腿质心偏置角 (rad)
theta_r0_val = 0;           % 右腿质心偏置角 (rad)

fprintf('  物理参数设置完成\n\n');

%% ======================== Step 3: LQR权重矩阵 ========================

fprintf('Step 3: 设置LQR权重矩阵...\n');

% Q矩阵: 状态权重
%        X_b^h  dX_b^h  phi   dphi  theta_l dtheta_l theta_r dtheta_r theta_b dtheta_b
lqr_Q = diag([10,    1,      1,    1,    1000,    1,       1000,    1,       10000,  1]);

% R矩阵: 控制输入权重
%        T_wr   T_wl   T_r    T_l
lqr_R = diag([50,    50,    1,     1]);

fprintf('  Q矩阵 (状态权重):\n');
disp(lqr_Q);
fprintf('  R矩阵 (控制权重):\n');
disp(lqr_R);

%% ======================== Step 4: 数值代入并计算K矩阵 ========================

fprintf('Step 4: 计算LQR增益矩阵K...\n');

% 参数代入列表
param_syms = [g, R, R_w, m_b, I_b, l_b, theta_b0, I_yaw, ...
              m_wl, m_wr, I_wl, I_wr, ...
              m_l, m_r, I_l, I_r, l_l, l_r, l_l_d, l_r_d, theta_l0, theta_r0];

param_vals = [g_val, R_val, R_w_val, m_b_val, I_b_val, l_b_val, theta_b0_val, I_yaw_val, ...
              m_wl_val, m_wr_val, I_wl_val, I_wr_val, ...
              m_l_val, m_r_val, I_l_val, I_r_val, l_l_val, l_r_val, l_l_d_val, l_r_d_val, theta_l0_val, theta_r0_val];

% 代入数值
A_num = double(subs(A, param_syms, param_vals));
B_num = double(subs(B, param_syms, param_vals));

fprintf('  数值A矩阵:\n');
disp(A_num);
fprintf('  数值B矩阵:\n');
disp(B_num);

% 检查可控性
Co = ctrb(A_num, B_num);
rank_Co = rank(Co);
fprintf('  可控性矩阵秩: %d (需要 = 10)\n', rank_Co);
if rank_Co < 10
    fprintf('\n  ⚠ 系统不完全可控 (秩=%d < 10)\n', rank_Co);
    fprintf('  物理原因: X_b_h(水平位置) 和 phi(yaw角) 是积分器状态\n');
    fprintf('           机器人可以在任意位置/朝向平衡，这两个状态不影响动力学\n');
    fprintf('  解决方案: 这是正常的! LQR仍然可以计算可控子空间的增益\n\n');
end

% 计算LQR增益矩阵
try
    [K, S, e] = lqr(A_num, B_num, lqr_Q, lqr_R);
    fprintf('\n  LQR增益矩阵K (4x10):\n');
    disp(K);
    
    fprintf('  闭环特征值:\n');
    disp(e);
    
    % 检查稳定性 (允许0特征值对应不可控的积分器状态)
    stable_eigs = real(e) < 1e-6;  % 允许非常小的正实部
    if all(stable_eigs)
        fprintf('  ✓ 闭环系统稳定!\n\n');
    else
        warning('闭环系统不稳定!');
    end
catch ME
    fprintf('  LQR计算失败: %s\n', ME.message);
    K = [];
end

%% ======================== Step 5: 输出格式化结果 ========================

fprintf('========================================\n');
fprintf('格式化输出 (可直接复制到C代码):\n');
fprintf('========================================\n\n');

if ~isempty(K)
    % 输出K矩阵的C代码格式
    fprintf('// LQR增益矩阵 K[4][10]\n');
    fprintf('// 控制律: u = -K * x\n');
    fprintf('// u = [T_wr, T_wl, T_r, T_l]^T\n');
    fprintf('// x = [X, dX, phi, dphi, theta_l, dtheta_l, theta_r, dtheta_r, theta_b, dtheta_b]^T\n');
    fprintf('float K[4][10] = {\n');
    for i = 1:4
        fprintf('    {');
        for j = 1:10
            if j < 10
                fprintf('%.6gf, ', K(i,j));
            else
                fprintf('%.6gf', K(i,j));
            end
        end
        if i < 4
            fprintf('},\n');
        else
            fprintf('}\n');
        end
    end
    fprintf('};\n\n');
    
    % 输出单行格式 (便于复制)
    fprintf('// 单行格式:\n');
    K_str = sprintf([strjoin(repmat({'%.6g'},1,size(K,2)),',  ') '\n'], K.');
    fprintf('%s\n', K_str);
end

%% ======================== Step 6: 腿长拟合功能 ========================

fprintf('========================================\n');
fprintf('Step 6: 腿长拟合功能\n');
fprintf('========================================\n\n');

% 腿长数据集 (l_l, l_wl, l_bl, I_ll)
% 数据来源: 原HerKules脚本
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

enable_fitting = true;  % 设为true启用腿长拟合

if enable_fitting
    fprintf('正在计算不同腿长下的K矩阵...\n');
    
    % ========== 优化：将符号表达式转换为数值函数 ==========
    fprintf('  预编译数值函数 (这一步稍慢，但后续计算会很快)...\n');
    
    % 定义可变参数（腿长相关）和固定参数
    var_params = [l_l, l_r, l_l_d, l_r_d, I_l, I_r];
    
    % 先代入固定参数
    fixed_syms = [g, R, R_w, m_b, I_b, l_b, theta_b0, I_yaw, m_wl, m_wr, I_wl, I_wr, m_l, m_r, theta_l0, theta_r0];
    fixed_vals = [g_val, R_val, R_w_val, m_b_val, I_b_val, l_b_val, theta_b0_val, I_yaw_val, ...
                  m_wl_val, m_wr_val, I_wl_val, I_wr_val, m_l_val, m_r_val, theta_l0_val, theta_r0_val];
    
    A_partial = subs(A, fixed_syms, fixed_vals);
    B_partial = subs(B, fixed_syms, fixed_vals);
    
    % 转换为数值函数 (参数顺序: l_l, l_r, l_l_d, l_r_d, I_l, I_r)
    A_func = matlabFunction(A_partial, 'Vars', var_params);
    B_func = matlabFunction(B_partial, 'Vars', var_params);
    
    fprintf('  数值函数编译完成!\n');
    
    num_legs = size(Leg_data, 1);
    sample_size = num_legs^2;
    
    % K矩阵有4x10=40个元素
    K_sample = zeros(sample_size, 3, 40);  % [l_l, l_r, K_ij]
    
    idx = 0;
    tic_fit = tic;
    for i = 1:num_legs
        for j = 1:num_legs
            idx = idx + 1;
            
            % 提取左腿数据
            l_l_fit = Leg_data(i, 1);
            l_l_d_fit = Leg_data(i, 2);
            I_l_fit = Leg_data(i, 4);
            
            % 提取右腿数据
            l_r_fit = Leg_data(j, 1);
            l_r_d_fit = Leg_data(j, 2);
            I_r_fit = Leg_data(j, 4);
            
            % 使用预编译的数值函数计算 (非常快!)
            A_fit = A_func(l_l_fit, l_r_fit, l_l_d_fit, l_r_d_fit, I_l_fit, I_r_fit);
            B_fit = B_func(l_l_fit, l_r_fit, l_l_d_fit, l_r_d_fit, I_l_fit, I_r_fit);
            
            % 计算LQR
            try
                K_fit = lqr(A_fit, B_fit, lqr_Q, lqr_R);
                
                % 存储结果
                for k = 1:40
                    K_sample(idx, 1, k) = l_l_fit;
                    K_sample(idx, 2, k) = l_r_fit;
                    row = ceil(k/10);
                    col = mod(k-1, 10) + 1;
                    K_sample(idx, 3, k) = K_fit(row, col);
                end
            catch
                warning('LQR计算失败: l_l=%.2f, l_r=%.2f', l_l_fit, l_r_fit);
            end
        end
        % 每完成一行显示进度
        if mod(i, 5) == 0
            fprintf('  进度: %d/%d 行 (%.1f秒)\n', i, num_legs, toc(tic_fit));
        end
    end
    fprintf('  全部 %d 个样本计算完成! 耗时: %.2f秒\n', sample_size, toc(tic_fit));
    
    % 多项式拟合
    fprintf('\n正在进行多项式拟合...\n');
    K_Fit_Coefficients = zeros(40, 6);
    
    for n = 1:40
        try
            K_Surface_Fit = fit([K_sample(:,1,n), K_sample(:,2,n)], K_sample(:,3,n), 'poly22');
            K_Fit_Coefficients(n,:) = coeffvalues(K_Surface_Fit);
        catch
            warning('拟合失败: K元素 %d', n);
        end
    end
    
    fprintf('\n拟合多项式表达式: p(l_l,l_r) = p00 + p10*l_l + p01*l_r + p20*l_l^2 + p11*l_l*l_r + p02*l_r^2\n\n');
    
    fprintf('// 拟合系数 K_Fit_Coefficients[40][6]\n');
    fprintf('// 行: K矩阵元素 (K_11, K_12, ..., K_1_10, K_21, ..., K_4_10)\n');
    fprintf('// 列: [p00, p10, p01, p20, p11, p02]\n');
    K_Fit_str = sprintf([strjoin(repmat({'%.6g'},1,6),',  ') '\n'], K_Fit_Coefficients.');
    fprintf('%s\n', K_Fit_str);
    
    % 保存结果
    save('lqr_fitting_results.mat', 'K_sample', 'K_Fit_Coefficients', 'Leg_data');
    fprintf('拟合结果已保存到 lqr_fitting_results.mat\n');
end

%% ======================== 保存结果 ========================

save('lqr_results.mat', 'A_num', 'B_num', 'K', 'lqr_Q', 'lqr_R', 'e');
fprintf('\nLQR结果已保存到 lqr_results.mat\n');

toc
fprintf('\n计算完成!\n');
