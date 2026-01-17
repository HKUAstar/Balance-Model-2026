% compute_lqr.m
% 基于线性化状态空间模型计算LQR控制器增益矩阵
% 
% 依赖文件: linearized_system.mat (由 linearize_system_v2.m 生成)
%
% ========================================================================
%                          变量定义 (对应推导文档 §4.3)
% ========================================================================
%
% 状态向量 X (10维):
%   X = [X_b^h; V_b^h; phi; dphi; theta_l; dtheta_l; theta_r; dtheta_r; theta_b; dtheta_b]
%
%   序号   符号          物理意义                        单位
%   ─────────────────────────────────────────────────────────────────────
%    1    X_b^h        机体水平位置                      m
%    2    V_b^h        机体水平速度 (= dX_b^h/dt)        m/s
%    3    phi          偏航角                            rad
%    4    dphi         偏航角速度                        rad/s
%    5    theta_l      左腿与Z轴负方向夹角               rad
%    6    dtheta_l     左腿角速度                        rad/s
%    7    theta_r      右腿与Z轴负方向夹角               rad
%    8    dtheta_r     右腿角速度                        rad/s
%    9    theta_b      机体俯仰角                        rad
%   10    dtheta_b     机体俯仰角速度                    rad/s
%
% 控制向量 u (4维):
%   u = [T_{r→b}; T_{l→b}; T_{wr→r}; T_{wl→l}]
%
%   序号   符号          物理意义                        执行器
%   ─────────────────────────────────────────────────────────────────────
%    1    T_{r→b}      右腿对机体的扭矩                  右髋关节电机
%    2    T_{l→b}      左腿对机体的扭矩                  左髋关节电机
%    3    T_{wr→r}     右轮对右腿的扭矩                  右轮电机
%    4    T_{wl→l}     左轮对左腿的扭矩                  左轮电机
%
% 扭矩符号约定: T_{A→B} 表示物体A对物体B施加的扭矩
%
% LQR控制律:
%   u = -K * X
%
% 其中 K 为 4×10 增益矩阵
%
% ========================================================================
% 作者: 基于2026公式推导
% 日期: 2026/01/13
% ========================================================================

clear all; clc;
tic

%% ======================== Step 0: 加载线性化系统 ========================

fprintf('========================================\n');
fprintf('轮腿机器人LQR控制器计算\n');
fprintf('========================================\n\n');

fprintf('Step 0: 加载线性化状态空间模型...\n');

% 检查文件是否存在
if ~exist('linearized_system.mat', 'file')
    error('未找到 linearized_system.mat! 请先运行 linearize_system_v2.m');
end

load('linearized_system.mat', 'A_func', 'B_func', 'A_sym', 'B_state_sym', 'param_list');

fprintf('  ✓ 线性化系统加载成功\n');

% 获取函数句柄的输出维度 (22 个参数)
% param_list = [m_b, m_l, m_r, m_wl, m_wr, I_b, I_l, I_r, I_wl, I_wr, I_yaw, l_l, l_r, l_l_d, l_r_d, l_b, R, R_w, g, theta_l0, theta_r0, theta_b0]
test_params = [1, 1, 1, 0.2, 0.2, 0.1, 0.03, 0.03, 0.0002, 0.0002, 0.4, 0.2, 0.2, 0.084, 0.084, 0.05, 0.055, 0.1, -9.81, 0.1, 0.1, 0.05];
A_test = A_func(test_params);
B_test = B_func(test_params);
fprintf('  状态维度: %d\n', size(A_test, 1));
fprintf('  控制维度: %d\n', size(B_test, 2));

%% ======================== Step 1: 定义物理参数 ========================

fprintf('\nStep 1: 定义机器人物理参数...\n');

% linearize_system_v2.m 使用左右腿/轮参数分开模型，参数列表为:
%   [m_b, m_l, m_r, m_wl, m_wr, I_b, I_l, I_r, I_wl, I_wr, I_yaw, l_l, l_r, l_l_d, l_r_d, l_b, R, R_w, g, theta_l0, theta_r0, theta_b0]
% 其中左右腿/轮参数可以独立设置（支持不对称）

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
theta_b0 = 0;               % 质心偏移角度，单位：弧度

% ==================== 轮子参数 (分别定义左右) ====================
m_wl_val = 0.19;            % 左轮质量 (kg)
m_wr_val = 0.19;            % 右轮质量 (kg)
I_wl_val = 0.000207897;     % 左轮转动惯量 (kg·m²)
I_wr_val = 0.000207897;     % 右轮转动惯量 (kg·m²)

% ==================== 腿部参数 (分别定义左右, 默认腿长 0.20m) ====================
l_l_val = 0.20;             % 左腿长度 (m)
l_r_val = 0.20;             % 右腿长度 (m)
m_l_val = 1.62;             % 左腿质量 (kg)
m_r_val = 1.62;             % 右腿质量 (kg)
I_l_val = 0.0339;           % 左腿转动惯量 (kg·m²)
I_r_val = 0.0339;           % 右腿转动惯量 (kg·m²)
l_l_d_val = 0.084685385397954;  % 左腿质心到轮轴距离 (m)
l_r_d_val = 0.084685385397954;  % 右腿质心到轮轴距离 (m)
theta_l0 = 0.8152113973036387;  % 左腿偏移角度，单位：弧度
theta_r0 = 0.8152113973036387;  % 右腿偏移角度，单位：弧度

fprintf('  ✓ 物理参数设置完成\n');

%% ======================== Step 2: 数值代入 ========================

fprintf('\nStep 2: 代入数值参数...\n');

% 参数值向量 (顺序与 linearize_system_v2.m 中 param_list 一致)
%   param_list = [m_b, m_l, m_r, m_wl, m_wr, I_b, I_l, I_r, I_wl, I_wr, I_yaw, l_l, l_r, l_l_d, l_r_d, l_b, R, R_w, g, theta_l0, theta_r0, theta_b0]
param_vals = [m_b_val, m_l_val, m_r_val, m_wl_val, m_wr_val, I_b_val, I_l_val, I_r_val, I_wl_val, I_wr_val, I_yaw_val, ...
              l_l_val, l_r_val, l_l_d_val, l_r_d_val, l_b_val, R_val, R_w_val, g_val, theta_l0, theta_r0, theta_b0];

% 使用函数句柄计算数值矩阵
A_num = A_func(param_vals);
B_num = B_func(param_vals);

fprintf('  ✓ 数值代入完成\n');

% 显示矩阵
fprintf('\n  数值A矩阵 (10×10):\n');
disp(A_num);
fprintf('  数值B矩阵 (10×4):\n');
disp(B_num);

%% ======================== Step 3: 检查可控性 ========================

fprintf('Step 3: 检查系统可控性...\n');

Co = ctrb(A_num, B_num);
rank_Co = rank(Co);
fprintf('  可控性矩阵秩: %d (系统维度: 10)\n', rank_Co);

if rank_Co < 10
    fprintf('\n  ⚠ 系统不完全可控 (秩=%d < 10)\n', rank_Co);
    fprintf('  物理原因: X_b^h(水平位置) 和 phi(yaw角) 是积分器状态\n');
    fprintf('           机器人可以在任意位置/朝向平衡，这两个状态不影响动力学\n');
    fprintf('  解决方案: 这是正常的! LQR仍然可以计算可控子空间的增益\n\n');
else
    fprintf('  ✓ 系统完全可控\n\n');
end

%% ======================== Step 4: 设置LQR权重 ========================

fprintf('Step 4: 设置LQR权重矩阵...\n');

% Q矩阵: 状态权重
% 状态: [X_b^h, V_b^h, phi, dphi, theta_l, dtheta_l, theta_r, dtheta_r, theta_b, dtheta_b]
%        位置    速度           偏航  偏航速  左腿角   左腿速   右腿角   右腿速    俯仰角   俯仰速
lqr_Q = diag([100,    1,      4000,    1,      1000,     10,       1000,     10,       40000,    1]);

% R矩阵: 控制输入权重
% 控制: [T_{r→b}, T_{l→b}, T_{wr→r}, T_{wl→l}]
%        右髋扭矩   左髋扭矩   右轮扭矩   左轮扭矩
lqr_R = diag([1,       1,        10,        10]);

fprintf('  Q矩阵 (状态权重):\n');
fprintf('         X_b^h  V_b^h  phi   dphi  θ_l   dθ_l  θ_r   dθ_r  θ_b   dθ_b\n');
disp(lqr_Q);

fprintf('  R矩阵 (控制权重):\n');
fprintf('         T_{r→b}  T_{l→b}  T_{wr→r}  T_{wl→l}\n');
disp(lqr_R);

%% ======================== Step 5: 计算LQR增益 ========================

fprintf('Step 5: 计算LQR增益矩阵K...\n');

try
    [K, S, e] = lqr(A_num, B_num, lqr_Q, lqr_R);
    
    fprintf('\n  ✓ LQR增益矩阵K (4×10):\n');
    disp(K);
    
    fprintf('  闭环特征值:\n');
    disp(e);
    
    % 检查稳定性
    stable_eigs = real(e) < 1e-6;
    if all(stable_eigs)
        fprintf('  ✓ 闭环系统稳定!\n\n');
    else
        warning('闭环系统不稳定!');
    end
catch ME
    fprintf('  ✗ LQR计算失败: %s\n', ME.message);
    K = [];
end

%% ======================== Step 6: 格式化输出 ========================

fprintf('========================================\n');
fprintf('格式化输出 (可直接复制到C代码)\n');
fprintf('========================================\n\n');

if ~isempty(K)
    % K矩阵行列含义
    fprintf('// ═══════════════════════════════════════════════════════════════════════\n');
    fprintf('// LQR增益矩阵 K[4][10]\n');
    fprintf('// ═══════════════════════════════════════════════════════════════════════\n');
    fprintf('// 控制律: u = -K * X\n');
    fprintf('//\n');
    fprintf('// 状态向量 X (列向量 10×1):\n');
    fprintf('//   X[0] = X_b^h     机体水平位置 (m)\n');
    fprintf('//   X[1] = V_b^h     机体水平速度 (m/s)\n');
    fprintf('//   X[2] = phi       偏航角 (rad)\n');
    fprintf('//   X[3] = dphi      偏航角速度 (rad/s)\n');
    fprintf('//   X[4] = theta_l   左腿角 (rad)\n');
    fprintf('//   X[5] = dtheta_l  左腿角速度 (rad/s)\n');
    fprintf('//   X[6] = theta_r   右腿角 (rad)\n');
    fprintf('//   X[7] = dtheta_r  右腿角速度 (rad/s)\n');
    fprintf('//   X[8] = theta_b   机体俯仰角 (rad)\n');
    fprintf('//   X[9] = dtheta_b  机体俯仰角速度 (rad/s)\n');
    fprintf('//\n');
    fprintf('// 控制向量 u (列向量 4×1):\n');
    fprintf('//   u[0] = T_r_to_b   右髋扭矩 (右腿→机体) (Nm)\n');
    fprintf('//   u[1] = T_l_to_b   左髋扭矩 (左腿→机体) (Nm)\n');
    fprintf('//   u[2] = T_wr_to_r  右轮扭矩 (右轮→右腿) (Nm)\n');
    fprintf('//   u[3] = T_wl_to_l  左轮扭矩 (左轮→左腿) (Nm)\n');
    fprintf('//\n');
    fprintf('// K矩阵含义:\n');
    fprintf('//   K[i][j] 表示控制输入 u[i] 对状态 X[j] 的反馈增益\n');
    fprintf('//   K[0][*]: 右髋扭矩对各状态的增益\n');
    fprintf('//   K[1][*]: 左髋扭矩对各状态的增益\n');
    fprintf('//   K[2][*]: 右轮扭矩对各状态的增益\n');
    fprintf('//   K[3][*]: 左轮扭矩对各状态的增益\n');
    fprintf('// ═══════════════════════════════════════════════════════════════════════\n\n');
    
    fprintf('float K[4][10] = {\n');
    control_names = {'T_r_to_b', 'T_l_to_b', 'T_wr_to_r', 'T_wl_to_l'};
    for i = 1:4
        fprintf('    {%11.6ff, %11.6ff, %11.6ff, %11.6ff, %11.6ff, ', K(i,1), K(i,2), K(i,3), K(i,4), K(i,5));
        fprintf('%11.6ff, %11.6ff, %11.6ff, %11.6ff, %11.6ff}', K(i,6), K(i,7), K(i,8), K(i,9), K(i,10));
        if i < 4
            fprintf(',  // %s\n', control_names{i});
        else
            fprintf('   // %s\n', control_names{i});
        end
    end
    fprintf('};\n\n');
    
    % 单行格式
    fprintf('// 单行格式 (每行对应一个控制输入):\n');
    for i = 1:4
        fprintf('// K[%d] (%s): ', i-1, control_names{i});
        fprintf('%.6g, ', K(i,1:9));
        fprintf('%.6g\n', K(i,10));
    end
    fprintf('\n');
end

%% ======================== Step 7: 腿长拟合功能 ========================

fprintf('========================================\n');
fprintf('Step 7: 腿长拟合功能\n');
fprintf('========================================\n\n');

% 腿长参数查找表
% 格式: [腿长(m), 质心偏离角度(rad), 质心到髋关节距离(m), 转动惯量(kg·m²)]
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

enable_fitting = true;  % 设为 false 跳过腿长拟合

if enable_fitting
    fprintf('正在计算不同腿长下的K矩阵...\n');
    fprintf('  注意: 使用左右腿参数分开模型，支持左右腿长不同\n\n');
    
    % ========== 计算采样点 (二维网格) ==========
    num_legs = size(Leg_data, 1);
    sample_size_2d = num_legs^2;
    
    % K矩阵 4×10 = 40 个元素
    % 二维拟合: [l_l, l_r, K矩阵的40个元素]
    K_sample_2d = zeros(sample_size_2d, 44);  % [l_l, l_r, K矩阵的40个元素]
    
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
            
            % 构建参数向量
            % param_list = [m_b, m_l, m_r, m_wl, m_wr, I_b, I_l, I_r, I_wl, I_wr, I_yaw, l_l, l_r, l_l_d, l_r_d, l_b, R, R_w, g, theta_l0, theta_r0, theta_b0]
            param_fit = [m_b_val, m_l_val, m_r_val, m_wl_val, m_wr_val, I_b_val, I_l_fit, I_r_fit, I_wl_val, I_wr_val, I_yaw_val, ...
                         l_l_fit, l_r_fit, l_l_d_fit, l_r_d_fit, l_b_val, R_val, R_w_val, g_val, theta_l0_fit, theta_r0_fit, theta_b0];
            
            % 计算数值矩阵
            A_fit = A_func(param_fit);
            B_fit = B_func(param_fit);
            
            % 计算LQR
            try
                K_fit = lqr(A_fit, B_fit, lqr_Q, lqr_R);
                
                % 存储结果
                K_sample_2d(idx, 1) = l_l_fit;
                K_sample_2d(idx, 2) = l_r_fit;
                K_sample_2d(idx, 3:42) = K_fit(:)';  % 按行展开
            catch
                warning('LQR计算失败: l_l=%.2f, l_r=%.2f', l_l_fit, l_r_fit);
            end
            
            % 显示进度
            if mod(idx, 49) == 0
                fprintf('  进度: %d/%d (%.1f秒)\n', idx, sample_size_2d, toc(tic_fit));
            end
        end
    end
    
    fprintf('  ✓ %d 个样本计算完成! 耗时: %.2f秒\n', sample_size_2d, toc(tic_fit));
    
    % ========== 二维多项式拟合 ==========
    fprintf('\n正在进行二维多项式拟合...\n');
    
    % 拟合多项式: K_ij(l_l, l_r) = p00 + p10*l_l + p01*l_r + p20*l_l^2 + p11*l_l*l_r + p02*l_r^2
    K_Fit_Coefficients = zeros(40, 6);
    
    l_l_samples = K_sample_2d(:, 1);
    l_r_samples = K_sample_2d(:, 2);
    
    for n = 1:40
        K_values = K_sample_2d(:, n+2);
        try
            % 二维二次多项式拟合
            K_Surface_Fit = fit([l_l_samples, l_r_samples], K_values, 'poly22');
            coeffs = coeffvalues(K_Surface_Fit);
            K_Fit_Coefficients(n, :) = coeffs;  % [p00, p10, p01, p20, p11, p02]
        catch
            warning('二维拟合失败: 元素 %d', n);
        end
    end
    
    fprintf('  ✓ 拟合完成\n\n');
    
    % ========== 输出拟合系数 ==========
    fprintf('// ═══════════════════════════════════════════════════════════════════════\n');
    fprintf('// 腿长拟合系数 K_Fit_Coefficients[40][6] (左右腿可不同)\n');
    fprintf('// ═══════════════════════════════════════════════════════════════════════\n');
    fprintf('// 拟合多项式 (二维):\n');
    fprintf('//   K_ij(l_l, l_r) = p00 + p10*l_l + p01*l_r + p20*l_l^2 + p11*l_l*l_r + p02*l_r^2\n');
    fprintf('//\n');
    fprintf('// K元素排列 (按行优先):\n');
    fprintf('//   n=0~9:   K[0][0~9] → T_{r→b} 对各状态的增益\n');
    fprintf('//   n=10~19: K[1][0~9] → T_{l→b} 对各状态的增益\n');
    fprintf('//   n=20~29: K[2][0~9] → T_{wr→r} 对各状态的增益\n');
    fprintf('//   n=30~39: K[3][0~9] → T_{wl→l} 对各状态的增益\n');
    fprintf('//\n');
    fprintf('// 系数顺序: [p00, p10, p01, p20, p11, p02]\n');
    fprintf('// ═══════════════════════════════════════════════════════════════════════\n\n');
    
    fprintf('float K_Fit_Coefficients[40][6] = {\n');
    for n = 1:40
        row = ceil(n/10) - 1;  % 0-indexed
        col = mod(n-1, 10);    % 0-indexed
        fprintf('    {');
        for c = 1:6
            if c < 6
                fprintf('%12.6gf, ', K_Fit_Coefficients(n,c));
            else
                fprintf('%12.6gf', K_Fit_Coefficients(n,c));
            end
        end
        if n < 40
            fprintf('},  // K[%d][%d]\n', row, col);
        else
            fprintf('}   // K[%d][%d]\n', row, col);
        end
    end
    fprintf('};\n\n');
    
    % 保存拟合结果
    save('lqr_fitting_results.mat', 'K_sample_2d', 'K_Fit_Coefficients', 'Leg_data');
    fprintf('拟合结果已保存到 lqr_fitting_results.mat (左右腿独立参数版本)\n');
end

%% ======================== Step 8: 保存结果 ========================

fprintf('\n========================================\n');
fprintf('保存结果\n');
fprintf('========================================\n\n');

save('lqr_results.mat', 'A_num', 'B_num', 'K', 'lqr_Q', 'lqr_R', 'e');
fprintf('✓ LQR结果已保存到 lqr_results.mat\n');

elapsed_time = toc;
fprintf('\n计算完成! 总耗时: %.2f秒\n', elapsed_time);
