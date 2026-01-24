# 轮腿机器人动力学推导与LQR控制器设计

## 项目概述

本项目完成了轮腿机器人的完整动力学推导、线性化分析以及LQR控制器设计。

### 机器人模型

- **结构**：双轮差速驱动 + 双腿（左右独立）
- **自由度**：5个广义坐标
- **控制输入**：4个力矩（左右轮电机 + 左右髋关节电机）

---

## 文件结构

```
2026公式推导/
├── README.md                    # 本文件
├── 推导-v2.md                   # 完整动力学推导文档
│
├── step1_define_and_derive.m    # Step 1: 定义符号变量，建立动力学方程
├── step1_results.mat            # Step 1 输出结果
│
├── simplify_dynamics.m          # 动力学方程化简（旧版）
├── simplify_dynamics_v2.m       # 动力学方程化简（新坐标系）
├── dynamics_v2.mat              # 化简后的动力学方程
│
├── apply_kinematics_v2.m        # Step 2: 代入运动学约束，提取M、B、g矩阵
├── dynamics_new_coords.mat      # M、B、g矩阵结果
│
├── linearize_system.m           # 线性化脚本（初版）
├── linearize_system_v2.m        # 线性化脚本（优化版，已合并至compute_lqr_v2.m）
├── linearized_system.mat        # 线性化A、B矩阵结果
│
├── compute_lqr.m                # LQR控制器计算（旧版，不含平衡点求解）
├── compute_lqr_v2.m             # ★ 主脚本：平衡点求解 + 线性化 + LQR控制器
├── lqr_results.mat              # LQR计算结果（单腿长）
└── lqr_fitting_results.mat      # LQR拟合结果（变腿长）
```

---

## 工作流程

### Step 1: 动力学方程推导

**文件**: `推导-v2.md`, `simplify_dynamics_v2.m`

1. 定义所有符号变量（质量、惯量、几何参数等）
2. 基于牛顿-欧拉方法建立各刚体的动力学方程
3. 利用牛顿第三定律消去内力
4. 得到5个最终动力学方程

**广义坐标** $\mathbf{q}$ (5维):
$$\mathbf{q} = [X_b^h, \phi, \theta_l, \theta_r, \theta_b]^T$$

- $X_b^h$ — 机体水平位置
- $\phi$ — 偏航角 (Yaw)
- $\theta_l$ — 左腿角度
- $\theta_r$ — 右腿角度  
- $\theta_b$ — 机体俯仰角 (Pitch)

### Step 2: 代入运动学约束

**文件**: `apply_kinematics_v2.m`

将运动学约束代入动力学方程，整理为矩阵形式：
$$M(\mathbf{q}) \cdot \ddot{\mathbf{q}} = B(\mathbf{q}) \cdot \mathbf{u} + \mathbf{g}(\mathbf{q}, \dot{\mathbf{q}})$$

**输出**:
- $M$ — 5×5 质量矩阵
- $B$ — 5×4 控制矩阵
- $\mathbf{g}$ — 5×1 重力+离心力项

### Step 3: 线性化

**文件**: `linearize_system_v2.m`（已合并至 `compute_lqr_v2.m`）

在平衡点处线性化，得到状态空间模型：
$$\dot{\mathbf{x}} = A \mathbf{x} + B_c \mathbf{u}$$

**注意**：平衡点不再假设所有角度为零，而是通过求解 $\mathbf{g}(\mathbf{q}) = 0$ 得到。详见 Step 4。

**状态向量** $\mathbf{x}$ (10维):
$$\mathbf{x} = [X_b^h, \dot{X}_b^h, \phi, \dot{\phi}, \theta_l, \dot{\theta}_l, \theta_r, \dot{\theta}_r, \theta_b, \dot{\theta}_b]^T$$

**控制向量** $\mathbf{u}$ (4维):
$$\mathbf{u} = [T_{r\to b}, T_{l \to b}, T_{wr \to r}, T_{wl \to l}]^T$$

- $T_{wr}$ — 右轮电机力矩
- $T_{wl}$ — 左轮电机力矩
- $T_r$ — 右髋关节电机力矩
- $T_l$ — 左髋关节电机力矩

### Step 4: 平衡点求解与LQR控制器设计

**文件**: `compute_lqr_v2.m` ★

#### 4.1 平衡点求解

平衡点条件：所有速度、加速度为零，控制输入为零时，重力项也必须为零：
$$\mathbf{g}(\mathbf{q}_{eq}) = 0$$

由于 $\mathbf{g}$ 向量中的方程是独立的 $\sin(\theta)$ 形式，可以符号化直接求解：
- $\theta_b^{eq}$ 由 $g_2 = 0$ 求得
- $\theta_r^{eq}$ 由 $g_3 = 0$ 求得  
- $\theta_l^{eq}$ 由 $g_4 = 0$ 求得

平衡点偏移量（equilibrium offset）会随腿长变化，需要输出到控制器中使用。

#### 4.2 LQR 控制律

计算LQR增益矩阵 $K$，使得控制律：
$$\mathbf{u} = -K (\mathbf{x} - \mathbf{x}_{eq})$$

其中 $\mathbf{x}_{eq}$ 是平衡点状态（速度为零，角度为平衡点值）。

最小化代价函数：
$$J = \int_0^\infty (\mathbf{x}^T Q \mathbf{x} + \mathbf{u}^T R \mathbf{u}) dt$$

---

## 快速开始

### 单腿长LQR计算

1. 打开 `compute_lqr_v2.m`
2. 修改物理参数（Part 5 区域）
3. 修改Q、R权重矩阵（Part 11 区域）
4. 确保 `enable_fitting = false`
5. 运行脚本

```matlab
matlab -batch "run('compute_lqr_v2.m')"
```

输出：
- 控制台打印K矩阵和**平衡点偏移量**（可直接复制到C代码）
- `lqr_results.mat` 保存完整结果

### 变腿长拟合

1. 设置 `enable_fitting = true`
2. 修改 `Leg_data` 数据集（不同腿长对应的参数）
3. 运行脚本

输出：
- `lqr_fitting_results.mat` 包含拟合系数
- K矩阵拟合公式：$K_{ij}(l_l, l_r) = p_{00} + p_{10} l_l + p_{01} l_r + p_{20} l_l^2 + p_{11} l_l l_r + p_{02} l_r^2$
- 平衡点偏移拟合公式：$\theta_{eq}(l_l, l_r)$ 同样用二次多项式拟合

---

## 参数说明

### 机体参数

| 参数 | 符号 | 单位 | 示例值 |
|:---|:---:|:---:|:---:|
| 轮子半径 | $R$ | m | 0.055 |
| 轮距/2 | $R_w$ | m | 0.225 |
| 机体质量 | $m_b$ | kg | 15.564 |
| 机体转动惯量 | $I_b$ | kg·m² | 0.212 |
| 机体质心距离 | $l_b$ | m | 0.043 |
| yaw轴转动惯量 | $I_{yaw}$ | kg·m² | 0.276 |

### 轮子参数

| 参数 | 符号 | 单位 | 示例值 |
|:---|:---:|:---:|:---:|
| 轮子质量 | $m_w$ | kg | 0.187 |
| 轮子转动惯量 | $I_w$ | kg·m² | 2.08e-4 |

### 腿部参数（可变）

| 参数 | 符号 | 单位 | 示例值 |
|:---|:---:|:---:|:---:|
| 腿长 | $l_{leg}$ | m | 0.13~0.40 |
| 腿质量 | $m_{leg}$ | kg | 1.205 |
| 腿转动惯量 | $I_{leg}$ | kg·m² | 0.013~0.017 |
| 腿质心距离 | $l_{leg}^d$ | m | 随腿长变化 |

### LQR权重矩阵

**Q矩阵**（状态权重）:
```
        X    dX   phi  dphi  θ_l  dθ_l  θ_r  dθ_r  θ_b  dθ_b
Q = diag[10,  1,   1,   1,   1000,  1,  1000,  1,  10000,  1]
```

**R矩阵**（控制权重）:
```
        T_wr  T_wl  T_r  T_l
R = diag[50,   50,   1,   1]
```

---

## C代码使用示例

### 固定腿长

```c
// LQR增益矩阵 K[4][10]
// 控制律: u = -K * (x - x_eq)
float K[4][10] = {
    {-0.278925f, -1.03799f, -0.0997239f, -0.129885f, -0.11638f, -0.057609f, 0.253744f, -0.0369642f, -4.76675f, -0.39819f},
    // ... (从MATLAB输出复制)
};

// 平衡点偏移量
float theta_l_eq = 0.0f;  // 从MATLAB输出复制
float theta_r_eq = 0.0f;  // 从MATLAB输出复制
float theta_b_eq = 0.0f;  // 从MATLAB输出复制
float X_eq[10] = {0, 0, 0, 0, theta_l_eq, 0, theta_r_eq, 0, theta_b_eq, 0};

// 状态向量
float x[10] = {X, dX, phi, dphi, theta_l, dtheta_l, theta_r, dtheta_r, theta_b, dtheta_b};

// 计算控制力矩: u = -K * (x - x_eq)
float u[4];
for (int i = 0; i < 4; i++) {
    u[i] = 0;
    for (int j = 0; j < 10; j++) {
        u[i] -= K[i][j] * (x[j] - X_eq[j]);
    }
}
// u[0] = T_r_to_b, u[1] = T_l_to_b, u[2] = T_wr_to_r, u[3] = T_wl_to_l
```

### 变腿长（使用拟合）

```c
// K矩阵拟合系数 K_coef[40][6]
// 第n个K元素: K_n = p00 + p10*l_l + p01*l_r + p20*l_l^2 + p11*l_l*l_r + p02*l_r^2
float K_coef[40][6] = { /* 从MATLAB输出复制 */ };

// 平衡点偏移拟合系数 Offset_coef[3][6]
// [0]: theta_l_eq, [1]: theta_r_eq, [2]: theta_b_eq
float Offset_coef[3][6] = { /* 从MATLAB输出复制 */ };

// 根据腿长计算平衡点偏移
float theta_l_eq = Offset_coef[0][0] + Offset_coef[0][1]*l_l + Offset_coef[0][2]*l_r 
                 + Offset_coef[0][3]*l_l*l_l + Offset_coef[0][4]*l_l*l_r + Offset_coef[0][5]*l_r*l_r;
float theta_r_eq = Offset_coef[1][0] + Offset_coef[1][1]*l_l + Offset_coef[1][2]*l_r 
                 + Offset_coef[1][3]*l_l*l_l + Offset_coef[1][4]*l_l*l_r + Offset_coef[1][5]*l_r*l_r;
float theta_b_eq = Offset_coef[2][0] + Offset_coef[2][1]*l_l + Offset_coef[2][2]*l_r 
                 + Offset_coef[2][3]*l_l*l_l + Offset_coef[2][4]*l_l*l_r + Offset_coef[2][5]*l_r*l_r;

// 根据腿长计算K矩阵
float K[4][10];
for (int n = 0; n < 40; n++) {
    int row = n / 10;
    int col = n % 10;
    K[row][col] = K_coef[n][0] 
                + K_coef[n][1] * l_l 
                + K_coef[n][2] * l_r
                + K_coef[n][3] * l_l * l_l 
                + K_coef[n][4] * l_l * l_r 
                + K_coef[n][5] * l_r * l_r;
}
```

---

## 注意事项

1. **可控性**: 系统可控性矩阵秩为8（而非10），因为位置和yaw角是可积分状态。这是正常的，LQR仍能稳定系统。

2. **平衡点**: 平衡点通过求解 $\mathbf{g} = 0$ 得到，不再假设所有角度为零。控制律使用 $\mathbf{u} = -K(\mathbf{x} - \mathbf{x}_{eq})$ 形式，需要减去平衡点偏移。

3. **参数一致性**: 确保MATLAB中的参数与实际机器人一致，特别是：
   - 腿部惯量随腿长变化
   - 腿质心位置随腿长变化

4. **符号约定**: 
   - 角度以Z轴负方向为0参考
   - 力矩正方向遵循右手定则

---

## 参考文献

- 2023上交轮腿开源: https://bbs.robomaster.com/forum.php?mod=viewthread&tid=22756
- 哈工程轮腿开源: https://zhuanlan.zhihu.com/p/563048952

---

## 更新日志

- **2026/01/15**: 新增 `compute_lqr_v2.m`，整合平衡点求解、线性化和LQR计算；支持平衡点偏移输出和拟合
- **2026/01/09**: 完成新动力学模型推导，实现LQR控制器计算和腿长拟合功能