# BATSRUS 核心求解器深度分析报告

> **Block Adaptive Tree Solar-wind Roe-type Upwind Scheme**
>
> 分析日期：2026-05-01 | 分析师：batsrus_solver（deepseek_pro 分身）
>
> 源仓库：https://github.com/SWMFsoftware/BATSRUS

---

## 目录

1. [总体架构](#1-总体架构)
2. [时间推进方案](#2-时间推进方案)
3. [Roe 近似黎曼求解器](#3-roe-近似黎曼求解器)
4. [高阶重构（MUSCL/TVD/WENO）](#4-高阶重构muscltvdweno)
5. [时间步约束（CFL 条件）](#5-时间步约束cfl-条件)
6. [守恒性与精度](#6-守恒性与精度)
7. [八波方案的 divB 处理](#7-八波方案的-divb-处理)
8. [特殊数值技术](#8-特殊数值技术)

---

## 1. 总体架构

### 1.1 求解器全称与历史

**BATSRUS** 全称 **B**lock **A**daptive **T**ree **S**olar-wind **R**oe-type **U**pwind **S**cheme，由密歇根大学开发（Powell et al., J. Comp. Phys., 154, 284, 1999），是 Space Weather Modeling Framework (SWMF) 的核心 MHD 模块。

### 1.2 核心数据结构

```
BATSRUS 采用自相似块结构：
┌─────────────────────────────────────┐
│  网格：Block-Adaptive Tree (BAT)    │
│  - 每个块包含 nI×nJ×nK 个物理单元   │
│  - 自相似：所有块尺寸相同            │
│  - 分层嵌套：从根块逐层细化          │
│  - 并行：MPI + OpenMP + OpenACC 混合 │
└─────────────────────────────────────┘
```

### 1.3 顶层调用流程

```
BATS_setup          → 初始网格建立、初条件设置
    ↓
BATS_init_session   → 初始化物理模块、选择步进策略
    ↓
BATS_advance        → 单时间步推进（被 SWMF 或 standalone 主循环调用）
    ├── set_global_timestep       → 全局时间步计算
    ├── load_balance_blocks       → 块类型选择 + 负载均衡
    ├── advance_explicit           → 显式推进（默认主路径）
    │   └── STAGELOOP (iStage=1..nStage)
    │       ├── calc_face_value   → MUSCL 重构面值
    │       ├── calc_face_flux    → 面通量计算（Roe/HLLD/LF等）
    │       ├── calc_source       → 源项计算（含 divB 源项）
    │       ├── calc_timestep     → 局部时间步计算
    │       └── update_state      → 状态更新
    ├── advance_semi_impl         → 半隐式处理（可选）
    └── prepare_amr + do_amr      → 网格自适应（按频率触发）
```

---

## 2. 时间推进方案

### 2.1 多种时间步策略

BATSRUS 实现了极其丰富的时间推进方案，通过 `ModMain.f90` 的逻辑变量控制：

| 方案 | 控制变量 | 适用场景 |
|------|---------|---------|
| **显式多级** | `iStage`/`nStage` | 时间精确 + 定常 |
| **局部时间步** | `UseLocalTimeStep` | 定常加速收敛 |
| **FLIC** | `UseFlic` | 混合模拟 |
| **半隐式** | `UseSemiImplicit` | 刚性源项（辐射、热传导） |
| **全隐式** | `UseFullImplicit` | 极端刚性 |
| **部分定常** | `UsePartSteady` | 部分区域已达定常 |

### 2.2 显式多级推进（核心路径）

**源文件**：`ModAdvanceExplicit.f90`，子程序 `advance_explicit`

```fortran
STAGELOOP: do iStage = 1, nStage
    ! 1. 计算面重构值 (calc_face_value)
    ! 2. 应用边界条件 (set_face_boundary)
    ! 3. 计算面通量 (calc_face_flux)
    ! 4. 应用守恒通量修正 (apply_cons_flux)
    ! 5. 计算源项 (calc_source)
    ! 6. 计算时间步 (calc_timestep)，第一阶段
    ! 7. 更新状态 (update_state)
    ! 8. 消息传递 (exchange_messages)
end do STAGELOOP
```

**关键发现**：

1. **默认是 2 级格式**：`UseHalfStep = .true.` 意味着经典的 Dt/2, Dt 两级方案（等价于二阶 Runge-Kutta / Heun 方法）

2. **nStage 可配置**：通过 `nOrder` 和 `nStage` 变量，支持更多级数

3. **时间推进格式**（非 FLIC）：
   ```
   t_{sim} = t_{sim} + Dt / nStage
   ```
   每级推进 Dt/nStage

4. **FLIC 格式**（3 级）：
   ```
   Stage 1: Dt/2，计算电磁场
   Stage 2: Dt/2，带电粒子推进（Boris 步）
   Stage 3: Dt，最终更新
   ```

### 2.3 局部时间步（Local Time Stepping）

**逻辑**：`UseLocalTimeStep = .true.`

- 每个块有独立的 `DtMax_B(iBlock)`
- 通过 `rLocalTimeStep` 半径控制分区域时间步
- 时间精确模式下用于准定常外区域加速
- `ModLocalTimeStep.f90` 处理子循环逻辑

### 2.4 部分定常（Partially Steady State）

**源文件**：`ModPartSteady.f90`

五种块类型：
```
SkippedBlock_  = 0  → 未使用
SteadyBlock_   = 1  → 已达定常，不推进
SteadyBoundBlock_=2 → 定常块边界，需推进（信息传播）
ExplBlock_     = 3  → 显式推进
ImplBlock_     = 4  → 隐式推进
```

- 仅推进 `ExplBlock_` 和 `SteadyBoundBlock_`
- 大幅减少计算量

### 2.5 半隐式方案

**源文件**：`ModSemiImplicit.f90`

处理类型（`TypeSemiImplicit`）：
- `'radiation'` — 辐射扩散
- `'radcond'` — 辐射 + 热传导
- `'cond'` — 热传导
- `'parcond'` — 平行热传导
- `'hall'` — Hall 效应

线性求解器：GMRES + BILU 预处理（Hypre/Jacobi/BlockJacobi/GS/DILU）

---

## 3. Roe 近似黎曼求解器

### 3.1 通量函数层次结构

**源文件**：`ModFaceFlux.f90` → `get_numerical_flux`（行 1156-1513）

BATSRUS 支持多种通量格式，通过 `TypeFlux` 参数选择：

| TypeFlux 值 | 方法 | 说明 |
|------------|------|------|
| `'Roe'` | 7 波 Roe 求解器 | **推荐默认**，=`UseRS7` |
| `'RoeOld'` | 经典 8 波 Roe | 旧版保留 |
| `'Rusanov'` | Lax-Friedrichs | 最简单 |
| `'Linde'` | HLL | Harten-Lax-van Leer |
| `'LFDW'` | 主导波 Lax-Friedrichs | 省计算 |
| `'HLLDW'` | 主导波 HLL | 省计算 |
| `'HLLD'` | HLLD | 五波近似 |
| `'Sokolov'` | Artificial Wind | Sokolov 方案 |
| `'Simple'` | 中心通量 | 无迎风 |

### 3.2 Roe 求解器的实际计算

**关键代码**（`ModFaceFlux.f90` 行 1529-1540）：

```fortran
subroutine roe_solver_new
    ! Roe flux = (F_L + F_R)/2 - Dissipation
    Flux_V(1:p_) = &
         0.5*(FluxLeft_V(1:p_) + FluxRight_V(1:p_)) &
         + DissipationFlux_V(1:p_)
    Flux_V(Energy_) = &
         0.5*(FluxLeft_V(Energy_) + FluxRight_V(Energy_)) &
         + DissipationFlux_V(p_+1)
    CmaxDt = Cmax
end subroutine roe_solver_new
```

**耗散通量**在 `ModCharacteristicMhd.f90` 中计算：
```
DissipationFlux = -0.5 Σ_k |λ̃_k| α_k r_k
```
其中：
- `λ̃_k` — 经 Harten 熵修正的特征值
- `α_k` — 波强度（特征变量中的跳跃）
- `r_k` — 右特征向量

### 3.3 七波 + 一发散波结构

**源文件**：`ModCharacteristicMhd.f90` 行 20-24

8 个波（命名索引）：
```
EntropyW_  = Rho_   → 熵波（对流）
AlfvenRW_  = Ux_    → Alfvén 右行波
AlfvenLW_  = Uy_    → Alfvén 左行波
SlowRW_    = Uz_    → 慢磁声右行波
FastRW_    = Bx_    → 快磁声右行波
SlowLW_    = By_    → 慢磁声左行波
FastLW_    = Bz_    → 快磁声左行波
DivBW_     = p_     → 散度波
```

注意：`UseRS7 = .true.` 时只有 7 个物理波（无散度波），但通过 `diffBn_D` 修正面法向磁场以保证跨面连续性 —— 这是 Sokolov 的 7 波方案改进。

### 3.4 特征速度计算

**源文件**：`ModCharacteristicMhd.f90` 行 386-409

```fortran
! 声速
A = sqrt(gamma * p / rho)

! Alfvén 速度
Ca = |Bn| / sqrt(rho)

! 横向磁场平方
BTang2 = |B_tang|^2

! 快磁声速（精确公式）
Cf = 0.5 * (sqrt((A-Ca)^2 + BTang2/rho) + sqrt((A+Ca)^2 + BTang2/rho))

! 慢磁声速
Cs = Ca * A / Cf
```

七波特征速度：
```
λ_1  = U_n - Cf      (快磁声左行)
λ_2  = U_n - Ca      (Alfvén 左行)
λ_3  = U_n - Cs      (慢磁声左行)
λ_4  = U_n           (熵波 / 对流)
λ_5  = U_n + Cs      (慢磁声右行)
λ_6  = U_n + Ca      (Alfvén 右行)
λ_7  = U_n + Cf      (快磁声右行)
```

### 3.5 Roe 平均

**源文件**：`ModCharacteristicMhd.f90` 行 178-211

Roe 平均公式（带密度加权）：

```
ρ_H = sqrt(ρ_L) · sqrt(ρ_R)
U_H = (sqrt(ρ_L)·U_L + sqrt(ρ_R)·U_R) / (sqrt(ρ_L) + sqrt(ρ_R))
B_H = (sqrt(ρ_L)·B_R + sqrt(ρ_R)·B_L) / (sqrt(ρ_L) + sqrt(ρ_R)) + 0.5(BnL+BnR)·n̂
```

声速平均包含额外项（用于处理激波）：
```
a²_H = (sqrt(ρ_L)·a²_L + sqrt(ρ_R)·a²_R)/(sqrt(ρ_L)+sqrt(ρ_R)) 
     + γ·X_H + (γ-1)·[Xn_H + 0.5·Σ(dU)²·ρ_H·(sqrt(ρ_L)+sqrt(ρ_R))⁻²]
```
其中 `X_H` 和 `Xn_H` 分别是法向和切向磁场跳跃的贡献项。

### 3.6 Harten 熵修正

**源文件**：`ModCharacteristicMhd.f90` 行 415-447

```fortran
subroutine get_fixed_abs_eigenvalue(Eigenvalue_V, EigenvalueL_V, EigenvalueR_V, &
                                    LambdaB0, EigenvalueFixed_V, CMax, IsBoundary)
    Eps_V = max(LambdaB0, abs(λ_fastR - λ_fastL) * 0.05)
    do iWave = 1, nVar-1
        λ = Eigenvalue_V(iWave)
        Eps_V(iWave) = max(λ_R - λ, λ - λ_L, Eps_V(iWave))
        EigenvalueFixed_V(iWave) = max(|λ|, Eps_V(iWave))
    end do
    CMax = max(EigenvalueFixed_V(FastRW_), EigenvalueFixed_V(FastLW_))
end subroutine
```

- `LambdaB0`：基于 B0 变化的附加耗散（处理 B0 尖锐梯度）
- `Eps_V`：确保 `|λ̃| ≥ max(|λ|, ε)`，其中 `ε = max(0.05*|λ_fastR-λ_fastL|, LambdaB0)`
- 可通过 `Climit > 0` 限制最大波速以减少数值耗散

### 3.7 七波方案的 Bn 连续性修正

**源文件**：`ModFaceFlux.f90` 行 1326-1344

```fortran
if(UseRS7 .or. UseLindeFix)then
    ! 计算面法向磁场跳跃
    DiffBn_D = n̂ * 0.5*Σ[(B_R - B_L)·n̂]
    
    ! 消除跳跃（强制跨面连续性）
    B_L = B_L + DiffBn_D
    B_R = B_R - DiffBn_D
    
    ! 修正能量跳跃
    DiffE = 0.5*Σ[(B_R + B_L)·DiffBn_D]
    DiffBb = Σ(DiffBn_D²)
end if
```

并通过 `modify_flux` 相应调整动量通量：
```fortran
Flux(RhoUx_:RhoUz_) += 0.5 * DiffBb * n̂
Flux(Energy_)       += Un * DiffBb
```

---

## 4. 高阶重构（MUSCL/TVD/WENO）

### 4.1 重构子程序

**源文件**：`ModFaceValue.f90`

主要重构函数：

| 函数 | 类型 | 精度 |
|------|------|------|
| `get_face_tvd` | TVD/MUSCL | 二阶 |
| `get_face_accurate1d` | 一维高精度 | 高阶 |
| `get_face_accurate2d` | 二维高精度 | 高阶 |
| `get_face_accurate3d` | 三维高精度 | 高阶 |

### 4.2 默认配置

```fortran
TypeLimiter = 'minmod'        ! 默认限制器
LimiterBeta = 1.0             ! 限制器压缩参数
UseAccurateResChange = .false. ! 分辨率变化使用精确重构（默认关闭）
UseTvdResChange = .true.      ! 分辨率变化使用 TVD（默认开启）
nGUsed = nG                   ! 使用的 ghost cell 数量
```

### 4.3 可用限制器

BATSRUS 支持多种限制器，默认 `minmod`：
- **minmod** — 最耗散、最鲁棒
- **mc** (monotonized central) — 中等耗散
- **vanleer** — 平滑限制
- **superbee** — 最小耗散（易振荡）
- **mp** (monotonized central-5) — 五阶精度限制器

### 4.4 CWENO5（中心 WENO 五阶）

```fortran
UseCweno = .false.   ! 默认关闭
TypeLimiter5 = 'mp'  ! 五阶备选限制器
```

### 4.5 对数限制器

MHD 密度和压强经常跨越数量级，BATSRUS 实现了对数重构：

```fortran
UseLogLimiter      ! 对密度和压强用对数限制
UseLogRhoLimiter   ! 密度对数限制
UseLogPLimiter     ! 压强的对数限制
```

对极大密度梯度区域（如行星际空间），对数重构大幅提升鲁棒性。

### 4.6 低阶回退

```fortran
UseLowOrder = .false.      ! 部分面使用低阶
UseLowOrderRegion = .false. ! 区域低阶
UseAdaptiveLowOrder = .false. ! 自适应低阶
IsLowOrderOnly_B(:)         ! 每块的低阶标记
```

在激波附近或高梯度区域可以自动回退到低阶通量以保证鲁棒性。

---

## 5. 时间步约束（CFL 条件）

### 5.1 单单元时间步计算

**源文件**：`ModTimeStepControl.f90`，子程序 `calc_timestep`（行 122-200）

```fortran
! 每个单元的时间步
Vdt = max(Flux_VXI(Vdt_,i,j,k), Flux_VXI(Vdt_,i+1,j,k))   ! X 方向
if(nJ > 1) Vdt += max(Flux_VYI(Vdt_,i,j,k), Flux_VYI(Vdt_,i,j+1,k))  ! Y 方向
if(nK > 1) Vdt += max(Flux_VZI(Vdt_,i,j,k), Flux_VZI(Vdt_,i,j,k+1))  ! Z 方向
DtMax_CB(i,j,k,iBlock) = CellVolume / Vdt
```

其中 `Vdt_` 是 CmaxDt × Area（在 `get_numerical_flux` 末尾设置）：

```fortran
! 总波速 = 面波速 + 耗散贡献（电阻率、粘性、热传导、辐射扩散）
CmaxDt = CmaxDt + 2*(DiffCoef + ViscoCoeff + Eta)*InvDxyz
```

### 5.2 全局时间步

```fortran
! DtMax_B(iBlock) = Cfl * min(DtMax_CB(:,:,:,iBlock))
! Dt = min(DtMax_B(1:nBlock))

if(IsTimeAccurate)then
    Dt = DtFixed              ! 时间精确：固定 Dt
else
    Dt = Cfl * minval(DtMax_CB)  ! 定常：局部 CFL
end if
```

### 5.3 CFL 强制

```fortran
DoEnforceCfl = .false.        ! 默认不强制
```

在半隐式热传导导致局部温度骤升时，可开启强制 CFL 修正时间步。

### 5.4 时间步控制

**源文件**：`ModTimeStepControl.f90` 行 42-51

自适应时间步调整：
```fortran
RejectStepLevel1   = 0.3    ! 变化 > 30% → 拒绝
RejectStepLevel2   = 3.0    ! 变化 > 300% → 拒绝（局部时间步）
RejectStepFactor   = 0.50   ! 拒绝后时间步乘以 0.5
ReduceStepLevel1   = 0.6    ! 变化 > 60% → 减小
ReduceStepLevel2   = 1.5    ! 变化 > 150% → 减小
ReduceStepFactor   = 0.95   ! 减小因子
IncreaseStepLevel1 = 0.8    ! 变化 < 80% → 可增大
IncreaseStepLevel2 = 1.2    ! 变化 < 120% → 可增大
IncreaseStepFactor = 1.05   ! 增大因子
```

---

## 6. 守恒性与精度

### 6.1 空间精度

| 组件 | 方法 | 精度 |
|------|------|------|
| 重构 | MUSCL + minmod/mc/van Leer 限制器 | 二阶 |
| 通量 | Roe 近似黎曼求解器（7 波方案）| 迎风型 |
| 源项 | 单元中心有限体积 | 二阶 |
| CWENO5 | 中心 WENO 五阶 | 五阶（可选） |

### 6.2 时间精度

| 模式 | 方法 | 精度 |
|------|------|------|
| 时间精确 | Dt/2, Dt 两级 Runge-Kutta | 二阶 |
| 定常 | 局部时间步 + 多级 | 一阶（仅稳态） |
| BDF2 | 二阶后向差分（隐式） | 二阶 |

### 6.3 守恒性

1. **通量守恒**：`DoConserveFlux = .true.`（默认启用在分辨率变化处）
   - 粗-细界面通量相等
   - 通过 `save_cons_flux`/`apply_cons_flux` 在 `ModConserveFlux.f90` 中实现

2. **总变差递减（TVD）**：迎风格式 + minmod 限制器保证标量方程 TVD

3. **正性保持**：`check_positivity` 子程序在每次更新后检查密度和压强非负

### 6.4 精度限制因素

1. **MHD 散度控制**：8 波方案不严格满足 ∇·B=0，引入了额外耗散
2. **B0 场处理**：内禀磁场 B0 的梯度被熵修正感知，增加了额外数值耗散
3. **Climit**：限速可降低耗散但也可能造成激波附近的非物理振荡

---

## 7. 八波方案的 divB 处理

### 7.1 散度控制选项

**源文件**：`ModMain.f90` 行 206-214

```fortran
UseDivbSource     = UseB    ! 8 波源项（Powell 方案，默认开）
UseDivbDiffusion  = .false.  ! 抛物型散度扩散
UseProjection     = .false.  ! 投影法（椭圆修正）
UseConstrainB     = .false.  ! 约束传输（CT，Yee 网格）
UseHyperbolicDivb = .false.  ! 双曲/抛物清洗（Dedner GLM 方案）
SpeedHypDim = -1.0            ! 双曲清洗波速
HypDecay = 0.1                ! 散度衰减率
```

### 7.2 八波方案原理

Powell 八波方案的核心：

1. **修正 MHD 方程**：在动量、能量和感应方程中添加与 ∇·B 成正比的源项
   ```
   ∂(ρu)/∂t + ∇·[ρuu + (p + B²/2)I - BB] = -B(∇·B)
   ∂e/∂t + ∇·[(e + p + B²/2)u - B(u·B)] = -(u·B)(∇·B)
   ∂B/∂t + ∇·(uB - Bu) = -u(∇·B)
   ```

2. **特征系统中引入第 8 个波**：散度波以流速 U_n 传播

3. **不修改黎曼求解器**：源项用分裂法单独处理

### 7.3 七波改进（当前 Roe 默认）

当前默认方案（`DoRoe = 'Roe'`，即 `UseRS7 = .true.`）：

1. 通过 `diffBn_D` 强制面法向磁场连续性
2. 通过 `modify_flux` 修正动量/能量通量
3. 等价于消除了第 8 个散度波，保留 7 个物理波

这种方式结合了 Powell 源项和 Sokolov 的 7 波方案，在保证 Galilean 不变性的同时提高了精度。

### 7.4 Linde 修正

`UseLindeFix` 提供额外保护：
```fortran
if(.not.UseRS7)then
    Flux_V(Bx_:Bz_)  = Flux_V(Bx_:Bz_) - Cmax*DiffBn_D  ! Lax-Friedrichs for Bn
    Flux_V(Energy_)  = Flux_V(Energy_) - Cmax*DiffE      ! 能量修正
end if
```

---

## 8. 特殊数值技术

### 8.1 Boris 修正

**源文件**：`ModBorisCorrection.f90`，`ModPhysicalFlux.f90` → `get_boris_flux`

在近相对论区域，动量通过 (1 + V_A²/c²) 因子修正：
```
ρu_Boris = ρu + (u B² - B(u·B))/c²
E_Boris  = E + E²/(2c²)
```

通量包含电场的 Poynting 通量贡献。

### 8.2 B0 内禀场处理

势场/无力场 B0 被分为内禀部分处理：
- 面重建时施加 B0 梯度
- 耗散通量中通过 `LambdaB0` 添加额外耗散以稳定 B0 尖锐梯度

### 8.3 网格自适应（AMR）

**源文件**：`ModAMR.f90`

- 块结构：每个块 nI×nJ×nK 个单元（自相似）
- 细化判据：几何（`geo`）、物理（`phy`）、混合（`all`）
- 初细化级别：`nInitialAmrLevel`
- 动态细化：`DoAutoRefine`
- 分辨率变化处的通量守恒
- 延长（prolongation）阶：`nOrderProlong = 1`

### 8.4 扩展物理模块

BATSRUS 支持以下扩展物理：
- Hall MHD（`UseHallResist`）
- 双流体/多离子 MHD（`UseMultiIon`）
- 各向异性压强（`UseAnisoPressure`）
- 电子压强方程（`UseElectronPressure`）
- 辐射扩散（`UseRadDiffusion`）
- 热传导（`UseHeatConduction`）
- 电阻 MHD（`UseResistivity` / `Eta_GB`）
- 阿尔文波压强（`UseWavePressure`）
- 多物种追踪（`UseMultiSpecies`）

---

## 总结

BATSRUS 求解器是一个功能极其丰富的 MHD 求解框架，其核心数值方法可概括为：

| 特性 | 实现 |
|------|------|
| **空间离散** | 有限体积 + 块自适应网格（BAT） |
| **重构** | MUSCL/TVD 二阶 / CWENO5 五阶 |
| **黎曼求解器** | 7 波 Roe 迎风（默认）/ HLLD / HLL / Lax-Friedrichs / Sokolov |
| **时间推进** | Dt/2-Dt 二阶显式 RK / BDF2 隐式 / FLIC 三级 |
| **定常加速** | 局部时间步 / 部分定常 / 半隐式 |
| **divB 控制** | 8 波源项 + 7 波方案 / 双曲清洗 / CT / 投影 |
| **并行** | MPI + OpenMP + OpenACC (GPU) |
| **精度** | 空间二阶（可达五阶），时间二阶，严格守恒 |

**核心设计哲学**：以模块化方式组合不同精度的方案，在不同物理区域自动选择最适方法 —— 从显式到隐式、从低阶到高阶、从全局时间步到局部时间步，一切皆可灵活配置。

---

> **引用**
>
> 1. Powell, K.G., Roe, P.L., Linde, T.J., Gombosi, T.I., De Zeeuw, D.L., "A Solution-Adaptive Upwind Scheme for Ideal Magnetohydrodynamics", J. Comp. Phys., 154, 284-309, 1999.
> 2. Tóth, G., et al., "A Parallel Explicit/Implicit Time Stepping Scheme on Block-Adaptive Grids", J. Comp. Phys., 2006.
> 3. BATSRUS GPU: "Faster-than-Real-Time Magnetospheric Simulations with a Block-Adaptive Grid Code", ApJ, 2025.
> 4. GitHub 仓库：https://github.com/SWMFsoftware/BATSRUS
