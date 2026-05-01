# BATSRUS 自适应网格系统 (AMR) 深度分析报告

> **分析者**: batsrus_grid (he我)
> **生成时间**: 2026-05-01
> **源代码仓库**: https://github.com/SWMFsoftware/BATSRUS
> **分析范围**: 21 个 Fortran 源文件，涵盖 AMR、数据结构、求解器接口、守恒性、负载均衡全局架构

---

## 一、系统概览

BATSRUS (**B**lock-**A**daptive-**T**ree **S**olar-wind **R**oe-type **U**pwind **S**cheme) 是密歇根大学 CSEM 中心开发的大规模三维磁流体力学 (MHD) 模拟代码，是 Space Weather Modeling Framework (SWMF) 的核心求解器。

### 关键数字

| 参数 | 说明 |
|---|---|
| 块尺寸 | 典型 $8\times8\times8$ 至 $16\times16\times16$（由 `BATL_size` 中 `nI, nJ, nK` 定义） |
| 块上限 | `MaxBlock`——每 MPI 进程可持有的最大块数 |
| 幽灵层 | `nG` 层（通常 $nG=2$，支持二阶精度） |
| 细化比 | $2:1$（每个方向分裂为 2） |
| 3D 细化 | $1\to 8$ 块（八叉树/octree） |

---

## 二、块自适应树数据结构

### 2.1 核心抽象：Block vs Node

BATSRUS 有两层数据结构：

**Block（块）**——承载物理数据的矩形盒子：
- `State_VGB(nVar, MinI:MaxI, MinJ:MaxJ, MinK:MaxK, MaxBlock)`——三维守恒变量数组
- 含 `nG` 层幽灵单元用于 MPI 通信和边界条件
- `ModBlockData.f90` 提供块级动态存储（`put_block_data` / `get_block_data`）

**Node（节点）**——块树的逻辑节点：
- `iTree_IA(Status_:Block_, MaxNode)`——定义在底层 `BATL_lib` 中
- 关键索引：`Status_`（Used/Unused）、`Proc_`（所属进程）、`Block_`（本地块 ID）
- `ModNodes.f90` 提供节点级统一编号（`iNodeLocal_NB`、`iNodeGlobal_NB`）

**块树结构总结**：

```
Node Tree (iTree_IA)           Block Array (State_VGB)
┌─────────────────┐            ┌──────────────────────┐
│ Node 0 (root)   │──映射到──→│ Block on Proc 0       │
│  ├─ Node 1      │──映射到──→│ Block on Proc 0       │
│  ├─ Node 2      │──映射到──→│ Block on Proc 1       │
│  │  ├─ Node 3   │──映射到──→│ Block on Proc 1       │
│  │  └─ ...      │            │ ...                   │
│  └─ ...         │            └──────────────────────┘
└─────────────────┘
```

### 2.2 八叉树细化

虽然 BATSRUS 源码中未见显式的八叉树代码（该逻辑封装在 `BATL_lib` 的 `regrid_batl` 中），但从 `ModAMR.f90` 的调用链可以推断：

- **细化比 = 2:1**：每个方向均分裂为 2
- **3D 场景**：$1 \to 8$ 个子块
- **约束**：相邻块之间的细化级别差 $\le 1$（通过 `DiLevel_EB` 跟踪）
- **`DiLevel_EB(6, MaxBlock)`**：6 个面各记录与邻居的层级差（$-1, 0, +1$）

### 2.3 空间填充曲线与 Morton 序

`ModLoadBalance.f90` 第 275 行的注释明确说明：
```
! Load balance grid using space filling (Morton) ordering of blocks
```

`regrid_batl` 调用使用 Morton 序（Z-order）将多维块树线性化，实现以下目标：
- 空间局部性保持
- 负载均衡重分布
- 按层级分别负载均衡（`DoBalanceEachLevelIn`）

---

## 三、AMR 控制流程

### 3.1 调用链

```
ModMain (主循环)
  └─ prepare_amr()          ← ModAMR.f90:171
       ├─ exchange_messages()   ← 二阶 ghost cell 交换
       ├─ amr_criteria()        ← 计算细化判据
       └─ set_amr_criteria()    ← 将判据送入 BATL_lib
  └─ do_amr()               ← ModAMR.f90:236
       ├─ regrid_batl()         ← 核心：执行细化/粗化/重分布
       ├─ set_batsrus_grid()    ← 更新几何信息
       ├─ count_true_cells()    ← 统计有效单元数
       ├─ clean_block_data()    ← 清理动态数据
       ├─ set_batsrus_state()   ← 修复新块/移动块的物理状态
       ├─ load_balance()        ← Morton 序负载均衡
       └─ exchange_messages()   ← 最终消息传递
```

### 3.2 细化/粗化频率控制

由 `AdaptGrid`（类型 `FreqType`）控制：
- `#AMR` 命令：`DnRefine`——每 N 步执行一次 AMR
- `#DOAMR` 命令：`DoAmr`、`DnAmr`、`DtAmr`——按步数或物理时间触发

### 3.3 层级控制

- `#AMRLEVELS`：`MinAmrLevel`（最粗层级）、`MaxAmrLevel`（最细层级）
- `#AMRRESOLUTION`：`CellSizeMin`、`CellSizeMax`——按物理单元尺寸控制
- 由 `iTree_IA(MinLevel_,:)` 和 `iTree_IA(MaxLevel_,:)` 实施

---

## 四、细化/粗化判据（AMR Criteria）

### 4.1 架构

判据通过 `#AMRCRITERIA`、`#AMRCRITERIALEVEL`、`#AMRCRITERIARESOLUTION` 等命令配置。核心子程序 `amr_criteria()`（`ModAMR.f90:391`）计算每个块的标量判据值，然后 `BATL_lib::set_amr_criteria()` 根据这些值决定细化/粗化。

### 4.2 14 种内置物理判据

| 判据名 | 物理量 | 说明 |
|---|---|---|
| `gradt` | $\|\nabla T\|$ | 温度梯度——追踪激波和边界层 |
| `gradlogrho` | $\|\nabla\log_{10}\rho\|$ | 对数密度梯度 |
| `gradp` | $\|\nabla p\|$ | 压力梯度 |
| `gradlogp` | $\|\nabla\log_{10}p\|$ | 对数压力梯度 |
| `-divu` | $-\nabla\cdot\mathbf{u}$ | 速度散度（负号捕获激波压缩区） |
| `-divudx` | $-\Delta x\|\nabla\cdot\mathbf{u}\|$ | 乘以单元尺寸的速度散度 |
| `pjumpratio` | $\max(p_{\text{neighbor}}/p_{\text{center}})$ | 相邻单元最大压力跳变比 |
| `j` | $\|\mathbf{J}\|$ | 电流密度模——追踪电流片 |
| `j2` | $\|\mathbf{J}\|^2$ | 电流密度平方——更敏感的电流片指示 |
| `currentsheet` | $\text{sign}(\mathbf{r}\cdot\mathbf{B})$ | 电流片检测——通过径向磁场反转 |
| `user` | 用户自定义 | 可扩展接口 |

### 4.3 瞬态/动态判据（trace_transient）

通过比较当前状态 `State_VGB` 与旧状态 `StateOld_VGB`：

| 判据名 | 公式 |
|---|---|
| `rho_dot` | $\|\rho-\rho^{\text{old}}\|/\max(\rho,\rho^{\text{old}}, \epsilon)$ |
| `p_dot` | $\|p-p^{\text{old}}\|/\max(p,p^{\text{old}}, \epsilon)$ |
| `T_dot` | $\|T-T^{\text{old}}\|/\max(T,T^{\text{old}}, \epsilon)$ |
| `rhoU_dot` | $\|\|\rho\mathbf{u}\|-\|\rho\mathbf{u}\|^{\text{old}}\|/\max(\dots)$ |
| `B_dot` | $\|\|\mathbf{B}\|-\|\mathbf{B}\|^{\text{old}}\|/\max(\dots)$ |
| `meanUB` | $(d\|\rho\mathbf{u}\|/\|\rho\mathbf{u}\|) \times (d\|\mathbf{B}\|/\|\mathbf{B}\|)$ |
| `Rho_2nd_1` | $(\|d^2\rho/dx^2\|+\|d^2\rho/dy^2\|+\|d^2\rho/dz^2\|)/\rho$ |
| `Rho_2nd_2` | $\|d^2\rho/dx^2+d^2\rho/dy^2+d^2\rho/dz^2\|/\rho$ |

### 4.4 判据工作流

```
amr_criteria() 对每个活跃块:
  1. 提取密度、动量、磁场、压力到局部数组
  2. 对每个判据:
     - calc_gradient() 计算梯度
     - 取块内 maxval / sqrt(maxval(grad²))
     - 量纲化（乘以 No2Io_V(Unit*)）
  3. Crit_IB(iCrit, iBlock) = 标量值
  4. set_amr_criteria() 根据判据+阈值决定细化/粗化
```

---

## 五、守恒性保证

### 5.1 通量修正法（Flux Correction）

BATSRUS 采用 Berger-Colella 式通量修正法保证粗细界面守恒。

**ModConserveFlux.f90**：

- `save_cons_flux(iBlock)`：在粗网格面保存通量（仅当 `DiLevel_EB == 1`，即与细网格相邻时）
  - 各方向独立处理：`save_corrected_flux_x/y/z`
  - 保存 `Flux_VXI/VYI/VZI`（$Vdt$）和磁场法向分量（$B_{nL}, B_{nR}$）
  - 对面面积加权（笛卡尔使用 `FaceRatio = 1/2^(nDim-1)`，曲线坐标使用 `FaceNormal_DDFB`）

- `apply_cons_flux(iBlock)`：将保存的粗通量替换细网格面上的通量

**原理**：
```
    粗网格面 F_c
    ┌────────┬────────┐
    │        │        │
    │ 粗块   │ 粗块   │
    │        │        │
    ├──┬──┬──┼──┬──┬──┤  ← 细网格面 f_1, f_2
    │  │  │  │  │  │  │
    │细│细│细│细│细│细│
    │  │  │  │  │  │  │
    └──┴──┴──┴──┴──┴──┘

    守恒条件: F_c = f_1 + f_2
    修正: f_1_corrected + f_2_corrected = F_c
```

### 5.2 守恒/非守恒混合方案

**ModConservative.f90**：

- `UseNonConservative`——全局开关
- `#CONSERVATIVECRITERIA`——基于物理判据选择单元更新方式
  - `r`（半径）：$r < r_{\text{Conserv}}$ 用非守恒
  - `parabola`：抛物面内用非守恒
  - `p`（Balsara/Ryu 开关 1）：压力跳变检测激波
  - `gradp`/`jumpp`（Balsara/Ryu 开关 2）：压力梯度
- 混合方案优缺点：守恒格式在激波处更准确，非守恒格式在平滑流中更稳定

---

## 六、负载均衡

### 6.1 架构

`ModLoadBalance.f90::load_balance()` 采用三级策略：

**第一层：块类型分类（位标记法）**

每个块根据多个属性被赋予一个整数类型码（bitmask）：

| 位 | 含义 |
|---|---|
| bit 0 | `NotSkippedBlock`——是否活跃 |
| bit 1 | `TrueBlock`——是否为实单元（非 body） |
| bit 2 | `PartImplicit`——部分隐式 |
| bit 3 | `SemiImplicit`——半隐式 |
| bit 4 | `PointImplicit`——点隐式 |
| bit 5 | `FieldLineThread`——场线追踪边界 |
| bit 6 | `PicBlock`——PIC 粒子 |
| bit 7 | `SteadyBlock`——稳态块 |
| bit 8 | `HighOrderBlock`——高阶格式 |
| bits 9+ | `UserType`——用户自定义 |
| 高位 | `SubCycleBlock`——子循环时间层级 |

**第二层：MPI 全局收集**

```fortran
MPI_allgatherv(iTypeAdvance_B → iTypeAdvance_BP)  ! 全局块类型视图
MPI_allreduce(iTypeBalance_A)                       ! 节点级类型汇总
```

**第三层：Morton 序重分布**

```fortran
regrid_batl(nVar, State_VGB, DtMax_B,        &
    iTypeBalance_A=iTypeBalance_A,            &  ! 按类型组重分布
    iTypeNode_A=iTypeAdvance_A,               &  ! 节点类型信息
    DoBalanceEachLevelIn=UseLocalTimeStep,    &  ! 按层级分别均衡
    DoBalanceOnlyIn=.true.,                   &  ! 仅均衡不细化
    nExtraData=nBuffer,                       &  ! 额外数据缓冲区
    pack_extra_data=pack_load_balance,        &  ! 序列化回调
    unpack_extra_data=unpack_load_balance)       ! 反序列化回调
```

### 6.2 数据打包/解包

`pack_load_balance` 和 `unpack_load_balance` 处理块迁移时的数据传递：
1. 标量标志（动态数据量、点隐式标志）
2. 面约束磁场 `BxFace_GB/ByFace_GB/BzFace_GB`
3. BDF2 隐式旧态 `ImplOld_VCB`
4. 场线追踪状态 `Trace_DSNB`
5. 用户自定义块数据

---

## 七、梯度计算与面值重构

### 7.1 梯度计算 (ModCellGradient.f90)

三种梯度计算接口：

- `calc_gradient1(iBlock, Var_G, Grad_DG)`——返回向量梯度（$nDim$ 分量）
- `calc_gradient3(iBlock, Var_G, GradX_C, GradY_C, GradZ_C)`——返回三个标量分量
- `calc_gradient_ghost`——梯度扩展到 ghost cells

**算法**：二阶中心差分（笛卡尔网格）
$$
\frac{\partial f}{\partial x}\bigg|_{i} = \frac{f_{i+1} - f_{i-1}}{2\Delta x}
$$

对于 body 单元（行星内部），使用单侧差分外推。

### 7.2 散度与旋度

- `calc_divergence(iBlock, Var_DG, nG, Div_G)`——在笛卡尔和广义坐标均可用
- `calc_cell_curl_ghost`——含一层 ghost 的旋度

### 7.3 LIMITER（斜率限制器）

BATSRUS 的 MUSCL-Hancock 格式中的斜率限制在面值计算中实现（在 `ModFaceFlux.f90` 和相关面值模块中，不属于本次分析范围）。从判据中对 `calc_gradient` 的调用可知系统使用二阶精度重建，加上适当的限制器防止振荡。

---

## 八、坐标系统支持

| 系统 | 标志 | 说明 |
|---|---|---|
| 笛卡尔 (Cartesian) | `IsCartesian` | 最简单，面积恒定 |
| 笛卡尔网格 | `IsCartesianGrid` | 正交笛卡尔 |
| 柱坐标 (RZ) | `IsRzGeometry` | 轴对称，$r$-$z$ 平面 |
| 球坐标 (r-lon-lat) | `IsRLonLat` | 全球磁层模拟标准 |
| 对数半径 | `IsLogRadius` | $\log(r)$ 拉伸 |
| 广义半径 | `IsGenRadius` | 用户定义 $r\to\xi$ 映射 |
| 柱对称轴 | `IsCylindricalAxis` | 柱坐标极轴 |
| 球对称轴 | `IsSphericalAxis` | 球坐标极轴 |

---

## 九、极轴特殊处理

### 9.1 轴粗化 (ModCoarseAxis.f90)

在球坐标极轴附近，经度方向单元退化。`coarsen_axis_cells()` 通过合并相邻 $j$ 方向的单元来缓解 Courant 条件限制：

```
nCoarseLayer = 3:
┌───────────────────────┐  第1层: 合并 8 个单元
│         axis          │
├───────────┬───────────┤  第2层: 合并 4 个单元
│           │           │
├─────┬─────┼─────┬─────┤  第3层: 合并 2 个单元
│     │     │     │     │
├──┬──┼──┬──┼──┬──┼──┬──┤  第4层: 原始单元
│  │  │  │  │  │  │  │  │
└──┴──┴──┴──┴──┴──┴──┴──┘
```

### 9.2 轴修正 (ModFixAxisCells.f90)

`fix_axis_cells()` 对极轴附近单元进行线性拟合修正：
- 收集全局各进程在轴的物理量总和
- 利用 `MPI_allreduce` 获得全局平均值
- 用受限制的梯度（`Beta` 参数控制斜率限制）重构轴单元值
- 支持一层轴修正 (`rFixAxis`) 和两层轴修正 (`r2FixAxis`)

---

## 十、消息传递与 Ghost Cells

### 10.1 核心接口

`ModMessagePass.f90::exchange_messages()` 提供：
- `UseOrder2In`——是否使用二阶 ghost cell 交换
- `DoResChangeOnlyIn`——仅更新分辨率变化面

### 10.2 AMR 中的消息传递时序

1. **AMR 前**（`prepare_amr`）：二阶 ghost cell 交换——确保判据计算有正确的邻居值
2. **AMR 后**（`do_amr`）：
   - 若网格未变：仅分辨率变化面 ghost cell 更新
   - 若网格已变：全量消息传递 + GPU 同步 + 负载均衡

### 10.3 Ghost Cell 用途

- 二阶梯度/散度计算需要 $nG=2$ 层 ghost
- 通量计算需要 $nG=1$ 层 ghost（Riemann 求解器）
- 粗细界面：`DiLevel_EB` 指示 ghost 填充是否需要插值

---

## 十一、GPU 支持

BATSRUS 通过 OpenACC 指令支持 GPU 加速（arXiv:2501.06717）：

| 模块 | GPU 特性 |
|---|---|
| `ModAMR.f90` | `sync_cpu_gpu('update on GPU', ...)`——CPU/GPU 状态同步 |
| `ModCoarseAxis.f90` | `!$acc parallel`、`!$acc loop`——极轴粗化并行化 |
| `ModUpdateStateFast.f90` | `sync_cpu_gpu`——快速状态更新接口 |
| 性能 | 1× A100 ≈ 270× AMD Rome CPU 核心 |

---

## 十二、架构亮点与局限

### 亮点

1. **块树 + Morton 序**——经典高效的并行 AMR 范式
2. **丰富物理判据**——14+ 种判据覆盖激波、电流片、不稳定性
3. **多坐标系统**——笛卡尔/柱/球坐标 + 对数/广义拉伸
4. **守恒通量修正**——Berger-Colella 式保证守恒
5. **GPU 支持**——通过 OpenACC 迁移到 GPU
6. **位标记负载均衡**——灵活的多维块分类

### 局限/可改进之处

1. **块树逻辑封闭**——`BATL_lib` 是预编译库，无法从 Fortran 源码级理解八叉树实现
2. **细化比固定为 2**——不支持各向异性细化
3. **块尺寸固定**——全网格使用相同 `nI×nJ×nK`
4. **负载均衡粒度**——以块为单位，无法在单元级均衡
5. **MPI 通信模式**——基于 CPU 的 MPI（即使 GPU 计算）
6. **无显式 LIMITER**——面值重构中的限制器在独立模块中

---

## 参考文献

1. Tóth et al. (2012), "Adaptive numerical algorithms in space weather modeling", *J. Comput. Phys.*
2. Gombosi et al. (2025), "BATSRUS GPU: Faster-than-Real-Time Magnetospheric Simulations", *ApJ*, arXiv:2501.06717
3. CCMC Model Description: https://ccmc.gsfc.nasa.gov/models/SWMF~GM=BATSRUS~20180525/
4. BATSRUS GitHub: https://github.com/SWMFsoftware/BATSRUS

---

*本报告基于对 BATSRUS 开源代码 (`src/` 目录下 21 个 Fortran 源文件) 的人工研读与深度分析。AMR 底层逻辑 (BATL_lib) 为预编译库，其八叉树实现细节未在 Fortran 源中暴露。*
