# BATSRUS MPI 并行架构与负载均衡深度分析

> 分析日期：2026-05-01  
> 分析者：batsrus_parallel（deepseek_pro 化身）  
> 源代码：https://github.com/SWMFsoftware/BATSRUS  
> BATL 库：https://github.com/SWMFsoftware/BATL

---

## 目录

1. [总体架构概览](#1-总体架构概览)
2. [MPI 域分解策略](#2-mpi-域分解策略)
3. [通信模式](#3-通信模式)
4. [负载均衡算法](#4-负载均衡算法)
5. [GPU 加速：OpenACC](#5-gpu-加速openacc)
6. [共享内存并行：OpenMP](#6-共享内存并行openmp)
7. [可扩展性](#7-可扩展性)
8. [代码架构总结](#8-代码架构总结)

---

## 1. 总体架构概览

BATSRUS（Block-Adaptive-Tree Solar-wind Roe-type Upwind Scheme）是 Space Weather Modeling Framework (SWMF) 中使用最广泛、计算资源消耗最大的模块。其核心并行架构基于三层设计：

```
┌─────────────────────────────────────────────────┐
│              BATSRUS 并行架构                    │
├─────────────────────────────────────────────────┤
│  MPI (分布式内存)   ←→  跨节点/处理器粗粒度并行    │
│        ↕                                           │
│  OpenACC (GPU)      ←→  细粒度数据并行            │
│        ↕                                           │
│  OpenMP (共享内存)  ←→  节点内多线程加速            │
└─────────────────────────────────────────────────┘
```

**核心库：BATL (Block-Adaptive Tree Library)**

2012 年，Tóth、van der Holst、Sokolov 等人将 BATS-R-US 的核心重构为独立的 BATL 库（发表于 J. Comp. Phys. 231, 870–903, 2012），提供：

- **AMR**（自适应网格细化）：块状八叉树结构
- **消息传递**：块间 ghost cell 通信
- **负载均衡**：基于空间填充曲线的块分配

BATSRUS 自身约 20 万行现代 Fortran 代码，通过调用 `BATL_lib` 模块的公共接口来使用这些并行能力。

---

## 2. MPI 域分解策略

### 2.1 基于块的空间分解

BATSRUS 采用**基于块的域分解**（Block-based Domain Decomposition），而非传统的基于网格单元的分解。整个计算域由若干**矩形块**（rectangular blocks）组成，这些块被组织在一棵自适应的八叉树中。

**关键数据结构**（定义于 `ModParallel.f90` 和 BATL 树）：

```
                    ┌──────────────┐
                    │  Root Block  │  (level 0)
                    └──────┬───────┘
              ┌────────────┼────────────┐
         ┌────┴────┐  ┌────┴────┐  ┌────┴────┐
         │ Block A │  │ Block B │  │ Block C │   (level 1)
         └────┬────┘  └─────────┘  └─────────┘
       ┌──────┼──────┐
   ┌───┴───┐ ┌───┴───┐ ┌───┴───┐
   │  A.1  │ │  A.2  │ │  A.3  │              (level 2)
   └───────┘ └───────┘ └───────┘
```

### 2.2 邻居关系编码

每个块的邻居信息存储在两个关键数组中（`ModParallel.f90`）：

| 变量 | 维度 | 含义 |
|------|------|------|
| `DiLevel_EB(1:6, MaxBlock)` | 6方向 × 最大块数 | 邻居块的AMR层级差（0=同级, -1=低一级, +1=高一级, Unset_=无邻居） |
| `jProc_IEB(4, 1:6, MaxBlock)` | 4子面 × 6方向 × 最大块数 | 6个方向中每个方向最多4个子块的进程号 |
| `jBlock_IEB(4, 1:6, MaxBlock)` | 4子面 × 6方向 × 最大块数 | 对应邻居的块索引 |

**关键设计**：由于相邻块之间最多允许一级的 AMR 层级差，因此每个方向最多有 1 个（同级或高级邻居）或 4 个（低一级邻居的子块）相邻块。

方向编码（1-6）：-X, +X, -Y, +Y, -Z, +Z

### 2.3 进程-块映射

BATL 维护一个全局八叉树 `iTree_IA(:,:)`，记录每个树节点的状态：

| 字段 | 含义 |
|------|------|
| `Status_` | 节点状态（Used_、Unused_ 等） |
| `Proc_` | 处理此节点的 MPI 进程号 |
| `Block_` | 节点对应的块索引 |
| `iTimeLevel_A` | 时间层级（用于子循环） |

每个进程只存储分配给它的块的实际数据。通过 `iTree_IA` 可以查询任意树节点所在的进程和块号。

### 2.4 全局通信优化

BATSRUS 使用 `MPI_allgatherv` 而非 `MPI_allgather` 来优化全局通信（`ModParallel.f90`）：

```fortran
! MPI_allgatherv 的优化版本
! MaxBlockDisp_P(jProc) = jProc * MaxBlock  — 固定位移
! nBlockMax_P(iProc)    = nBlockMax          — 每进程接收的最大数据量
```

这种优化的关键洞察：位移量始终等于 `MaxBlock`，因此无需每次重新计算，这减少了 MPI 内部的开销。

---

## 3. 通信模式

### 3.1 Ghost Cell 交换：message_pass_cell

核心通信函数是 BATL 提供的 `message_pass_cell`，封装在 `ModMessagePass.f90` 的 `exchange_messages` 子程序中。

**三种通信模式**（由 `TypeMessagePass` 参数控制）：

| 模式 | 策略 | 适用场景 |
|------|------|---------|
| `'all'` | 一次性传递所有 ghost cell 层 | 所有阶数均适用（简化实现） |
| `'normal'` (默认) | 分阶段：先传面、再传边和角 | 更高效，减少冗余通信 |
| `'reschangeonly'` | 仅传递 AMR 层级变化处的 ghost cell | 网格未变化时的优化 |

### 3.2 通信流程

```
exchange_messages 执行流程:
┌──────────────────────────────────────────────────┐
│ 1. sync_cpu_gpu('update on GPU')                 │ ← OpenACC 数据同步
│ 2. 应用物理边界条件 (set_cell_boundary)           │
│    - 如果 iTypeUpdate == UpdateFast_:            │
│      set_boundary_fast()                         │
│ 3. [可选] 周期楔形坐标变换 (IsPeriodicWedge)       │
│ 4. message_pass_cell (核心通信)                   │
│    ├─ 非阻塞 MPI 点对点通信                      │
│    ├─ Ghost cell 数据打包/解包                   │
│    └─ AMR 层级间的限制/延拓                      │
│ 5. fix_boundary_ghost_cells                      │
│ 6. [可选] 逆周期楔形变换                           │
│ 7. 重新应用边界条件 (确保角 ghost cell 正确)      │
│ 8. limit_pressure (物理约束)                     │
│ 9. fill_in_from_buffer (球形缓冲网格)            │
└──────────────────────────────────────────────────┘
```

### 3.3 Ghost Cell 层级处理

由于 AMR 允许相邻块具有不同的细化级别，ghost cell 填充需要特殊处理：

- **同级别邻居** (`DiLevel_EB = 0`)：直接拷贝
- **粗 → 细** (`DiLevel_EB = -1`)：使用**延拓**（prolongation），粗网格值插值到细 ghost cell
  - `nOrderProlong=1`：一阶延拓（直接拷贝，配合 face 处限制）
  - `nOrderProlong=2`：二阶延拓（更精确的插值）
- **细 → 粗** (`DiLevel_EB = +1`)：使用**限制**（restriction），细网格值平均或注入到粗 ghost cell
  - `DoRestrictFace` 控制 face 处的限制操作

### 3.4 缓冲网格（Buffer Grid）

`ModBuffer.f90` 实现了一个球形壳层缓冲网格，用于处理外部边界条件（如太阳风输入）：

- 在球形坐标系的固定半径范围（如 19-21 Rs）上定义二维经纬度网格
- 存储全局一致的边界条件状态（`BufferState_VG`）
- 支持坐标系统转换（`DoTransformCoord`）
- 通过 `fill_in_from_buffer` 将缓冲网格数据插值到计算块的 ghost cell 中
- 使用 OpenACC 指令（19 个 `!$acc`）实现 GPU 加速

### 3.5 全局通信模式总结

| 通信类型 | MPI 函数 | 用途 |
|---------|---------|------|
| 点对点 | 非阻塞 MPI（通过 BATL 封装） | Ghost cell 交换 |
| 全局收集 | `MPI_allgatherv` | 收集各进程块类型信息（`iTypeAdvance_BP`） |
| 全局归约 | `MPI_allreduce` | 确定 `nBlockMax`、`nBlockALL`、负载均衡权重 |
| 全局求和 | `MPI_allreduce(MPI_IN_PLACE, ...)` | 累加跨进程的 `iTypeBalance_A` |

---

## 4. 负载均衡算法

### 4.1 空间填充曲线排序

BATSRUS 使用**空间填充曲线**（Space-Filling Curve, SFC）将多维块映射到一维排序：

- **主要选择：Morton 排序**（Z-order curve）
  - 优点：计算快速（位交错操作）、较好的局部性
  - 论文指出："This ordering corresponds to the Morton space filling curve"
  
- **备选：Peano-Hilbert 曲线**
  - 原文记载："Another popular choice that we actually use in the original BATS-R-US code is the Peano-Hilbert space filling curve"
  - Morton 在当前实现中是默认选择

**负载均衡流程**：

```
┌─────────────────────────────────────────────────┐
│            load_balance 核心流程                  │
├─────────────────────────────────────────────────┤
│ 1. select_stepping: 确定块类型                    │
│    ├─ 隐式/显式/稳态/跳过的分类                   │
│    └─ MPI_allgatherv 收集全局类型分布              │
│                                                   │
│ 2. 块类型分类（多比特位编码）                      │
│    ├─ Bit 0: SkippedBlock（跳过块=0）             │
│    ├─ Bit 1: IsNoBody_B（真块）                   │
│    ├─ Bit 2: iPartImplBlock（部分隐式）            │
│    ├─ Bit 3: iSemiImplBlock（半隐式）              │
│    ├─ Bit 4: iPointImplBlock（点隐式）             │
│    ├─ Bit 5: iFieldLineThreadBlock（场线线程）      │
│    ├─ Bit 6: iPicBlock（PIC 区域）                │
│    ├─ Bit 7: iSteadyBlock（稳态块）                │
│    ├─ Bit 8: iHighOrderBlock（高阶格式）           │
│    ├─ Bit 9+: iUserTypeBlock（用户自定义）          │
│    └─ 高位：iSubCycleBlock（子循环时间级）          │
│                                                   │
│ 3. MPI_allreduce 全局累加类型权重                  │
│                                                   │
│ 4. 统计实际存在的块类型（IsTypeExist_I）            │
│    └─ 建立间接索引 iType_I: 可能类型 → 实际类型     │
│                                                   │
│ 5. regrid_batl: 核心重网格                         │
│    ├─ 按 Morton 曲线排序节点                       │
│    ├─ 按权重将节点分配给进程                        │
│    ├─ Restrict: 粗化 → 数据限制                    │
│    ├─ Load Balance: 重新分配块                     │
│    └─ Prolong: 细化 → 数据延拓                    │
│                                                   │
│ 6. pack/unpack_load_balance: 迁移块数据            │
│    ├─ 解变量 (State_VGB)                          │
│    ├─ ∇·B 约束面数据 (BxFace, ByFace, BzFace)     │
│    ├─ 隐式历史状态 (ImplOld_VCB)                  │
│    ├─ 场线追踪数据 (Trace_DSNB)                   │
│    └─ 用户自定义块数据 (block_data)                │
│                                                   │
│ 7. set_batsrus_grid: 重建邻居关系                  │
│    └─ 更新 DiLevel_EB, jProc_IEB, jBlock_IEB      │
│                                                   │
│ 8. sync_cpu_gpu_amr: GPU 数据同步                 │
└─────────────────────────────────────────────────┘
```

### 4.2 权重基础负载均衡

负载均衡的核心是**带权重的块分区**：

- 每种块类型被视为相同的计算权重
- 不同类型可能具有不同的权重（通过 `iTypeBalance_A` 编码）
- 使用 `iTypeAdvance_A` 区分被跳过/未使用的块
- 支持 `DoBalanceEachLevel`：是否在每个 AMR 层级内单独均衡（当使用局部时间步长时）

### 4.3 动态负载重均衡

`load_balance_blocks` 子程序在以下情况下触发动态重均衡：

| 条件 | 触发原因 |
|------|---------|
| `UseMaxTimeStep` | 子循环方案中的时间级变化 |
| `UsePartImplicit` | 隐式/显式块类型重新分配 |
| `UsePartSteady ∧ IsNewSteadySelect` | 稳态选择变化 |
| `nBlockSemi ≥ 0 ∧ DoBalanceSemiImpl` | 半隐式区域变化 |
| `DoBalancePointImplicit ∧ nIteration>1` | 点隐式动态变化 |
| `UsePic ∧ DoBalanceActivePicBlock` | PIC 活动区域变化 |

`IsDynamicSemiImpl` 和 `IsDynamicPointImplicit` 控制是否持续进行重均衡。

### 4.4 数据迁移机制

块在不同进程间迁移时，需要安全地传递所有相关数据（`pack_load_balance` / `unpack_load_balance`）：

```
Buffer_I 格式：
┌──────────────────────────────────────────────┐
│ [1] nDynamicData (标量)                      │
│ [2] UseUserPointImplicit_B (标量)            │
│ [3..] Bx/By/BzFace_GB (if UseConstrainB)    │
│ [...] B0_DGB (if DoMoveExtraData)            │
│ [...] ImplOld_VCB (if UseImplicit+BDF2)      │
│ [...] Trace_DSNB (if DoSendTrace)            │
│ [...] block_data (if nDynamicData > 0)       │
└──────────────────────────────────────────────┘
```

缓冲区大小 `nBuffer` 在 `init_load_balance` 中预先计算，以容纳最大可能的数据量。

---

## 5. GPU 加速：OpenACC

### 5.1 加速方案选择

与许多其他 HPC 代码采用 CUDA 不同，BATSRUS 选择了 **OpenACC** 指令式编程模型：

| 属性 | 值 |
|------|-----|
| API | OpenACC（非 CUDA） |
| 编译器 | nvfortran（NVIDIA HPC SDK） |
| 编译配置 | `./Config.pl -install -compiler=nvfortran,nvc` → `./Config.pl -acc` |
| 关键指令前缀 | `!$acc` |
| 重写代码量 | 约1%（核心求解器），其余通过 OpenACC 指令自动加速 |
| 参考系统 | Frontera (nvhpc-hpcx/24.5), Pleiades (nvhpc-nompi/24.3 + mpi-hpe/mpt) |

### 5.2 OpenACC 指令分布

| 源文件 | `!$acc` 指令数 | 主要内容 |
|--------|---------------|---------|
| `ModBuffer.f90` | 19 | 缓冲网格的 GPU 计算 |
| `ModMessagePass.f90` | 13 | Ghost cell 交换的 GPU 同步 |
| `ModParallel.f90` | 2 | 邻居信息数组的 GPU 声明 |
| `ModBatlInterface.f90` | 2 | `update device` 同步 |
| `ModLoadBalance.f90` | 1 | AMR 后 GPU 状态更新 |
| `ModBatsrusUtility.f90` | 1 | `routine seq` 声明 |

### 5.3 GPU 数据管理

BATSRUS 通过以下机制管理 CPU/GPU 数据一致性：

```fortran
! 在 GPU 上声明关键数组
!$acc declare create(DiLevel_EB)
!$acc declare create(jProc_IEB, jBlock_IEB)

! 数据同步函数
call sync_cpu_gpu('update on GPU', NameSub, State_VGB)
call sync_cpu_gpu_amr              ! AMR 操作后的同步
!$acc update device(DtMax_B)       ! 特定数组的更新
```

**关键模式**：在每次 `message_pass_cell` 调用前进行 `sync_cpu_gpu`，确保 GPU 上的数据最新；在负载均衡后调用 `sync_cpu_gpu_amr` 同步所有 AMR 相关数据。

### 5.4 GPU 通信架构

GPU 论文（arXiv:2501.06717, 2025）详细描述了多 GPU 通信方案：

> "To port it to multiple GPUs, we implement a new message passing algorithm to support its unique block-adaptive grid feature."

这意味着 **GPU 间的 MPI 通信逻辑被重写**，以适应 GPU 直接内存访问（GPU-Direct RDMA）和非连续数据传输的需求。

**多 GPU 通信模式**：
- 节点内 4 GPU：NVLink 直连 → 95% 并行效率
- 跨节点：MPI + GPU-Direct → 50-60% 效率（受硬件瓶颈限制）

### 5.5 GPU 性能数据

| 指标 | 数值 |
|------|------|
| 单 A100 vs CPU | 1 GPU ≈ 270 AMD "Rome" CPU 核心 (2.1 个 128 核节点) |
| 单 GPU 实时比 | 3.6× faster than real-time |
| 4 GPU 实时比 | 6.9× faster than real-time |
| 256 GPU 弱扩展 | 50-60% 并行效率 |
| 单节点 (4 GPU) | 95% 效率 |
| 总代码量 | 200,000+ 行 Fortran，GPU 专用改写仅 ~1% |

---

## 6. 共享内存并行：OpenMP

### 6.1 混合 MPI+OpenMP 架构

BATSRUS 实现了 **MPI + OpenMP 混合并行**（参见 Sciencedirect 2020 论文）：

- **MPI**：跨节点的分布式内存并行
- **OpenMP**：节点内的共享内存并行

主要优势：**减少每进程内存占用**。在大规模模拟中，纯 MPI 方案因内存限制无法使用更多网格单元，而混合方案允许更少的 MPI 进程持有更大的网格。

### 6.2 OpenMP 应用场景

OpenMP 在 BATSRUS 中的主要应用：

1. **ModThreadedLC.f90**：场线线程化边界条件
   ```fortran
   !$omp parallel do
   do iBlock = 1, nBlock
       ! 每个线程处理一个块
   end do
   !$omp end parallel do
   ```

2. **exchange_messages** 中的边界条件应用
   ```fortran
   !$omp parallel do
   do iBlock = 1, nBlock
       if (Unused_B(iBlock)) CYCLE
       call set_cell_boundary(...)
       call fill_in_from_buffer(iBlock)
       call limit_pressure(...)
   end do
   !$omp end parallel do
   ```

3. **线程私有数据**
   ```fortran
   real, allocatable :: Te_G(:,:,:)
   !$omp threadprivate(Te_G)
   ```
   避免线程间的数据竞争。

### 6.3 线程安全

- `!$omp threadprivate`：为每个线程创建独立的温度数组副本
- 块级操作天然适合 `!$omp parallel do`，因为块之间无数据依赖
- 边界通信（MPI）仍然是进程级别的操作，不在 OpenMP 区域内执行

---

## 7. 可扩展性

### 7.1 历史基准

| 系统 | 处理器数 | 性能 |
|------|---------|------|
| Cray T3E-1200 | 1,490 PEs | 342 GFlops 持续 |
| Cray T3E（早期） | ~1,500 处理器 | 近乎完美线性扩展 |
| 多种平台 | — | "nearly perfectly scales to 1,500 processors" |
| IBM SP2, SGI Powerchallenge | — | 完整移植支持 |

来源：SC98 论文、NASA 技术报告

### 7.2 现代 CPU 扩展

- **MPI+OpenMP 混合**（2020）：显著缓解内存限制
- 每核心内存需求降低（更少的 MPI 进程）
- 适合现代多核节点（每节点 64-128 核心）

### 7.3 GPU 扩展

| 规模 | 弱扩展效率 | 瓶颈 |
|------|-----------|------|
| 1 节点 (≤4 GPU) | 95% | — |
| 多节点 (≤256 GPU) | 50-60% | 硬件通信瓶颈（跨节点 NVLink 带宽不足） |
| 256 GPU + | 良好（小问题）或效率下降（大问题） | 问题规模与通信比的权衡 |

---

## 8. 代码架构总结

### 8.1 关键文件功能矩阵

| 文件 | 核心功能 | MPI | OpenACC | OpenMP |
|------|---------|:---:|:-------:|:------:|
| `ModParallel.f90` | 邻居关系数据结构 | ✓ | ✓ (2) | — |
| `ModMessagePass.f90` | Ghost cell 交换主控 | ✓ | ✓ (13) | ✓ |
| `ModLoadBalance.f90` | 负载均衡与块类型分类 | ✓ | ✓ (1) | — |
| `ModBuffer.f90` | 球形缓冲网格 | — | ✓ (19) | — |
| `ModBatlInterface.f90` | BATL ⇄ BATSRUS 桥接 | ✓ | ✓ (2) | — |
| `ModBatsrusUtility.f90` | 工具函数/GPU 同步 | — | ✓ (1) | — |
| `ModThreadedLC.f90` | 场线线程化边界 | — | — | ✓ |
| `ModFieldLineThread.f90` | 场线追踪 | — | — | ✓ |

### 8.2 依赖关系

```
ModMain
  ├── BATL_lib (外部库)
  │     ├── BATL_tree (八叉树管理)
  │     ├── message_pass_cell (核心 MPI 通信)
  │     ├── regrid_batl (重网格+负载均衡)
  │     └── find_test_cell (测试辅助)
  │
  ├── ModParallel (邻居关系)
  │     └── 由 set_batsrus_grid 填充
  │
  ├── ModMessagePass (Ghost cell 交换)
  │     ├── 调用 message_pass_cell
  │     ├── 调用 ModCellBoundary
  │     └── 调用 ModBuffer
  │
  ├── ModLoadBalance (负载均衡)
  │     ├── 调用 regrid_batl
  │     ├── 调用 ModBatlInterface
  │     └── pack/unpack_load_balance
  │
  └── ModBatlInterface (BATL 接口)
        └── set_batsrus_grid / set_batsrus_state
```

### 8.3 设计原则总结

1. **关注点分离**：BATL 库处理 AMR、消息传递和负载均衡的底层实现；BATSRUS 专注于物理求解
2. **渐进式并行**：MPI → MPI+OpenMP → MPI+OpenMP+OpenACC，每一层都是可选的增量加速
3. **块为基本单元**：所有操作（通信、负载均衡、细化）都以块为原子单位
4. **空间局部性优先**：空间填充曲线排序最大化数据局部性，减少通信
5. **类型感知负载均衡**：不仅仅是计数块，而是根据块的计算特性（隐式/显式/PIC/稳态等）分配权重

---

## 参考资料

1. Tóth, G., et al. (2012). "Adaptive numerical algorithms in space weather modeling." *J. Comp. Phys.*, 231, 870–903.
2. BATSRUS GPU 论文 (2025). "Faster-than-Real-Time Magnetospheric Simulations with a Block-Adaptive Grid Code." arXiv:2501.06717, *ApJ*, 981, 188.
3. SC98 论文: "Adaptive Parallel Computation of a Grand-Challenge Problem: Prediction of the Path of a Solar Coronal Mass Ejection."
4. SWMF 官方文档: https://ccmc.gsfc.nasa.gov/models/SWMF~GM=BATSRUS~20180525/
5. BATSRUS GitHub: https://github.com/SWMFsoftware/BATSRUS
6. BATL GitHub: https://github.com/SWMFsoftware/BATL
7. Hybrid MPI+OpenMP 论文 (2020), *J. Parallel Distrib. Comput.*
