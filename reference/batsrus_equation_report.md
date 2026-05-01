# BATSRUS 物理方程模块全景分析报告

> 生成者：batsrus_equation（deepseek_pro 之化身）
> 生成时间：2026-05-01
> 数据源：https://github.com/SWMFsoftware/BATSRUS
> 分析目录：`srcEquation/` (75 模块) · `srcUser/` (33 模块) · `srcInterface/` (14 文件)

---

## 一、总览

BATSRUS (Block-Adaptive Tree Solar-wind Roe Upwind Scheme) 是密歇根大学空间天气建模框架 (SWMF) 的核心 MHD 求解器。它采用 **有限体积法 + Roe 近似黎曼求解器 + 块自适应网格 (AMR)**，支持从理想 MHD 到多流体、多离子、辐射输运等多种物理模型。方程模块位于 `srcEquation/` 目录，每个模块定义一组守恒变量及对应的特征波结构。

### 统计摘要

| 指标 | 数值 |
|------|------|
| 方程模块总数 | **75** |
| 用户模块（应用案例） | **33** |
| SWMF 耦合接口 | **9**（另有 5 个辅助文件） |
| 最小 nVar | 1（ModEquationScalar） |
| 最大 nVar | 40（ModEquationOuterHelio2PUIPe） |
| 典型 nVar 范围 | 5-29 |
| 主要物理模型 | 理想 MHD、多流体、霍尔 MHD、电子压强、辐射、外日球层、彗星、行星、多矩 |

---

## 二、方程模块分类

### 2.1 核心 MHD 模块（理想 MHD + 扩展）

最基础的理想 MHD 方程模块为 `ModEquationMhd.f90`（nVar=8），采用 8 波结构（1 熵波 + 2 Alfvén 波 + 4 磁声波 + 1 散度波）。所有 MHD 衍生模块皆基于此扩展。

| 模块 | nVar | nWave | 类型 | 特性 |
|------|------|-------|------|------|
| **ModEquationMhd** | 8 | ? | 理想 MHD | 基准模块，8 波结构 |
| ModEquationMhdHyp | 9 | ? | MHD + 双曲散度清理 | hyperbolic divergence cleaning |
| ModEquationMhdPe | 9 | ? | MHD + 电子压强 | separate electron pressure |
| ModEquationMhdAnisoP | 9 | ? | MHD + 各向异性压强 | anisotropic pressure |
| ModEquationMhdEos | 9 | ? | MHD + 通用物态方程 | general EOS |
| ModEquationMhdSA | 9 | ? | MHD + 流线对齐 | stream-aligned |
| ModEquationMhdNonCons | 8 | ? | 非守恒 MHD | non-conservative formulation |
| ModEquationMhdEosRad | 9 | 1 | MHD + EOS + 辐射 | radiation |
| ModEquationMhdHpOp | 10 | ? | 双离子 MHD（H⁺+O⁺） | 2 Species, Tóth 2017 |
| ModEquationMhdHypPe | 10 | ? | 双曲 MHD + 电子压强 | hyperbolic + Pe |
| ModEquationMhdPeSA | 10 | ? | MHD + Pe + 流线对齐 | |
| ModEquationMhdPeAnisoPi | 10 | ? | MHD + 各向异性 Pi | anisotropic Pi |
| ModEquationMhdPeAniso | 11 | ? | MHD + 各向异性 Pe | |
| ModEquationMhdSwIono | 10 | ? | 太阳风-电离层双离子 | 2 Species (SwIono) |
| ModEquationMhdFluidsPeAniso | 17 | ? | 双流体 MHD + 各向异性 | 2-fluid with aniso ion/electron pressures |

### 2.2 MHD 波模块

专门用于 Alfvén 波追踪的方程模块：

| 模块 | nVar | nWave | 特性 |
|------|------|-------|------|
| ModEquationMhdWaves | 9 | 2 | MHD + Alfvén 波 |
| ModEquationMhdWavesHyp | 9 | 2 | + 双曲散度清理 |
| ModEquationMhdWavesHypPe | 10 | 2 | + Pe + 双曲散度清理 |

### 2.3 霍尔 MHD / 多离子模块

| 模块 | nVar | 类型 | 特性 |
|------|------|------|------|
| ModEquationMultiIon | 13 | 多离子 MHD | Multi-ion MHD |
| ModEquationMultiIonPe | 14 | 多离子 + Pe | H⁺, O⁺ and Pe |
| ModEquationSwh | 9 | 太阳风加热 | Solar wind protons with extra indexes |
| ModEquationSwhPui | 14 | SWH + PUI | Multi-ion MHD |
| ModEquationSwhPuiPe | 15 | SWH + PUI + Pe | |
| ModEquationSwh2PuiPe | 20 | SWH + 2 PUI + Pe | |
| ModEquationIonsHNO | 18 | 电离层离子 | H, N, O ions |
| ModEquationMhdHd | 13 | MHD + HD 耦合 | MHD and HD |

### 2.4 纯流体力学 (HD) 模块

| 模块 | nVar | nWave | 特性 |
|------|------|-------|------|
| ModEquationHd | 5 | ? | 纯 Euler HD（密度+动量+能量） |
| ModEquationHdEos | 6 | ? | HD + 通用 EOS |
| ModEquationHdEosRad | 6 | 1 | HD + 辐射 |
| ModEquationHdCrash | 6 | ? | HD + 电离能级 |
| ModEquationHdRadCrash | 6 | 1 | HD + 电离 + 辐射 |

### 2.5 CRASH 辐射模块

CRASH (Center for Radiative Shock Hydrodynamics) 相关模块，耦合辐射输运：

| 模块 | nVar | nWave | 特性 |
|------|------|-------|------|
| ModEquationCrash | 7 | 1 | HD + 电离 + 能级 + 电子能量 + 辐射 |
| ModEquationCrashTe | 8 | 1 | + 分离电子温度 |
| ModEquationMhdCrash | 11 | 1 | MHD + 电离 + 能级 + 电子能量 + 辐射 |

### 2.6 AWSoM（日冕/太阳风）模块

Alfvén Wave Solar Model (AWSoM)，由 van der Holst 等人开发，用于日冕与太阳风模拟：

| 模块 | nVar | nWave | 特性 |
|------|------|-------|------|
| ModEquationAwsom | 10 | 2 | 基准 AWSoM（MHD + Alfvén 波 + Pe） |
| ModEquationAwsomSA | 11 | 2 | + 流线对齐 |
| ModEquationAwsomWdiff | 11 | 2 | + 波扩散 |
| ModEquationAwsomSAWdiff | 12 | 2 | + SA + Wdiff |
| ModEquationAwsomAnisoPi | 11 | 2 | + 各向异性 Pi |
| ModEquationAwsomAnisoPiSA | 12 | 2 | + 各向异性 Pi + SA |
| ModEquationAwsomChargeState | 11 | 2 | + 电荷态 |
| ModEquationAwsomFluids | 10 | 2 | 多流体 AWSoM |
| ModEquationAwsomFluidsAnisoPi | ? | 2 | 多流体 + 各向异性 Pi |

### 2.7 外日球层模块（最大 nVar）

携带多个中性流体与拾取离子 (PUI) 的超大方程组：

| 模块 | nVar | nWave | 特性 |
|------|------|-------|------|
| ModEquationOuterHelio | 29 | ? | MHD + 4 中性流体 + level-set |
| ModEquationOuterHelio2d | 28 | ? | 2D 简化版 |
| ModEquationOuterHelioPUI | 34 | ? | SWH + PUI + 4 中性流体 |
| ModEquationOuterHelioPUIPe | 35 | ? | SWH + PUI + Pe + 4 中性流体 |
| ModEquationOuterHelioPuiBin | 35 | ? | + PUI 分 bin |
| **ModEquationOuterHelioAwsom** | **36** | **2** | AWSoM + PUI + 4 中性流体 |
| ModEquationOuterHelioAwsomPuiBin | 36 | 2 | AWSoM + PUI bin + 4 中性流体 |
| ModEquationOuterHelio2PUIPe | **40** | ? | SWH + 2 PUI + Pe + 4 中性流体（最大） |

### 2.8 行星模块

| 模块 | nVar | 类型 | 特性 |
|------|------|------|------|
| ModEquationMhdMars | 12 | Mars | MHD 火星 |
| ModEquationMhdMarsPe | 13 | Mars | + Pe |
| ModEquationMarsFluids | 23 | Mars | 多流体火星 |
| ModEquationMarsFluidsPe | 23 | Mars | 多流体 + Pe |
| ModEquationMarsFluidsSw | **28** | Mars | 多流体 + 太阳风 |
| ModEquationMhdTitan | 15 | Titan | MHD 土卫六 |
| ModEquationMhdTitanPe | 16 | Titan | + Pe |
| ModEquationMhdSaturn3sp | 11 | Saturn | 土星 3 物种 |
| ModEquationMhdEuropa2FluidsPe | 14 | Europa | 木卫二 2 流体 |
| ModEquationMhdEuropa3FluidsPe | 19 | Europa | 木卫二 3 流体 |
| ModEquationMhdComet | 14 | Comet | 彗星 MHD |
| ModEquationComet3FluidsPe | 19 | Comet | 3 流体彗星 |
| ModEquationCometCG3FluidsPe | 19 | Comet | 2 流体 + Pe + 1 中性 |
| ModEquationCometCG3FluidsPeHyp | 20 | Comet | + 双曲清理 |
| ModEquationMhdComet1Sp | 8 | Comet | 单种彗星 |

### 2.9 极风 / 其他模块

| 模块 | nVar | 特性 |
|------|------|------|
| ModEquationMhdPw | 11 | 极风边界 3 物种 (H⁺, O⁺, He⁺) |
| ModEquationMhdCorona | 9 | 太阳日冕 MHD |
| ModEquationScalar | 1 | 标量方程（测试用） |

### 2.10 多矩闭合模块

高阶矩闭合模块（超出理想 MHD 的五矩/六矩方程）：

| 模块 | nVar | 特性 |
|------|------|------|
| ModEquationFiveMoment | 17 | 五矩闭合 |
| ModEquationFiveMomentHyp | 18 | 五矩 + 双曲清理 |
| ModEquationSixMoment | 19 | 六矩闭合 + 各向异性压强 |
| ModEquationSixMomentHyp | 20 | 六矩 + 双曲散度 B 清理 |

### 2.11 多离子流体模块

| 模块 | nVar |
|------|------|
| ModEquationThreeIonFluidPe | 19 |
| ModEquationFourIonFluidPe | 24 |

---

## 三、方程参数分布统计

### 3.1 nVar 分布

| nVar 范围 | 模块数 | 占比 | 典型代表 |
|-----------|--------|------|----------|
| 1-5 | 2 | 2.7% | Scalar(1), HD(5) |
| 6-10 | 29 | 38.7% | 标准 MHD 及其基本扩展 |
| 11-15 | 15 | 20.0% | 多离子/彗星/极风 MHD |
| 16-20 | 15 | 20.0% | 多矩/多流体 MHD |
| 21-30 | 7 | 9.3% | 火星多流体 / 外日球层 |
| 31-40 | 6 | 8.0% | 外日球层超大方程组 |

**统计**: 约 60% 的模块 nVar ≤ 15，属于中小规模方程组；外日球层模块的 nVar 可达 28-40，远超其他类别。

### 3.2 nWave 分布

| nWave | 模块数 | 说明 |
|-------|--------|------|
| 1 | 7 | 纯 HD + 辐射模块 |
| 2 | 15 | AWSoM 及波追踪模块 |
| 未提取 | 53 | 未在头部显式定义（运行时动态确定） |

nWave 多数情况下未在模块头部显式定义——它由方程系统的特征结构动态确定。

### 3.3 物理模型谱系

```
理想 MHD (nVar=8)
├── + 双曲散度清理 (nVar=9) ─→ MhdHyp
├── + 电子压强分离 (nVar=9) ─→ MhdPe
├── + 各向异性压强 (nVar=9-11) ─→ MhdAnisoP, MhdPeAniso
├── + 通用 EOS (nVar=9) ─→ MhdEos
├── + 流线对齐 (nVar=9-10) ─→ MhdSA, MhdPeSA
├── + Alfvén 波 (nVar=9-10) ─→ MhdWaves → AWSoM 家族
├── + 多离子 (nVar=13-20) ─→ MultiIon, Swh, MhdHpOp
├── + 多流体 (nVar=17-28) ─→ MhdFluidsPeAniso, MarsFluids
├── + 辐射 (nVar=9-11) ─→ MhdEosRad, MhdCrash
├── + 多矩闭合 (nVar=17-20) ─→ FiveMoment, SixMoment
└── + 外日球层全耦合 (nVar=29-40) ─→ OuterHelio 家族

纯 HD (nVar=5-6)
├── + 辐射 (nVar=6) ─→ HdEosRad, HdRadCrash
└── + CRASH 电离 (nVar=7-8) ─→ Crash, CrashTe

标量 (nVar=1)
```

---

## 四、用户模块 (srcUser/) —— 33 个应用案例

| 模块 | 应用领域 | 备注 |
|------|----------|------|
| ModUserAwsom | 日冕/太阳风 | AWSoM 驱动 |
| ModUserCCMC | CCMC 运行 | 最大文件 (234 KB) |
| ModUserComet1Sp | 彗星 1 物种 | |
| ModUserComet3FluidsPe | 彗星 3 流体 | |
| ModUserComet6Sp | 彗星 6 物种 | |
| ModUserCometCG | 彗星 CG | |
| ModUserCometCGfluids | 彗星 CG 多流体 | |
| ModUserCometNeutralFluids | 彗星中性流体 | |
| ModUserCylOutflow | 柱状流出 | |
| ModUserEarthXray | 地球 X 射线 | |
| ModUserEe | 电动力学 | |
| ModUserFluxRope | 磁通量绳 | |
| ModUserGemReconnect | GEM 重联 | |
| ModUserJupiter | 木星 | |
| ModUserKelvinHelmholtz | K-H 不稳定性 | 最小模块 (3.8 KB) |
| ModUserLogAdvection | 对数平流 | |
| ModUserMars | 火星 | |
| ModUserMarsFluids | 火星多流体 | |
| ModUserMercury | 水星 | |
| ModUserMoonImpact | 月球撞击 | |
| ModUserNonGyro | 非回旋 | |
| ModUserOuterHelio | 外日球层 | 最大模块 (187 KB) |
| ModUserOuterHelio2d | 外日球层 2D | |
| ModUserPointImplicit | 点隐式 | |
| ModUserSaturn | 土星 | |
| ModUserStretchedDipole | 拉伸偶极场 | |
| ModUserSwIono | 太阳风-电离层 | |
| ModUserSwarm | 群 | |
| ModUserSwitchback | 折返 | |
| ModUserTitan | 土卫六 | |
| ModUserVenus | 金星 | |
| ModUserWaveReflection | 波反射 | |
| ModUserWaves | 波 | |

用户模块覆盖 **太阳系内几乎所有行星及卫星**（水星、金星、地球、火星、木星、土星、土卫六、木卫二、彗星）+ 太阳（日冕、外日球层）+ 基础等离子体过程（重联、K-H 不稳定性、波反射等）。

---

## 五、SWMF 耦合接口 (srcInterface/)

BATSRUS 作为 SWMF 的核心组件 (GM - Global Magnetosphere)，通过标准耦合接口与其他物理模块交互：

| 耦合接口 | 耦合目标 | 文件大小 |
|----------|----------|----------|
| GM_couple_ie | **IE** — 电离层电动力学 (Ridley Ionosphere) | 14.6 KB |
| GM_couple_im | **IM** — 内磁层 (Rice Convection Model) | 61.0 KB |
| GM_couple_pc | **PC** — 粒子编码 / PIC (MHD-EPIC) | 16.0 KB |
| GM_couple_ps | **PS** — 极风流出 (Polar Wind Outflow) | 6.0 KB |
| GM_couple_pt | **PT** — 粒子追踪 | 2.8 KB |
| GM_couple_pw | **PW** — 极风 (PWOM) | 18.7 KB |
| GM_couple_rb | **RB** — 辐射带 | 15.6 KB |
| GM_couple_ua | **UA** — 上层大气 (GITM) | 5.2 KB |
| GM_couple_ih | **IH** — 内日球层 | 13.6 KB |
| GM_wrapper | 包装器 / 路由器 | 16.9 KB |
| ModGridDescriptor | 网格描述符 | 9.6 KB |

### 耦合架构

```
                    ┌─────────────┐
                    │  IH (内日球层) │ ←── 太阳风输入
                    └──────┬──────┘
                           │
  ┌────────────┐    ┌─────┴─────┐    ┌──────────────┐
  │ PS (极风流出) │◄───│           │───►│ PT (粒子追踪)  │
  └────────────┘    │           │    └──────────────┘
                    │  GM/BATSRUS │
  ┌────────────┐    │  (全局磁层) │    ┌──────────────┐
  │ IE (电离层)  │◄───│           │───►│ IM (内磁层)    │
  └────────────┘    └─────┬─────┘    └──────────────┘
                           │
  ┌────────────┐    ┌─────┴─────┐    ┌──────────────┐
  │ UA (上层大气) │◄───│ PC (PIC)   │───►│ RB (辐射带)   │
  └────────────┘    └────────────┘    └──────────────┘
```

GM_wrapper 负责将内外边界条件、场线追踪结果等在不同物理组件间路由。

---

## 六、架构洞察

### 6.1 模块化设计原则

1. **一物理一模块**：每个物理模型（方程集合）对应一个 ModEquation 模块，定义变量数、波结构、通量函数、特征速度等
2. **继承链**：复杂模块继承简单模块的变量集合（如 MhdPe 在 Mhd 的 8 变量上增加电子压强）
3. **用户模块 = 方程模块 + 边界/初始条件 + 源项**：用户模块指定使用哪个方程模块，并提供具体场景的设置

### 6.2 变量增长路径

```
理想 MHD (8)
  ├ + 电子压强 → 9
  ├ + 双曲散度 → 9
  ├ + Alfvén 波 → 9-10
  ├ + 辐射 → 9-11
  ├ + 多离子 → 13-20
  ├ + 多流体 → 17-28
  ├ + 彗星中性流体 → 19-20
  ├ + 外日球层中性流体 → 29-36
  └ + 多 PUI + 多中性流体 → 35-40
```

### 6.3 计算复杂度

- **nVar × 网格点数** 决定状态向量存储
- **nWave 特征波数** 决定黎曼求解器复杂度（每个界面需解 nWave 个特征方程）
- AWSoM 家族 (nWave=2) 使用 Alfvén 波分解，舍去完整的磁声波特征结构——这大大降低了黎曼求解器开销
- 外日球层模块虽 nVar 极大 (28-40)，但耦合多个中性流体与 PUI，物理上丰富

---

## 七、结论

BATSRUS 的方程模块体系展现了 **"一套框架、万般物理"** 的设计哲学：

1. **75 个方程模块** 覆盖从标量输送到 40 变量多流体+多中性流体+多拾取离子的全耦合系统的完整谱系
2. **33 个用户模块** 实现太阳系内所有主要行星/卫星及基础等离子体过程的模拟
3. **9 个耦合接口** 使 BATSRUS 可无缝嵌入 SWMF 多物理框架
4. **AWSoM 子家族 (nWave=2)** 是天体物理模拟中的特色设计——以 Alfvén 波分解换取计算效率
5. **外日球层模块** 是最大的方程组（nVar ≤ 40），代表了日球层-星际介质相互作用的最高保真度模型

---

*报告完毕。*
