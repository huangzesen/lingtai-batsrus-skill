# BATSRUS 案例目录

> 从 `Param/` 下 36 个案例目录中精选的代表性 PARAM.in 分析
> 运行：`make test_<case>` 或 `make -j test_<case> NP=4`

---

## 一、标准验证案例

### 1. Brio-Wu 激波管 (`Param/SHOCKTUBE/`)

**物理场景**：MHD 激波管问题，含复合激波、稀疏波、接触间断和慢激波

| 参数 | 值 |
|------|-----|
| 方程模块 | 理想 MHD（默认 8 变量） |
| 网格 | 单块 256×4×4，沿 Z 轴旋转 |
| 时间推进 | nStage=2（Dt/2-Dt 二阶 RK） |
| CFL | 0.8 |
| 通量格式 | Roe（Roe 型黎曼求解器） |
| 限制器 | mc3（monotonized-centered），β=1.5 |
| 精度 | nOrder=2 |
| 边界 | X: shear, Y: shear, Z: float |
| 终止 | tSimulationMax=25.6 |
| 运行 | `make test_shocktube` |

### 2. Orszag-Tang 涡旋

**物理场景**：2D MHD 湍流过渡——激波-激波相互作用

| 参数 | 值 |
|------|-----|
| 方程模块 | 理想 MHD |
| 特点 | 标准 MHD 代码验证基准 |
| 测试 | 激波捕捉、散度控制 |
| 参考文献 | Orszag & Tang, JFM 1979 |

---

## 二、行星磁层案例

### 3. 地球磁层 (`Param/EARTH/`)

**物理场景**：太阳风-地球磁层相互作用，准稳态

| 参数 | 值 |
|------|-----|
| 根块 | 2×1×1 |
| 网格范围 | x:[-224,32] y:[-64,64] z:[-64,64] |
| AMR 判据 | 球体 AMR（r=5.25, 17.5）+ 多个 box 嵌套 |
| 最小分辨率 | 1/8（内磁层外层） |
| 起始时间 | 1998-05-01 00:00 |
| 时间推进 | 非时间精确（定常） |
| 保存 | 每 1000 步存重启文件 |
| 运行 | `make test_earth` |

### 4. 木星磁层 (`Param/JUPITER/`)

**物理场景**：快速旋转行星磁层（10h 周期）+ 环面等离子体

| 参数 | 值 |
|------|-----|
| 根块 | 6×4×4 |
| 网格范围 | x:[-640,128] y:[-256,256] z:[-256,256] |
| AMR | 初始区域(3.0) + 球体(1.5 at r=4) + 环带(3/8 at h=8) |
| 旋转 | T（周期=10h） |
| rBody | 3.0 |
| 太阳风 | n=0.2/cc, T=200000K, Ux=-400 km/s, Bz=-1nT |
| 最大块数 | 1440 |

### 5. 金星 (`Param/VENUS/`)

**物理场景**：无内禀磁场的行星-太阳风相互作用（离子拾取）

| 参数 | 值 |
|------|-----|
| 归一化 | SOLARWIND |
| 网格层 | nLevel=2 + initial shape |
| CFL | 0.2 |
| 离子种类 | H⁺, CO₂⁺, O₂⁺, O⁺ |
| rBody | 1.023 |
| 模式 | 测试（`read_inputs amr`） |
| 边界 | 盒子: outflow(前)/inflow(后)，内边界: user |

### 6. 火星 (`Param/MARS/`)

**物理场景**：多流体火星电离层-太阳风作用

| 参数 | 值 |
|------|-----|
| 离子种类 | H⁺, O₂⁺, O⁺, CO₂⁺ |
| 通量格式 | Linde (HLL) |
| 限制器 | minmod |
| CFL | 0.8 |
| rBody | 1.029 |
| 太阳风 | n=4/cc, T=3.5e5K, Ux=-500 km/s, Parker 螺旋 B |
| 重力 | 开启（X 方向） |
| 内边界 | user（火星用户模块覆盖密度） |

### 7. 土卫六泰坦 (`Param/TITAN/`)

**物理场景**：多离子（6 种）大气逃逸模拟

| 参数 | 值 |
|------|-----|
| 离子种类 | L⁺, M⁺, H₁⁺, H₂⁺, MHC⁺, HHC⁺, HNI⁺（7 种） |
| 网格层 | nLevel=3 |
| 最大块数 | 700 |
| CFL | 0.8 |
| 通量格式 | Linde（HLL） |
| 限制器 | minmod |
| 测试 | `init_hall_resist`（霍尔电阻率测试） |
| 太阳风 | n=2.9/cc, T=3.11e5K, Ux=120 km/s |

---

## 三、彗星案例

### 8. 哈雷彗星 (`Param/COMET/`)

**物理场景**：彗星-太阳风相互作用，光致电离产生拾取离子

| 参数 | 值 |
|------|-----|
| 彗星种类 | H₂O⁺, H⁺, H₃O⁺, OH⁺, O⁺, CO⁺（6 种） |
| 产率 Qprod | 7×10²⁹ s⁻¹ |
| 中性速度 | 1 km/s |
| 网格 | 1×1×1 根块 [-0.8,0.2]×[-0.5,0.5]³ |
| 太阳风 | n=8/cc, T=1.8e5K, Ux=400 km/s, B=(3.4,3.4,0)nT |
| 通量格式 | Linde（HLL） |
| 限制器 | minmod |
| CFL | 0.8 |
| 最大块数 | 600 |

---

## 四、日冕与日球层

### 9. 日冕 AWSoM (`Param/CORONA/`)

**物理场景**：Alfvén Wave Solar Model——从色球层到日冕的太阳风加速

| 可用 PARAM 文件 | 说明 |
|----------------|------|
| `PARAM.in.Awsom` | 标准 AWSoM |
| `PARAM.in.Awsom.GPU` | GPU 加速版 |
| `PARAM.in.Awsom.large.GPU` | 大网格 GPU 版 |
| `PARAM.in.Awsom.signb` | 含磁场符号 |
| `PARAM.in.AwsomChargeState` | 含电荷态 |
| `PARAM.in.AwsomR` | 球坐标版 |
| `PARAM.in.1Dwedge` / `2Dwedge` | 1D/2D 楔形域 |
| `PARAM.in.bvector` / `awsom.bvector` | 含磁场矢量 |

> 注：此目录还包含磁图数据（HMI/MDI）、辐射冷却表等资源文件。

### 10. 日球层 (`Param/HELIOSPHERE/`)

**物理场景**：大尺度外日球层结构

| 子目录 | 说明 |
|--------|------|
| `Arcade/` | 磁拱模拟 |
| `CME_fluxrope/` | CME 磁通量绳 |
| `Heliosphere/` | 全球日球层 |
| 附加文件 | `EARTH_TRAJ.in`（地球轨道）、`SATELLITE.in` |

### 11. GEM 磁重联挑战 (`Param/GEMRECONNECTION/`)

**物理场景**：标准 Geospace Environment Modeling 磁重联基准

| 可用 PARAM 文件 | 说明 |
|----------------|------|
| `PARAM.in.MhdHypPe` | 双曲型 MHD + 电子压强 |
| `PARAM.in.double.CS.MHDEPIC` | 双电流片 MHD-EPIC 耦合 |
| `PARAM.in.sixmoment` | 六矩模型 |

---

## 五、其他案例速览

| 目录 | 物理场景 | 要点 |
|------|---------|------|
| AMR | 自适应网格测试 | 纯网格功能验证 |
| ANISOPRESSURE | 各向异性压强 | 离子回旋各向异性 |
| B0 | 背景磁场 | 偶极场初始条件 |
| COMET3FLUIDSPE | 彗星3流体+电子压强 | 多流体彗星 |
| CRASH | 辐射激波 | 辐射流体力学 |
| CURRENT | 电流片 | 电流片不稳定性 |
| EUROPA | 木卫二 | 2/3 流体模型 |
| FIVEMOMENT | 五矩模型 | 高阶矩闭合 |
| GANYMEDE | 木卫三 | 固有磁层 |
| MERCURY | 水星 | 弱内禀磁场 |
| MOONIMPACT | 月球撞击 | 撞击产生瞬态大气 |
| MULTIFLUID/MULTIION | 多流体/多离子 | 通用框架 |
| OUTERHELIO | 外日球层 | 中性流体+拾取离子 |
| SATURN | 土星 | 环面+快速旋转 |
| SIXMOMENT | 六矩模型 | 最高阶矩闭合 |
| VISCOSITY | 粘性测试 | 物理粘性项 |
| MAGNETOMETER | 虚拟磁强计 | 卫星数据同化 |
| ROTATINGFRAME | 旋转坐标系 | 共转系统 |
| FLUXEMERGENCE | 磁通浮现 | 太阳表面磁通 |
| ROSETTA | Rosetta 任务 | 彗星 67P |
| REGION | 区域分解 | 多区域耦合 |
| RUNGEKUTTA | Runge-Kutta | 时间格式测试 |
| MHDNONCONS | 非守恒MHD | 非守恒公式 |
| CYLOUTFLOW | 柱状出流 | 恒星风 |
| 2BODYPLOT | 二体绘图 | 可视化辅助 |

---

## 六、如何运行

```bash
# 编译后运行单个测试
make -j test_shocktube

# 指定并行度
make -j test_earth NP=4

# 串行+OpenMP
make -j test_mars MPIRUN= NTHREAD=4

# GPU 测试
make -j test_small_gpu NP=1
```

所有测试产生 `.diff` 文件——空（0 字节）表示通过。
