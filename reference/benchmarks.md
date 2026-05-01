# BATSRUS Benchmark & Validation 汇总

> 验证测试套件、性能基准与已知验证状态

---

## 一、标准 MHD 验证测试

### 1. Brio-Wu 激波管

| 属性 | 值 |
|------|-----|
| 类型 | 1D MHD 黎曼问题 |
| 物理 | 复合激波 + 稀疏波 + 接触间断 + 慢激波 |
| 参数 | nStage=2, CFL=0.8, Roe+mc3, 256×4×4 |
| 验证方式 | 与参考解（Brio & Wu, JCP 1988）比较密度/磁场剖面 |
| BATSRUS 测试 | `make test_shocktube` |
| 状态 | ✅ 标准验证 |

### 2. Orszag-Tang 涡旋

| 属性 | 值 |
|------|-----|
| 类型 | 2D 不可压缩→可压缩 MHD 湍流 |
| 物理 | 激波形成 + 激波-激波相互作用 |
| 参数 | 周期边界, γ=5/3 |
| 验证方式 | 密度/压强等值线 + 动能/磁能时间演化 |
| 参考文献 | Orszag & Tang, JFM 1979 |
| 状态 | ✅ 社区标准基准 |

### 3. MHD 转子 (Rotor)

| 属性 | 值 |
|------|-----|
| 类型 | 2D MHD |
| 物理 | 刚性旋转圆盘嵌入静止磁化等离子体 |
| 验证方式 | Mach 数分布 + 磁力线扭曲 |
| 参考文献 | Balsara & Spicer, JCP 1999 |
| 状态 | ✅ 社区标准基准 |

### 4. GEM 磁重联挑战

| 属性 | 值 |
|------|-----|
| 类型 | 2D 无碰撞磁重联 |
| 物理 | Harris 电流片 + 初始扰动 |
| BATSRUS 实现 | MhdHypPe / double.CS.MHDEPIC / sixmoment |
| 验证方式 | 重联率 vs 时间（Birn et al. 2001 基准） |
| 参考文献 | Birn et al., JGR 2001 |
| 状态 | ✅ 有 MHD-EPIC 耦合版 |

### 5. Sod 激波管

| 属性 | 值 |
|------|-----|
| 类型 | 1D 纯流体力学 |
| 验证方式 | 密度/速度/压强剖面 |
| 状态 | ✅ HD 模块标准测试 |

---

## 二、性能 Benchmark

### GPU 加速性能（arXiv:2501.06717）

| 指标 | 数据 |
|------|------|
| GPU 方案 | OpenACC（nvfortran 编译） |
| 代码修改量 | ~1% |
| 单 A100 等效 | ≈ 270 AMD Rome 核心 |
| 磁层模拟 | 3.6× faster-than-real-time |
| 单节点 4 A100 | 95% 强扩展效率 |
| 256 A100 弱扩展 | 50-60% 效率 |

**测试案例**：地球磁层（fastwave + earth）

### 历史扩展性数据

| 时期 | 平台 | 处理器数 | 性能 |
|------|------|---------|------|
| 1999 | Cray T3E | 1,490 PEs | 342 GFlops |
| 2005 | Columbia (SGI Altix) | 1,024 cores | — |
| 2012 | Pleiades | 10,000+ cores | — |
| 2024 | Frontera/Pleiades GPU | 256 A100 | 3.6× real-time |

---

## 三、物理验证体系

### 太阳风/磁层

| 验证项 | 参考数据 |
|--------|---------|
| 弓激波位置 | 卫星统计（Farris & Russell 1994） |
| 磁层顶位置 | Shue et al. 1998 模型 |
| 极盖电位 | Weimer 2005 / AMIE 同化 |
| 等离子体片密度 | Geotail 统计 |
| Dst 指数 | 地磁暴事件 |

### 日冕/太阳风

| 验证项 | 参考数据 |
|--------|---------|
| 冕洞形态 | SDO/AIA + STEREO EUV 成像 |
| 太阳风速 | Ulysses/ACE 原位测量 |
| 磁场极性 | PFSS 外推 |
| 日冕温度 | 光谱诊断 |

### 行星

| 行星 | 验证数据 |
|------|---------|
| 火星 | MAVEN/MEX 原位测量 |
| 金星 | VEX 磁场+等离子体 |
| 木星 | Galileo/Juno 磁场模型 |
| 土星 | Cassini 原位 |
| 土卫六 | Cassini TA-TB 飞掠 |
| 彗星 | Giotto/Rosetta |

---

## 四、代码验证方法

### 测试套件

```bash
# 完整测试（2 MPI + 2 OpenMP）
make -j test

# GPU 测试
make -j test_gpu

# 小规模 GPU 预检
make -j test_small_gpu NP=1

# 列表所有测试
make test_help
```

### 验证流程

1. 编译 + 运行 → 产生输出文件
2. 与参考输出比较 → `.diff` 文件
3. 空 `.diff` = 通过
4. 参考文件含在仓库中（`.ref.gz`）

---

## 五、已知局限

| 问题 | 影响 | 缓解措施 |
|------|------|---------|
| divB 增长 | 长时磁重联 | 8波Powell源项 + 投影法 |
| 低速流误差 | 太阳风亚声速区 | 低马赫数修正 |
| 负压强 | 强激波/磁重联 | 压强限幅器 + Boris 修正 |
| GPU 浮点精度 | 与 CPU 微小差异 | 双精度（默认） |
