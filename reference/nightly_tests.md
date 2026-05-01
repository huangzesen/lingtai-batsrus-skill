# BATSRUS Nightly Test Suite（完整测试列表）

> 来源：`Makefile.test`（`make test_help` 输出）
> 运行：`make -j test`（串行运行所有测试，用 `NP=4` 指定并行度）
> 机制：每个 `make test_xxx` = `compile → run → diff` 三步，`.diff` 空文件 = 通过

---

## 完整测试列表（50+ 项）

### GPU 测试（OpenACC）

| 测试 | 说明 |
|------|------|
| `test_earth` / `test_earth_large_gpu` | 地球磁层 GPU 测试 |
| `test_fastwave` / `test_fastwave_2d` | 3D/2D 快波传播 |
| `test_shocktube_1d` | 1D 激波管 GPU |
| `test_awsom_gpu` / `test_awsom_large_gpu` | AWSoM 日冕 GPU |
| `test_awsom_bvector` | AWSoM + 矢量B场 |
| `test_bvector` | B矢量日冕 |

### 标准 MHD 验证

| 测试 | 说明 |
|------|------|
| `test_shocktube` | Brio-Wu 激波管 |
| `test_shockround` | 圆管激波管 |
| `test_shockramp` | 激波斜台 |
| `test_hallmhd` | 霍尔 MHD |
| `test_twofluidmhd` | 双流体 MHD |
| `test_multifluid` | 多流体 |
| `test_multiion` | 多离子 |
| `test_mhdnoncons` | 非守恒 MHD |
| `test_viscosity` | 粘性测试 |
| `test_partsteady` | 部分定常 |
| `test_fluxemergence` | 磁通浮现 |

### AMR 测试

| 测试 | 说明 |
|------|------|
| `test_amr` | 笛卡尔 AMR |
| `test_amrsph` | 球坐标 AMR |
| `test_timewarp_1d` / `test_timewarp_2d` | 时间翘曲 AMR |

### 高阶矩模型

| 测试 | 说明 |
|------|------|
| `test_fivemoment_alfven` | 五矩 Alfvén 波 |
| `test_fivemoment_shock` | 五矩激波 |
| `test_fivemoment_langmuir` | 五矩 Langmuir 波 |

### 日冕/太阳风

| 测试 | 说明 |
|------|------|
| `test_awsom` | 标准 AWSoM |
| `test_awsomr` | AWSoM-R（球坐标） |
| `test_stitch` | STITCH 拼接 |
| `test_earthsph` | 球坐标地球 |
| `test_mercurysph` | 球坐标水星 |
| `test_bx0` | 3D Bx0 等值面 |

### 外日球层

| 测试 | 说明 |
|------|------|
| `test_outerhelio` | 3D 外日球层 |
| `test_outerhelio2d` | 2D 外日球层 |
| `test_outerhelio_1d` | 1D 外日球层 |
| `test_outerhelioawsom` | 外日球层 + 湍流 |
| `test_outerheliopui` | 外日球层 + PUI |
| `test_outerheliopui_1d` | 1D 外日球层 + PUI |

### 行星

| 测试 | 说明 |
|------|------|
| `test_mars` / `test_mars_start` / `test_mars_restart` | 火星全链路 |
| `test_marsfluids` | 火星多流体 |
| `test_ccmc_mars` | 火星 CCMC 版 |
| `test_venus` | 金星 |
| `test_titan` / `test_titan_start` / `test_titan_restart` | 土卫六全链路 |
| `test_jupiter` | 木星 |
| `test_saturn` | 土星 |

### 彗星

| 测试 | 说明 |
|------|------|
| `test_comet` | 彗星 |
| `test_cometCGhd` | 彗星 CG HD |
| `test_cometCGfluids` | 彗星 CG 流体 |

### 其他

| 测试 | 说明 |
|------|------|
| `test_anisotropic` | 各向异性压强 |
| `test_spectrum` | 光谱后处理 |
| `test_magnetometer` | 虚拟磁强计 |
| `test_moonimpact` | 月球撞击 |
| `test_2bodyplot` | 二体绘图 |
| `test_region2d` | 2D BATL 区域 |
| `test_L1toBC` | L1 到边界条件传播 |

### Extra 测试（需 `srcUserExtra/` 存在）

| 测试 | 说明 |
|------|------|
| `test_awsomchargestate` | 日冕 + 电荷态 |
| `test_eosgodunov` | EOS/Godunov |
| `test_europa2fluidspe` | 木卫二 2 流体 + Pe |
| `test_europa3fluidspe` | 木卫二 3 流体 + Pe |
| `test_ganymede` | 木卫三 |
| `test_graydiffusion` | 灰体扩散 |
| `test_laserpackage` | CRASH 激光 |

---

## 运行机制（三步流水线）

```makefile
test_xxx: test_xxx_compile → test_xxx_rundir → test_xxx_run → test_xxx_check
```

1. **编译**——生成 `BATSRUS.exe`
2. **建运行目录**——复制 PARAM.in + 数据文件到 `run_test_xxx/`
3. **运行**——`mpiexec -n NP BATSRUS.exe`
4. **检查**——`DiffNum.pl` 对比输出与参考解 → `.diff` 文件
5. **Bless**：`BLESS=YES` 将当前输出升级为新的参考解

## GPU 测试子集

```bash
# 完整 GPU 套件
make -j test_gpu

# 小规模预检（单 GPU）
make -j test_small_gpu NP=1

# 默认：2 MPI + 2 OpenMP 线程
make -j test
```
