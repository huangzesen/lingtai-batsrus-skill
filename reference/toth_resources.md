# Gábor Tóth — BATSRUS 相关关键资源

> 来源：https://public.websites.umich.edu/~gtoth/
> （主页受 Cloudflare 保护，Papers/ 子目录可直访）

---

## 关键论文（均可在 Papers/ 目录获取）

| 论文 | 年份 | 主题 | DOI/链接 |
|------|------|------|---------|
| **Adaptive Numerical Algorithms in Space Weather Modeling** | 2012 | BATSRUS + BATL 综述 | [Toth2012review.pdf](https://public.websites.umich.edu/~gtoth/Papers/Toth2012review.pdf) |
| **The ∇·B=0 Constraint in Shock-Capturing MHD Codes** | 2000 | divB 处理的经典文献 | [Toth2000_divb.pdf](https://public.websites.umich.edu/~gtoth/Papers/Toth2000_divb.pdf) |
| **AMR with Conservative Transport** | 2002 | AMR 通量修正（Berger-Colella） | JCP 180, 736-750 |
| **Parallel Hall MHD Scheme** | 2008 | 霍尔 MHD + 隐式求解器 | [Toth2008_hall.pdf](https://public.websites.umich.edu/~gtoth/Papers/Toth2008_hall.pdf) |
| **BATSRUS GPU** | 2025 | OpenACC GPU 加速，3.6× real-time | [arXiv:2501.06717](https://arxiv.org/abs/2501.06717) |

## 软件

| 项目 | 说明 | 仓库 |
|------|------|------|
| **BATL** | Block Adaptive Tree Library——BATSRUS 底层网格库 | [github.com/SWMFsoftware/BATL](https://github.com/SWMFsoftware/BATL) |
| **BATSRUS** | 主求解器 | [github.com/SWMFsoftware/BATSRUS](https://github.com/SWMFsoftware/BATSRUS) |
| **READAMR** | 读取/插值 BATSRUS 输出的辅助库 | 含在 BATL 中 |

## BATL 库细节

- **作者**：Gábor Tóth, Bart van der Holst, Lars Daldorff
- **语言**：Fortran 90 + C
- **许可证**：仅限密歇根大学开发的代码使用（BATSRUS、READAMR 等）
- **核心功能**：
  - AMR（八叉树）管理：`iTree_IA` 数据结构
  - 消息传递：`message_pass_cell`（Ghost cell 交换）
  - 负载均衡：基于 Morton 空间填充曲线
  - I/O：HDF5 并行输出

## 附加演示/资料

- Blue Waters 2017 研讨会：BATSRUS 空间天气模拟（含日冕/太阳风/磁层耦合）——[幻灯片](https://bluewaters.ncsa.illinois.edu/liferay-content/document-library/2017+Symposium/toth-space.pdf)
- GEM 重联挑战：BATSRUS + MHD-EPIC 耦合

## 如何从 Tóth 的网站获取更多

Papers 目录结构（受 Cloudflare 保护，可通过 GitHub→SWMF 链绕过）：
```
https://public.websites.umich.edu/~gtoth/Papers/
├── Toth2000_divb.pdf
├── Toth2002_amrct.pdf
├── Toth2008_hall.pdf
├── Toth2012review.pdf
└── ... (更多论文)
```

大多数 BATSRUS 相关论文已被镜像到 arXiv 或可在 ResearchGate 获取。
