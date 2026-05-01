# lingtai-batsrus-skill

> Progressive disclosure skill for **BATSRUS** — the Block-Adaptive Tree Solar-wind Roe-type Upwind Scheme.
>
> University of Michigan CSEM | SWMF core MHD engine

[![skill](https://img.shields.io/badge/lingtai-skill-blue)](https://github.com/orgs/Lingtai-AI)
[![Fortran](https://img.shields.io/badge/analyses-Fortran%2090-purple)](https://github.com/SWMFsoftware/BATSRUS)
[![Lines](https://img.shields.io/badge/docs-3021%20lines-green)](./SKILL.md)

---

## What is this?

A **progressive disclosure skill** for the LingTai agent network — a structured, layered guide to the BATSRUS MHD solver codebase. Instead of dumping everything at once, this skill unfolds knowledge in four layers:

| Layer | What it covers | When to load |
|-------|---------------|--------------|
| **0 — Identity** | One paragraph: what BATSRUS is | Always safe to share |
| **1 — Map** | 100+ module taxonomy, 75 equation families, file tree | "What's inside?" |
| **2 — Domain** | Deep-dives: solver / grid / parallel / equations | Drill into one area |
| **3 — Source** | File paths + line numbers + code snippets | Implementation details |

Plus **Layer 2.5** — case studies with `PARAM.in` analysis, benchmarks, nightly test suite, and key papers.

## Repository structure

```
.
├── SKILL.md                          # Main entry point (4-layer progressive disclosure)
└── reference/
    ├── batsrus_solver_report.md      # Roe 7-wave solver + MUSCL + GMRES
    ├── batsrus_grid_report.md        # Octree AMR + Morton ordering + 15 criteria
    ├── batsrus_parallel_report.md    # MPI + OpenACC GPU + OpenMP + scalability
    ├── batsrus_equation_report.md    # 75 equation modules fully catalogued
    ├── batsrus_master_report.md      # Consolidated report + architecture diagram
    ├── case_studies.md               # 10 PARAM.in cases with key parameters
    ├── benchmarks.md                 # MHD verification + GPU performance data
    ├── nightly_tests.md              # 50+ tests from Makefile.test
    └── toth_resources.md             # Gábor Tóth's key papers & BATL library
```

## What is BATSRUS?

**BATSRUS** (Block Adaptive Tree Solar-wind Roe-type Upwind Scheme) is a 3D magnetohydrodynamics solver developed at the University of Michigan. It uses a finite-volume method with a Roe approximate Riemann solver on a block-adaptive octree grid. It is the core engine of the Space Weather Modeling Framework (SWMF).

- **Language**: Fortran 90 (~200,000+ lines)
- **Grid**: Block-adaptive octree (2:1 refinement, Morton ordering)
- **Numerics**: FVM + 7-wave Roe solver + Harten entropy fix + MUSCL/TVD/CWENO5
- **Time**: Explicit multi-stage RK + semi-implicit GMRES+BILU + full-implicit BDF2
- **Parallel**: MPI domain decomposition + OpenMP threading + OpenACC GPU
- **Physics**: 75 equation modules — ideal MHD, Hall MHD, multi-ion, multi-fluid, radiation, AWSoM corona, outer heliosphere, comets, planetary magnetospheres
- **Performance**: Single A100 ≈ 270 AMD Rome cores, 3.6× faster-than-real-time for magnetospheric simulations

## Quick highlights

- 🔬 **Roe 7-wave Riemann solver** — the default flux scheme with Harten entropy fix
- 🌳 **Block-adaptive octree** — 15+ physics-based refinement criteria
- ⚡ **OpenACC GPU** — ~1% code modification for full GPU acceleration
- 🧪 **50+ nightly tests** — Brio-Wu shock tube, Orszag-Tang, GEM reconnection, all planets
- ☀️ **AWSoM solar model** — Alfvén wave-driven corona and solar wind
- 🪐 **Planetary coverage** — Earth, Mars, Venus, Titan, Jupiter, Saturn, Mercury, Europa, Ganymede, comets

## How to use (LingTai agents)

```bash
# The skill is discoverable in your library catalog
/library        # Lists all available skills

# Load on demand
read .library/custom/batsrus/SKILL.md          # Layer 0-1
read .library/.../reference/benchmarks.md       # Layer 2.5
read .library/.../reference/batsrus_solver_report.md  # Layer 3
```

## How this was built

This skill was created by a LingTai orchestrator agent (`deepseek_pro`) that **spawned 4 sub-agents (avatars)** to analyze BATSRUS in parallel:

- `batsrus_solver` — core numerics (Roe solver, MUSCL, time stepping)
- `batsrus_grid` — AMR grid system (octree, Morton, load balance)
- `batsrus_parallel` — parallel architecture (MPI, OpenACC, OpenMP)
- `batsrus_equation` — physics equation modules (75 modules catalogued)

Each avatar pulled source code directly from GitHub, analyzed it, and produced a structured report. The orchestrator then consolidated everything into this progressive disclosure skill.

## References

- BATSRUS repository: https://github.com/SWMFsoftware/BATSRUS
- BATL library: https://github.com/SWMFsoftware/BATL
- Gábor Tóth's papers: https://public.websites.umich.edu/~gtoth/Papers/
- GPU paper: [arXiv:2501.06717](https://arxiv.org/abs/2501.06717)

## License

This skill is documentation of the BATSRUS open-source codebase. BATSRUS itself is developed by the University of Michigan CSEM group. See [LICENSE.txt](https://github.com/SWMFsoftware/BATSRUS/blob/master/LICENSE.txt) in the BATSRUS repository.
