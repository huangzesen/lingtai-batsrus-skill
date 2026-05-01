---
name: batsrus
description: Progressive disclosure guide to BATSRUS — the Block-Adaptive Tree Solar-wind Roe-type Upwind Scheme (MHD solver, SWMF core). Four layers: overview → module maps → deep-dive references → source code. Load this skill when you need to understand, use, or extend BATSRUS. Never dump everything at once — unfold layers only as the task demands.
version: 1.0.0
---

# BATSRUS — Progressive Disclosure

> BATS-R-US: **B**lock **A**daptive **T**ree **S**olar-wind **R**oe-type **U**pwind **S**cheme
>
> University of Michigan CSEM | SWMF core MHD engine
>
> GitHub: https://github.com/SWMFsoftware/BATSRUS

---

## How This Skill Works

This is a **layered, progressive disclosure** skill. Do **not** dump all layers at once. Start at Layer 0 and only descend to the next layer when the current context demands it.

| Layer | Content | When to load |
|-------|---------|---------------|
| **0 — Identity** | One-paragraph "what is BATSRUS" | Always safe to share |
| **1 — Map** | Module taxonomy, file tree, flow overview | User asks "what's inside?" or "how is it organized?" |
| **2 — Domain** | Specific domain deep-dive (solver / grid / parallel / equations) | User needs detail on one area |
| **2.5 — Case Studies** | Running examples with PARAM.in analysis (36 cases) | User asks "how do I run this?" or "show me an example" |
| **3 — Source** | File paths, line numbers, code snippets, benchmarks | User needs implementation-level detail or performance data |

The `reference/` directory contains the full deep-dive reports, case studies, and benchmarks from which all layers are drawn. Load them directly for Layer 2.5/3 needs.

---

## Layer 0: Identity

**BATSRUS** is a 3D magnetohydrodynamics (MHD) solver developed at the University of Michigan. It uses a finite-volume method with a Roe approximate Riemann solver on a block-adaptive tree (octree) grid. It is the core engine of the Space Weather Modeling Framework (SWMF), supporting simulations from ideal MHD to multi-fluid, multi-ion, radiative, and outer-heliosphere models across planetary magnetospheres, the solar corona, solar wind, and comets.

- **Language**: Fortran 90 (~200,000+ lines)
- **Grid**: Block-adaptive octree (self-similar rectangular blocks, 2:1 refinement)
- **Numerics**: FVM + Roe solver (7-wave + Harten entropy fix) + MUSCL/TVD reconstruction
- **Time**: Explicit multi-stage RK (default Dt/2-Dt second-order); optional semi-implicit (GMRES+BILU) and full-implicit (BDF2)
- **Parallel**: MPI domain decomposition + OpenMP threading + OpenACC GPU (nvfortran)
- **Physics**: 75 equation modules — ideal MHD, Hall MHD, multi-ion, multi-fluid, radiation, AWSoM coronal model, outer heliosphere with pickup ions, comets, planets (Mars, Venus, Titan, Jupiter, Saturn, Europa, Ganymede, Mercury)
- **Coupling**: 9 SWMF interfaces (ionosphere IE, inner magnetosphere IM, PIC, radiation belt RB, etc.)

---

## Layer 1: Module Map

### Top-Level Directory Structure

```
BATSRUS/
├── src/              ← Core solver (100+ .f90 modules)
├── srcEquation/      ← Physics equation modules (75)
├── srcUser/          ← User/application modules (33)
├── srcInterface/     ← SWMF coupling interfaces (14 files)
├── srcPostProc/      ← Post-processing (spectra, interpolation)
├── Param/            ← Parameter files (36+ applications)
├── Doc/              ← LaTeX documentation (USERMANUAL.pdf)
├── Scripts/          ← Perl, IDL, Tecplot utilities
├── PARAM.XML         ← 475KB XML — all input parameters
├── Config.pl         ← Build configuration
└── Makefile          ← Top-level build
```

### Core Solver Modules (`src/`)

| Category | Key Modules |
|----------|-------------|
| **Driver** | `ModMain.f90`, `ModBatsrusMethods.f90` |
| **Time Integration** | `ModAdvance.f90`, `ModAdvanceExplicit.f90` |
| **Riemann Solver** | `ModFaceFlux.f90` (4309 lines!), `ModPhysicalFlux.f90`, `ModCharacteristic.f90` |
| **Reconstruction** | `ModFaceValue.f90` (MUSCL/TVD/WENO), `ModFaceGradient.f90` |
| **AMR Grid** | `ModAMR.f90`, `ModGeometry.f90`, `ModBlockData.f90`, `ModNodes.f90` |
| **Parallel** | `ModParallel.f90`, `ModMessagePass.f90`, `ModLoadBalance.f90` |
| **GPU** | `ModBuffer.f90`, `ModBatlInterface.f90` (OpenACC directives) |
| **Implicit** | `ModSemiImplicit.f90` (GMRES+BILU), `ModPointImplicit.f90` |
| **divB** | `ModCleanDivB.f90`, `ModConstrainDivB.f90`, `ModProjectDivB.f90` |
| **Physics Sources** | `ModResistivity.f90`, `ModHeatConduction.f90`, `ModRadiativeCooling.f90`, `ModHallResist.f90`, `ModIonElectron.f90`, `ModLaserHeating.f90` |
| **I/O** | `ModIO.f90`, `ModHdf5.f90`, `ModRestartFile.f90`, `ModWritePlot.f90` |

### Physics Equation Families (`srcEquation/`)

75 modules, nVar ranging from 1 to 40:

| Family | Count | nVar Range | Key Modules |
|--------|-------|------------|-------------|
| **Ideal MHD** | ~30 | 8–17 | `ModEquationMhd`, `ModEquationMhdPe`, `ModEquationMhdAnisoP` |
| **AWSoM (Solar Corona)** | ~10 | 10–12 | `ModEquationAwsom`, `ModEquationAwsomSA`, `ModEquationAwsomFluids` |
| **Outer Heliosphere** | ~8 | 28–40 | `ModEquationOuterHelio`, `ModEquationOuterHelioPUIPe`, `ModEquationOuterHelioAwsom` |
| **HD (Euler)** | ~5 | 5–6 | `ModEquationHd`, `ModEquationHdEos`, `ModEquationHdEosRad` |
| **Multi-Moment** | ~4 | 5–6 moments | `ModEquationFiveMoment`, `ModEquationSixMoment` |
| **Multi-Ion** | ~4 | 13–15 | `ModEquationMultiIon`, `ModEquationSwh` |
| **CRASH (Radiative)** | ~5 | 7–11 | `ModEquationCrash`, `ModEquationMhdCrash` |
| **Planetary** | ~12 | 8–28 | `ModEquationMhdMars`, `ModEquationMhdTitan`, `ModEquationMhdComet` |
| **Solar Wind** | ~3 | 9–20 | `ModEquationSwh`, `ModEquationSwhPui` |
| **Scalar** | 1 | 1 | `ModEquationScalar` |

### User Applications (`srcUser/`)

33 modules covering: Mars, Venus, Titan, Jupiter, Saturn, Mercury, Europa, Ganymede, Comet (multiple variants), Outer Heliosphere, Flux Rope, Kelvin-Helmholtz, Gem Reconnection, Solar Corona, Waves, etc.

### SWMF Coupling (`srcInterface/`)

9 coupling interfaces: IE (ionosphere electrodynamics), IM (inner magnetosphere), PC (PIC), PS (plasmasphere), PT (PIC transport), PW (polar wind), RB (radiation belt), UA (upper atmosphere), plus wrapper and grid descriptor.

---

## Layer 2: Domain Deep-Dives

Ask the user which domain they need, then load the corresponding reference:

- **Solver**: `reference/batsrus_solver_report.md` — time stepping, Roe 7-wave solver, MUSCL reconstruction, CFL, divB schemes, semi-implicit
- **Grid**: `reference/batsrus_grid_report.md` — octree structure, AMR criteria, Morton ordering, load balance, Berger-Colella flux correction, coordinate systems
- **Parallel**: `reference/batsrus_parallel_report.md` — MPI domain decomposition, ghost cell exchange, load balancing, OpenACC GPU, OpenMP threading, scalability data
- **Equations**: `reference/batsrus_equation_report.md` — all 75 modules enumerated, equation families, nVar/nWave tables, SWMF coupling

Each reference report is a self-contained ~500-line Markdown with code citations and line numbers.

---

## Layer 2.5: Case Studies & Benchmarks

When the user wants to see running examples or performance data:

- **Case Studies**: `reference/case_studies.md` — 36 `Param/` directories analyzed. Key examples:
  - **Brio-Wu Shock Tube** — standard MHD verification, `make test_shocktube`
  - **Earth Magnetosphere** — realistic solar wind driving, multi-level AMR
  - **Mars / Venus / Titan** — multi-fluid planetary ionospheres
  - **Jupiter / Saturn** — rotating magnetospheres
  - **Comet Halley** — 6-species coma + pickup ions
  - **Solar Corona (AWSoM)** — GPU-ready Alfvén wave solar model
  - **GEM Reconnection** — standard reconnection challenge (MHD-EPIC)
- **Benchmarks**: `reference/benchmarks.md` — standard MHD verification suite (Brio-Wu, Orszag-Tang, Rotor, GEM), GPU performance data (arXiv:2501.06717), historical scalability, validation methodology
- **How to run**: `make test_<case>` / `make -j test NP=4` / `make -j test_gpu`

## Layer 3: Source-Level

All reference files cite exact file paths from the GitHub repository. To read actual source:

```bash
# Fetch any file from the repo
curl -s https://raw.githubusercontent.com/SWMFsoftware/BATSRUS/master/<path>

# Example: the Roe solver characteristic decomposition
curl -s https://raw.githubusercontent.com/SWMFsoftware/BATSRUS/master/src/ModCharacteristic.f90
```

Key file → purpose mapping is in the module overview above (Layer 1).

---

## Quick Reference: Technical Highlights

1. **7-wave Roe solver** with Harten entropy fix — the default Riemann solver (`src/ModCharacteristic.f90` for wave decomposition, `src/ModFaceFlux.f90` for `roe_solver_new`)
2. **Block-adaptive octree** with 2:1 refinement ratio, 15+ physics-based AMR criteria (gradT, gradP, ∇·u, J, current sheet, B·∇, etc.)
3. **Berger-Colella flux correction** across grid resolution changes for conservation (`src/ModConserveFlux.f90`)
4. **Morton (Z-order) space-filling curve** for load balancing across MPI ranks (`src/ModLoadBalance.f90`)
5. **OpenACC GPU** — ~1% of solver code annotated with `!$acc` directives; single A100 ≈ 270 AMD Rome cores; 3.6× faster-than-real-time for magnetospheric simulations
6. **Semi-implicit GMRES+BILU** for stiff source terms (radiation, heat conduction, Hall effect)
7. **8-wave Powell scheme** for divB control — source term proportional to ∇·B
8. **MUSCL + minmod** default reconstruction; optional CWENO5, logarithmic limiter, and various slope limiters
9. **AWSoM solar corona model** uses nWave=2 Alfvén wave decomposition instead of full Riemann solver — astrophysics-specific optimization
10. **75 physics equation modules** span from 1-variable scalar to 40-variable outer heliosphere with 2 PUIs + 4 neutral fluids + electrons

---

## Master Summary

For the full consolidated report covering all four domains: `reference/batsrus_master_report.md` (407 lines, ~20KB). This is suitable as a comprehensive reference document.
