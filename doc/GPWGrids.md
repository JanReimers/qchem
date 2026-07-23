# GPW grid inventory — every grid, what sizes it, knob vs algorithm

Written 2026-07-21 in answer to the user story in `doc/GPWPlan.md` §0e ("table of all grid usages, how
grid range/spacing is decided, and if it is a user knob or decided by a sensible algo").  Verified
against the code as of branch `gpw-0e-pp-local-split`; the run-start `[GPW grid]` diagnostic
(`GPW_Evaluator::ReportGrids`) prints every row of this table live, so the doc and a log can always be
cross-checked.  CP2K comparison numbers are the NaF deck (`doc/CP2Kresults.md`).

## 0. "Ecut for *what*?" — the disambiguation

Every `Ecut` in the periodic code is one of exactly four things.  A G-space cutoff always means the
BALL \(\{G : \tfrac12|G|^2 < E_{cut}\}\); the raster N is derived from a ball, never chosen directly.

| name in code | what it bounds | who has one |
|---|---|---|
| `Ecut` (PW) | the ORBITAL plane-wave basis ball | `PlaneWave_IBS` only.  **GPW has NO orbital cutoff** (Gaussians are analytic). |
| `densityEcut` | the DENSITY/collocation reference ball (row 1 below) | GPW's ONLY grid cutoff |
| ladder `ecut_L` | per-level sub-balls derived from `densityEcut` (row 4) | derived, never user-set |
| `relCutoff` | a MULTIPLIER on the ρ ball for the Vxc grid (row 3) | functional-derived; ==1 for LDA |

## 1. The inventory

**Row 1 — the FFT ρ↔G reference grid** (`GPW_Evaluator::itsFFT_R_G_Grids`, a `PW_Grid_Evaluator`).
- **Ball:** `densityEcut`, a THREE-MODE knob: `<0` = **automatic** (recommended): floor =
  `cutoffFactor·α_max` from the basis (algorithm); `=0` = DFT tier off; `>0` = explicit Hartree value,
  honoured with a `cerr` warning if below the floor.  `cutoffFactor` (default 8) is itself a calibrated
  constant, not a per-run knob: 8 = 4 (resolve a Gaussian of exponent p to XC-grade cleanliness needs
  Ecut≈4p, measured by the negCharge probe: 2p leaves −0.77 e of negative lobes, 4p is clean) × 2 (the
  density is a product of two orbitals: tightest product exponent = 2·α_max).
- **Raster:** `AutoGrid` = the smallest N resolving the DIFFERENCE set {G−G′} alias-free (algorithm),
  then padded to powers of two (`FFTGrid`) — a constraint of our radix-2 FFT engine, *not* of FFTs
  (CP2K's FFTW-class mixed-radix runs N=36=2²·3² at full efficiency; at Ecut=160 we pad to 128³ = 45×
  CP2K's points).
- **Consumers:** collocation target, G-space Poisson (4π/G²), XC quadrature raster, ladder level 0,
  the fit-basis factories (rows 2–3).

**Row 2 — the CD (ρ) fit basis `{G}_ρ`** (`GPW_IBS::CreateCDFitBasisSet` → `PlaneWaveFit_IBS`).
- **Policy: `{G}_ρ` == the reference ball (row 1), by construction.**  The factory on the IBS is the
  seam (the evaluator no longer decides — since the 2026-07-20 fit-grid thread-through,
  `MakeRepulsion3C/MakeOverlap3C(c)` build the tensor over the REQUESTED fit basis's grid; the old
  silent override is fixed and `GPW.SiliconGammaConverges` et al. gate it).
- **Is "FFT grid == CD fit basis" justified?  Split answer (this session's measurement):**
  - For **Hartree** — yes.  The Poisson kernel is diagonal in G; projecting ρ̃ onto the ball is a
    variational truncation in the Coulomb metric, exponentially convergent, adjoint-exact.  This is the
    legitimate "Ecut is a resolution dial" case (the G-direction pin).
  - For **XC** — ALSO yes, now MEASURED (the §0f increment-0 probe, 2026-07-21): tripling the ball
    (160→480 Ha, nG 16145→83659, raster fixed 128³) moved the converged NaF energy by **0.15 mHa** —
    the ball-truncated ρ is converged at the auto calibration, and the once-suspected Gibbs/ball
    mechanism is FALSIFIED.  (The 4.26 Ha that motivated the suspicion decomposed into a MOM-pinned
    excited state + a basis mismatch; GPW == CP2K to 0.45 mHa on the same basis — §0f.)  The
    raw-raster-XC design item is DEAD; this row is a validated policy.

**Row 3 — the Vxc fit basis `{G}_vxc`** (`GPW_IBS::CreateVxcFitBasisSet`).
- **Policy:** `relCutoff · {G}_ρ` where `relCutoff` comes from the FUNCTIONAL
  (`GridCutoffFactor()`, algorithm): LDA == 1, so today `{G}_vxc` IS row 1's grid; a GGA wants > 1
  (∇ρ bandwidth) and is guarded out (`assert(relCutoff<=1)`) until the denser grid is wired.
- The v_xc that reaches the KS matrix is the fine-raster v_xc spectrally restricted to each ladder
  level and gathered analytically per pair (the integrate-back) — the fit ball is the transport
  representation.  If lead (1) above lands, this row's role gets re-examined together.

**Row 4 — the REL_CUTOFF multigrid ladder** (`GPW_Evaluator::BuildLevels`, cached per block).
- **Levels (all algorithm):** L0 = the reference grid (row 1); coarser levels at ×¼ Ecut down to the
  most-diffuse-pair floor `Ecut·α_min/α_max`, each kept only if its SPACING still resolves the sharpest
  pair assignable to it (`h²·p_max ≤ 1`, the Poisson/trapezoid bound); plus a TOP COMPLETION RUNG at
  `RelCutoffSafety()·Ecut` appended last, gated on the energy calibration (skipped when the reference
  already sits at/above `RelCutoffSafety·cutoffFactor·α_max`).
- **Pair→level assignment (algorithm):** coarsest level with
  `ecut_l ≥ relScale·kRelSafety·ecut_L[0]·(αᵢ+αⱼ)/(2α_max)`; `kRelSafety=2` (charge-vs-energy
  calibration), `relScale=1` for the smooth fields.  CP2K's rule is ABSOLUTE:
  `needed = (αᵢ+αⱼ)·REL_CUTOFF` (30 Ha in the NaF deck) — ~7.5× stiffer than ours for that system.
- **Verification-only env knobs** (never interfaces): `GPW_MGRID_ECUTS` (explicit sub-level list —
  reproduces CP2K's progression-3 ladder), `GPW_RELCUTOFF` (the absolute CP2K rule).
- **Consumers:** density collocation + integrate-back (full ladder), local PP (row 5, base sub-ladder).

**Row 5 — the local-PP integration "grid"** (`MakeLocalPP`) — *why a grid at all?*
- It is NOT a separate mesh: it is the SAME ladder objects (row 4, base sub-ladder = no top rung), used
  as the quadrature of the FIELD.  V_loc's form factor is assembled in G-space, band-limited to each
  level's ball, inverse-FFT'd to that level's raster, and each orbital pair gathers it ANALYTICALLY on
  its own level (only V is ever sampled, never the orbital product) — the same integrate-back seam as
  V_H/v_xc.
- **Where r_loc enters (the user's instinct is right):** V_loc is spectrally BROAD (the erf/r deep well,
  width r_loc), so the smooth-field assignment under-resolves it; the stiffening is `relCutoffScale=3`
  (default; was 6, Q1 re-calibration) via `GPW_LOCALPP_SCALE`/`GPW_LOCALPP_FULL` for scale-convergence
  checks.  So r_loc is honoured only INDIRECTLY (a calibrated scale) — the honest endpoint (§0e-PP)
  removes this grid usage entirely: analytic short range (built, dormant) + Gaussian-core-charge long
  range folded into the Poisson solve, which is exactly why **CP2K has no local-PP grid to match**.

**Row 6 — the mesh-path KB quadrature** (`PPMeshParams`: uniform cell mesh, eCut = `densityEcut`).
- Fallback only, for a separable-PP model that does not expose the closed Gaussian form
  (`BetaGaussian`); HGH/GTH models all do, so Si/NaF use the ANALYTIC KB (no mesh at all).  Gate:
  `GPW.AnalyticSeparablePPMatchesMesh` (4.6e-11).

**Row 7 — KMesh** (BZ sampling; `MakeKMesh(shift)`).  A k-grid, not a spatial grid — listed only
because it is the one remaining `(points, weights)` set a user chooses (N×N×N + shift).  Γ-only for the
NaF work.

## 2. Direct answers to the story's questions

- **"Is FFT-grid==CD-fit fully justified?"** — For Hartree yes (variational ball, diagonal kernel);
  for XC it is the current lead suspect for the CP2K gap (see row 2).  The fit-basis FACTORY is the
  right owner of that policy, and since the thread-through fix it is actually honoured downstream.
- **"Why a grid for the V_local integral? Same as the FFT grid? Should it be?"** — It is the ladder
  (same objects), used to sample only the FIELD; it should ultimately be NO grid (analytic split,
  §0e-PP), at which point the question dissolves.  You are right that r_loc should size it — today it
  only does so through the calibrated `relCutoffScale`.
- **"Ecut for what?"** — see §0: four meanings, one knob (`densityEcut`), everything else derived.
- **"No fake −39 basin if grids are right?"** — supported by the evidence so far: the basin has only
  ever appeared on under-resolved-ρ configurations, and the CP2K-stiff assignment even removed the
  aliasing flattery at Ecut=40.  The remaining CP2K-vs-us cost asymmetry is the RASTER (radix-2 pow2
  padding: 128³ vs CP2K's mixed-radix 36³ at the same, validated answer) — an efficiency lever, not an
  accuracy one (the §0f probes).
- **POSTSCRIPT (2026-07-21, end of session): the grid system is VALIDATED.**  With the basis matched
  (SR2) and aufbau occupations, GPW == CP2K 2026.1 to **0.45 mHa** on NaF at matched grids (and both
  SCFs converge cleanly — the old misery was the near-singular SR basis).  Every suspicion raised in
  the user story above was either confirmed-and-fixed earlier (the ignored fit basis), answered
  benign (FFT==CD-fit, the V_loc ladder), or falsified by measurement (the ball/Gibbs story).  The
  honest open items are basis completeness-vs-conditioning (GPWPlan §1) and the two SCF-strategy
  lessons (MOM cross-grid guard; the εH/εL line masking a non-aufbau hole).
