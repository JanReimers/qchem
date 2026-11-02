# C2 — the grid ε window: ONE accuracy target, every grid tolerance derived (design note, 2026-10-30)

Lane A, chunk C2 of `doc/QueueCampaignPlan.md`.  **Design first, no code yet** (controller, 2026-10-30).  The user's rulings are in
`QueueCampaignPlan.md` §7 and are NOT re-argued here: AliasFree as the raster default *if* it costs ≲ 1.3× wall; R1.0r = the
INTERSECTION group; §K = one tolerance-sized ball for ρ AND v_xc plus an XC-raster audit; ONE user-facing accuracy target from
which every internal tolerance is derived and printed.  Added 2026-10-30: R1.5 (the mixer's G-space, V-MIXER-GBALL) moved here
from B3.  Evidence: `Records/GridSharpTerms.md` (which terms are sharp, the 2026-10-25 measurements, the audit) and
`Records/GPWGrids.md` (the grid inventory; two of its defaults have since changed: C 8 → 2, raster → `BallOnly`).

This note is a plan under execution: it lives at the top of `doc/` until C2 lands, then becomes a record.

## 0. Where the tolerances are today (code facts, 2026-10-30)

| knob | where | default | what it bounds | sized from |
|---|---|---|---|---|
| `densityEcut` / `cutoffFactor` C | `GPW_IBS` / `GPW_Evaluator` ctor | auto = C·α_max, C = 2 | the reference ball (ladder L0) | the basis only (not the PP: GridSharpTerms §3's latent case) |
| `RasterPolicy` | `GPW_IBS` default `BallOnly`; `GPW_Evaluator` member default `AliasFree` (two defaults!) | `BallOnly` in production | raster N per level: G_N² = 2Ecut (BallOnly) or 8Ecut (AliasFree) | a name, not a tolerance |
| `ladderFactor` | `GPW_Evaluator` ctor | 4 | the multigrid level spacing (Ecut ×¼ per level) | constant |
| pair → level rule | `PG_Cart_MnD::PairLevel` | `kRelSafety`·Ecut₀·max(p, β_field)/(2α_max), `kRelSafety` 2 | which level a pair is collocated/integrated on | calibrated constant (charge, then Si SR ±5 mHa) |
| `fieldSharp` | `GPWTolerances` | 2/3 | β_field = fieldSharp·α_max in the rule above | derived (v_xc ~ ρ^{1/3}), pinned by ρ_lost/N |
| `vlocEps` | `GPWTolerances` | 1e-5 | V_loc-long ball: harmonic rule 2 ln(1/ε) pβ/(p+β) | **already an ε** |
| `localPPRelCutoff` κ | `GPWTolerances` | 30 Ha | local-PP sweeps: absolute rule κ·max(p, β_PP), tail e^{−κ/2} | **already an ε** (e^{−15}; CP2K REL_CUTOFF) |
| `screenEps` | `GPWTolerances` | 1e-10 | analytic 1E lattice sums + collocation box reach | an ε (prefactor drop) |
| `densityEps` | `GPWTolerances` | 1e-10 | collocation task list (the screener floor) | an ε |
| `relCutoff`, `mgridEcuts` | `GPWTolerances` | 0, {} | CP2K-parity overrides | stated choices — **stay** |
| Vxc fit ball | `GPW_IBS::CreateVxcFitBasisSet` | relCutoff_functional·{G}_ρ, = {G}_ρ for LDA | v_xc transport + **today the mixer's G-space** | the functional |

Two of nine already follow the ε principle (vlocEps, κ); two are magnitude screens (screenEps, densityEps); the RASTER — the one the
2026-10-25 failure was in — has no tolerance at all, only a name.

## 1. The one target and the derivation

**The target.**  ε = the error the user accepts **per Hamiltonian matrix element** (Ha; user ruling 2026-10-30: the name is ε, the user's own -- τ is the symmetry code's fractional translation).  Deck `solid.accuracy`, code member `accuracyEps` (beside the per-knob `screenEps`/`densityEps`/`vlocEps` it derives) (proposed default
below, §6).  Every mechanism converts ε into its own internal tolerance by a STATED formula, and the run banner prints ε and every
derived value on one line (`[accuracy] ε=… → raster G_N²/p=…, ball …, κ=…, screenEps=…, densityEps=…`).  Each formula is judged by
the error it DELIVERS on its term, measured (§5), not by its form: a term far tighter than ε is wasted cost, one looser is a defect
(user ruling: "no point in one term accurate to 1e-10 and another only to 1e-4").

**The scale constant.**  Every rule below has the shape *tail amplitude × term scale ≤ ε*, i.e. L(ε) with S_term·L³e^{−L} = ε (M2b; first proposed as `ln(S_term/ε)`) where S_term is the
size of the term's field at its sharpest (the factor that turns a relative spectral tail into Ha).  Proposed scales — each to be
CALIBRATED once against the closed-form oracles (`GPW.TightPrimitiveVlongOracle` / `…CollocationOracle`, M2 below) and then fixed:

| term | S_term (proposal) | check against GridSharpTerms §3 |
|---|---|---|
| V_loc long | the erf well depth Z_ion·2√(β/π) of the sharpest species | Mn q7, β 1.22: S ≈ 19 Ha; BallOnly Ecut/p ≈ 1 → tail e^{−2} ≈ 0.14 → 2.6 Ha predicted, ±1.5–2.7 measured ✓ |
| Hartree | V_H's depth ≈ the same well scale (the density's own) | to calibrate (M2) |
| XC (uniform raster) | max v_xc ≈ ρ_max^{1/3} at the sharpest atom | to calibrate (M2) |
| local-PP sweeps | the short-range C₁…C₄ magnitude | κ = 30 ⇒ ln(S/ε) = 15 already |

**RULED 2026-10-31 (user, via the controller): the RIGOROUS bound and the HEURISTIC calibration are isolated in code** -- "the
rigorous math aspect well isolated from the heuristic part; the latter could easily become an advanced/expert user knob."
1. **Rigorous (the bound).**  Per term class, the oracle-fitted tail model (S_term, n) and L(ε) = the root of S·Lⁿ·e^{−L} = ε: pure
   math, gated by `GPW.RasterAliasCalibration`, no fudge inside it.  At η = 1 it guarantees the single-pair worst case.
2. **Heuristic (the calibration).**  ONE explicit factor per term class, η_term = how far the delivered error undershoots the bound,
   MEASURED by `gridAudit` on the real systems as η = delivered/bound ≤ 1; the sizing uses **L(ε/η)** (corrected 2026-10-31: η·Bound(x) ≤ ε ⇔ Bound(x) ≤ ε/η -- a measured η < 1 LOOSENS the raster; the first text said L(η·ε), the wrong direction).  Its defaults live in ONE place with a comment citing the
   gridAudit measurement that set them, exposed as an expert deck key (e.g. `solid.accuracyCalibration`) defaulting to the measured
   value; η = 1 = "rigorous mode" (geometry/forces validation runs ask for it).
3. **Reporting.**  The `[accuracy]` line prints L_rigorous AND L_used with η; the grid audit reports the delivered error against both,
   so a material where η is too optimistic shows at once.  Nothing in the sizing code may mix the two (no constant buried in a formula).
4. **η by material class (user 2026-10-31: "different values of η for insulators vs metals, or covalent versus ionic"): MEASURE FIRST.**
   The M2c gridAudit set spans the classes -- Si k222 (covalent), NaF (ionic), MnO / β-MnO₂ (magnetic TMO), Al fcc k222 smeared (metal)
   -- and reports η per term per system.  Flat within ~2× ⇒ one default per term class.  If it varies, prefer a SELF-CALIBRATING η (a
   pilot gridAudit at the seed or after a few SCF iterations measures THIS run's η and sizes the final raster) over a material taxonomy
   (a metal/insulator/ionic lookup is a classifier that can be wrong); the expert key still overrides.

**RULED 2026-10-31 (user, via the controller): ε follows the RUN'S OUTPUT -- one ε, and the user never picks a metric.**  A single-point
ENERGY run sizes the rasters on the energy error (the pilot fits each raster term's scale on |ΔE_T|, target ε × the atom count); any run
that PRODUCES forces, a geometry or a surface (a `SolidSession`, PES `ConvergeAt`, a relaxation, an explicit forces request) sizes on
max|ΔH_ij,T| automatically (target ε per matrix element; verdict 3: the energy error cancels inside Tr D ΔH).  It is a capability the run
ANSWERS, named for the consumer (`ClientConsumes {Energy, MatrixElements}`, `SolidCalculation::Consumes()`: a standalone calculation
answers Energy, a point a `SolidSession` builds answers MatrixElements), not a deck flag; an expert override can come later if asked.
Tests' default: ε = 1e-5 Ha per atom (energy) and 1e-5 Ha per element (H-graded).  η = 1 (stated, no pilot) stays the rigorous mode.
**RULED as `doc/Pins.md` pin 31** (user, 2026-11-01; the wording as proposed): *"One accuracy target, judged by what the client consumes.  Every GPW
raster is sized from the run's one ε by the error it DELIVERS on the measure its client consumes -- a single-point energy on each raster
term's |ΔE| ≤ ε per atom, anything that produces forces, a geometry or a surface on max|ΔH_ij| ≤ ε per element (the energy error cancels
in Tr D ΔH: MnO V_long 0.26 Ha in H against 5 μHa in E).  The run answers the question; the user never picks a metric.  The rigorous
bound (the fitted alias tail) and the measured calibration (the pilot audit) stay isolated in code, and the XC raster, which has no bound,
is reported as MEASURED.  Earned by NaF: BallOnly's −941 μHa was the XC raster, and the F site's 'curvature' was the grid's egg-box."*

**The derived tolerances.**  **L(ε) is the root of S_term·L³·e^{−L} = ε** (MEASURED, M2b: the alias tail of an s+p+d pair carries a
Hermite prefactor ~ x³; S_Vlong = 2.35 Ha, S_ρ = 0.35 ⇒ L(1e-5) ≈ 21.6 / 19.4).  ⚠ NOT the bare ln(S_term/ε) (≈ 14.5 at 1e-5), which
under-sizes every raster by ~7 units of x; the bare form survives only where a row says so (the balls, row 3, are pure Gaussians in G):

1. **Raster N, per multigrid level** (the one that failed).  A pair p sampled on a raster of inscribed Nyquist radius G_N aliases
   with amplitude e^{−G_N²/p} (GridSharpTerms §2).  So level l's raster is the smallest FFT-friendly N with
   **G_N,l² ≥ p_eff,max(l) · L(ε)**, p_eff = max(p, β_field) of the tightest pair ROUTED to l.  Today, on the finest level at C = 2 (its
   tightest pair p = 2α_max = Ecut): BallOnly G_N² = 2Ecut = 2p (alias e^{−2}: the ±2 Ha of GridSharpTerms §3), AliasFree 8p (e^{−8}).  ⚠ CORRECTED 2026-11-01 (lane A): "AliasFree 8p" holds only for an UNROUTED pair.  AliasFree sized G² = 4·p_h·ln(1e6) with p_h = pβ/(p+β); under HartreeXC routing (β = p/3) that is ≈ 13.8p, and the 5-smooth FFT menu rounds it UP (β-MnO₂ p5f3_dD: x ≈ 16 on 90³).  A deck that states x = 8 as "the old AliasFree" is half its alias exponent (that deck's first iterate was +4.2e4 Ha; OpenWork_History5 "Dense Mn d midpoints").
   ⚠ The S·e^{−G_N²/p} bound is PESSIMISTIC away from small G_N²/p: it predicts 19·e^{−8} ≈ 6e-3 Ha for AliasFree at the floor and
   19·e^{−11} ≈ 3e-4 Ha at Ecut/p = 5.5, against ≤ 3e-4 and ≤ 1e-5 Ha measured (the alias image carries the pair's overlap with the
   field, not the field's full depth).  So M2 calibrates an EFFECTIVE L per term from the oracles; the form stays.  On every level
   whose tightest routed pair is diffuser than the level's own ball, the ε-sized N is SMALLER than AliasFree's.  **`RasterPolicy {BallOnly, AliasFree}` is deleted**; the grid audit becomes an assert that each level met its own row.
2. **Pair → level routing.**  The coarsest level whose raster meets row 1 for THIS pair: G_N,l² ≥ p_eff(i,j)·L(ε).  It replaces
   `kRelSafety`·Ecut₀·p/(2α_max) (a calibration tied to C) with the same alias condition the raster is built to, so the ladder and the
   routing cannot disagree.  The ABSOLUTE rule κ·max(p, β_PP) is this rule with L = κ/2 — the local-PP sweeps keep their form,
   **κ = 2·L(ε)** (30 Ha ⇔ L = 15).  The level spacing (ladderFactor 4) stays a cost choice, not an accuracy one: any spacing works
   because each pair goes to a level that meets ITS requirement.
3. **Balls.**  A pair p against a field β sums e^{−G²/4p_h}, p_h = pβ/(p+β) (β = 0: p_h = p).  Ball edge **Ecut_b ≥ 2 p_h · L(ε)**.
   V_loc long already uses exactly this with `vlocEps` ⇒ **vlocEps = ε/S_Vlong**.  Hartree's 4π/G² suppresses its tail by a further
   1/G_b², so its ball may stay looser (report, as today; judge it once M2 shows by how much).  §K: XC uses the SAME ball (§3).
4. **Magnitude screens.**  `screenEps` drops r-space lattice-sum terms whose Gaussian prefactor is below it; the number of dropped
   terms in the tail shell grows like the shell volume, but the reach grows only as √ln(1/ε), so a tight value is nearly free:
   **screenEps = ε·10⁻³** (proposal; today 1e-10 ⇔ ε = 1e-7).  `densityEps` (the collocation floor) is a density, so its error on a
   matrix element is ε_ρ × |v|max × the point weight: **densityEps = ε / S_term** of the field it is integrated against.  Both are to be
   confirmed by the delivered-error measurement (M2): a screen is a cost lever only where it binds.

## 2. V_loc long: analytic, or folded into Hartree

Both routes remove the term with the largest measured alias error (GridSharpTerms §3) from the raster.

- **(a) CP2K's route — fold the ion Gaussians into the Poisson solve.**  Add Σ_A −Z_A(β_A/π)^{3/2}e^{−β_A|r−R_A|²} to ρ̃ and solve
  ONE Hartree problem for the neutral total; V_loc long is then part of V_H.  Needs: the Gaussian self-energy Σ_A Z_A²√(β_A/2π)
  subtracted, and the SHORT-range ion-ion overlap correction Σ_{A<B} Z_AZ_B erfc(√(β_Aβ_B/(β_A+β_B))·R_AB)/R_AB (CP2K's E_ovrl) —
  both analytic; the Ewald/E_αZ bookkeeping changes with it (pin 29 territory: image lists by PAIR distance).  Cost: ~free (the
  Hartree raster exists).  Accuracy: the combined ρ̃ is as sharp as the sharpest ion Gaussian β_max, so row 1 above must also see
  β_max as a "pair" — the latent GridSharpTerms §3 case becomes an explicit row.
- **(b) Analytic Ewald split.**  erf(√β r)/r = [erf(√β r) − erf(√η r)]/r + erf(√η r)/r: the first part short-ranged (an analytic
  real-space lattice sum, like the GTH short part already is), the second so smooth that Σ_G ρ_ij(G)V(G) needs only a tiny ball with
  ANALYTIC pair transforms — no raster at all.  Cost: a new integral family (erf-attenuated nuclear attraction; libcint has it).
  Gain: exact force derivatives from this term (no egg-box).
- **RULED 2026-10-30 (user): route (a) first** -- "so far CP2K has not steered us wrong".  M4 measures its delivered error and wall
  at the C2 raster; route (b) waits until forces need exact derivatives from this term.  Either way Hartree and XC stay on the raster.

## 3. §K: the XC ball IS the pair ball, plus an N-vs-2N XC-raster audit

- ⟨i|v_xc|j⟩ = Σ_G (φ_iφ_j)*(G) v_xc(G), and the pair spectrum is ≤ ε/S outside the density ball (row 3), so v_xc components outside
  it do not reach H_ij: a v_xc fit on the density ball is perfect for the Fock matrix to ε.  `GPW_IBS::CreateVxcFitBasisSet` already
  returns {G}_ρ for LDA (GPWGrids row 3); C2 makes "the XC ball = the density ball" the stated rule, and a GGA's relCutoff > 1 becomes a
  RASTER requirement (∇ρ bandwidth on the sampling), not a bigger fit ball.
- **The real knob is the XC raster** (v_xc(ρ(r)) is not band-limited: ρ^{1/3} where ρ is small).  Audit, once per run at the first
  converged density (and on demand): evaluate E_xc and the v_xc matrix on N and on 2N, print max|ΔH_ij| and ΔE_xc beside ε; if
  max|ΔH_ij| > ε, a consumer that needs ε escalates (P3) or refuses with the number.  Cost: one XC evaluation on 8× points per run.
  QE has the same parameter (its default XC grid was 107 μHa off on β-MnO₂, CLAUDE.md).  The Becke route stays the alternative
  quadrature, with its own radial/angular audit (not this one).

## 4. R1.0r (the intersection group) and R1.5 (the mixer's G-space)

- **R1.0r.**  `DetectPointOps` returns the crystal group; `FoldGrid` already computes which ops leave the k-mesh invariant.  Expose
  that test and star-average ρ under crystal ∩ mesh.  Γ and isotropic meshes: no change (nothing banked moves); anisotropic meshes
  (Si k211-class gates, any `k2x1x1`) move — in THIS window, re-judged once.  Forces need it: the energy must be stationary in ρ at
  the given sampling for Hellmann–Feynman to be dE/dR.
- **R1.5 (V-MIXER-GBALL).**  The mixer's G-space is the DENSITY's ball: Kerker's metric is Hartree's 4π/G², and the mixer mixes ρ̃.
  Today it borrows the Vxc fit set.  Facts that make this cheap:
  - At LDA {G}_vxc == {G}_ρ (relCutoff 1), so swapping the source is **bit-identical at LDA** if the raster shadow comes along.
  - B3 R1.4 (`02de08b2`) put the G-space in ONE place: `KerkerPreconditioner`'s constructor (fit basis + reciprocal lattice).
  - The obstacle is the TYPE: the density fit is `FIT_CD_ABS` (the Coulomb face), while the mixer reads ρ̃ through
    `FourierDensity::GetFourierDensity(const cFIT_SF_ABS&)` and `GetRhoOnGrid`.  Design: a face named for what the mixer consumes —
    the density's **G-space representation** (its ball's G set, the reciprocal lattice, its raster), answered by the density/basis that
    owns {G}_ρ — and `GetFourierDensity`/`GetRhoOnGrid` take THAT, not a Vxc fit.  The mixer then stops importing
    `Orbital_DFT_IBS::CreateVxcFitBasisSet`.  At GGA the two balls differ and the mixer keeps the density's, which is the point.
  - Lands in the same re-bank as the ball sizing, because the density ball itself moves under row 3.

## 5. Measurements before code (each a short probe; numbers go into this note)

- **M1 — AliasFree's wall cost** on the `doc/Benchmark.md` §5a rows, serial, COPYING the §5a commands (ruling: ≲ 1.3× ⇒ AliasFree
  default; dearer ⇒ BallOnly stays default and P3 escalates).  Under §1 row 1 the question becomes "the ε-sized N's cost at the
  default ε", which is ≤ AliasFree's.
  **MEASURED 2026-10-30 (lane-a, serial, `scripts/bench` + `rundeck --set solid.raster=…`; rows and logs in
  `~/Code/qchem6-runs/{Si,NaF,MnO}/C2M1_*`):**

  | deck | iters BO/AF | wall BallOnly | wall AliasFree | × | collocate s/call BO→AF | RSS BO→AF | E_AF − E_BO |
  |---|---|---|---|---|---|---|---|
  | `bench/Si_Gamma_parity` | 16/16 | 3.85 s | 18.61 s | **4.8** | 0.063→0.334 (5.3×) | 47→71 MB | +8.5 μHa |
  | `bench/Si_k222_parity` | 15/13 | 11.82 s | 55.54 s | **4.7** | 0.060→0.365 (6.1×) | 62→138 MB | +13.6 μHa |
  | `bench/Si_k222s_parity` | 14/14 | 10.32 s | 55.87 s | **5.4** | 0.054→0.340 (6.2×) | 62→141 MB | +14.8 μHa |
  | `bench/NaF_SR2_Gamma_parity` | 16/15 | 10.57 s | 66.05 s | **6.2** | 0.171→1.242 (7.3×) | 94→213 MB | **−941 μHa** |
  | `MnO_AFM2_free` (Becke XC) | 15/15 | 7:42 | 16:49 | **2.2** | 1.35→11.88 (8.8×) | 542→1204 MB | −22.6 μHa |

  **Verdict: AliasFree costs 2.2–6.2× wall (5–9× per collocate/gather call), far above the ≲1.3× ruling ⇒ BallOnly stays the
  default and P3 escalates.**  MnO reads only 2.2× because its Becke XC sampling (~300 s, raster-independent) dilutes the raster cost.
  Iteration counts are unchanged, so the whole cost is per-call raster size.  The ΔE column is NOT an acceptance measure (pin 3); it is
  recorded for M2/M3 to explain.  The SIGN splits on the grid audit: Si's BallOnly audit is clean and AliasFree sits +9..+15 μHa
  HIGHER; NaF and MnO BallOnly print `[grid audit] UNDER` (density raster alias e^{−G_N²/p} = 1.5e-2 for NaF's p = 80 Na pair on
  36³, 9.4e-3 for MnO's p = 72 on 60³; V_loc-long 1.4e-3 / 5.0e-3) and AliasFree (audit OK) sits LOWER, −941 / −23 μHa.  NaF's
  per-term shifts (ΔEkin +5.9, ΔEen −7.1, ΔEee +1.3, ΔExc −1.1 mHa) are D redistributing on a flat surface -- see M2a.
- **M2 — the scale constants and delivered errors.**  FIRST ITEM (controller 2026-10-30): explain NaF SR2's −941 μHa (M1, E_AF − E_BO; BallOnly audit UNDER at alias 1.5e-2) -- which
  term, which pairs, and its delivered max|ΔH_ij| on BallOnly vs AliasFree.  Then: on the TightPrimitive oracles (closed forms at any exponent), sweep G_N²/p and
  the ball, fit S_term per term, then confirm on Si / NaF / MnO / β-MnO₂ that every term's delivered max|ΔH_ij| sits within a factor
  ~3 of ε (neither far tighter nor looser).
  **M2a MEASURED 2026-10-31 (NaF SR2's −941 μHa explained; runs in `~/Code/qchem6-runs/NaF/C2M2a_*`):**
  (1) it is a FIXED-DENSITY quadrature error, not relaxation: the AliasFree converged state restored into a BallOnly run gives
  E_BO at the first iterate = −24.4303944, already E_BO's own minimum (−24.4303950) to 0.6 μHa, so E_BO[D_AF] ≈ E_BO[D_BO] and
  E_AF[D] − E_BO[D] ≈ −941 μHa at fixed D (variational: second order in ΔD).  (2) it is ~92 % the XC RASTER: the same pair with XC on
  the Becke mesh (`policy.cp2kCompat false`, same deck otherwise) gives BO −24.4304034 vs AF −24.4304801, −77 μHa -- so ≈ −864 μHa
  is v_xc evaluated on the aliased ρ sampled at 36³ (density audit alias 1.5e-2 at the p = 80 Na pair), and ≈ −77 μHa is Hartree +
  V_loc long (alias 1.5e-2 / 1.4e-3).  So for NaF the binding row is §3's XC raster (M3), not the Hartree/V_long rows; and the
  §5a NaF parity gap vs CP2K (−24.4312134) is mostly ours: BallOnly +0.82 mHa, AliasFree −0.12 mHa.  (The Becke runs use the deck's
  nR = 20 mesh, which warns it is under-resolved; the SAME mesh in both, so the difference isolates the raster terms.)
  **M2b MEASURED 2026-10-30 (29c20b8a, `GPW.RasterAliasCalibration`, UTLattice_BS):** one Mn s+p+d primitive, α = 10 / 36 / 76.68, the
  raster swept (BallOnly at densityEcut = c·p, c = 0.5…3; AliasFree at the AUTO floor), x = G_N²/p read off the run's own grid audit.
  **The tail is not the bare e^{−x}**: err·e^{x} GROWS with x (the d pair's Hermite prefactor is a polynomial in G²/p).  A fit
  ln err = ln S + n ln x − x gives n = 2.3–2.6 (V_long), 3.3–3.7 (ρ̃) at EVERY α, so the scale is α-independent as §1 assumed:

  | term | n = 3 (the gate) | n = 2 | x for ε = 1e-4 / 1e-5 / 1e-6 / 1e-7 | AliasFree delivers |
  |---|---|---|---|---|
  | V_long ⟨χ\|V\|χ⟩ (Ha) | S = 2.35 Ha | S ≈ 20 Ha (= §1's Mn erf well depth 19 Ha ✓) | 18.9 / 21.6 / 24.2 / 26.8 | x 16.5–19.7: 1.3e-4 … 1.8e-5 Ha |
  | ρ̃ at low G (relative) | S = 0.35 | — | 16.6 / 19.4 / 22.0 / 24.7 | x 16.5–19.7: 4e-5 … 2e-6 |

  Gate (x ≥ 4, below it the image is not a small tail): err ≤ S·x³·e^{−x}, S_Vlong = 2.5 Ha, S_ρ = 0.4.  **Consequences for §1 row 1:**
  (i) L(ε) is the root of S·L³e^{−L} = ε, not ln(S/ε) -- at ε = 1e-5, L ≈ 21.6 (V_long) vs ln(20/1e-5) = 14.5 from the bare form;
  (ii) AliasFree (x ≈ 16–20 on the finest raster) does NOT meet 1e-5 Ha on a tight pair's V_long -- a reason, beside the cost (M1),
  for route (a) taking V_long off the raster; (iii) these are WORST-CASE single-pair numbers -- a production H_ij weights the tight
  primitive by its tiny density share, which is what M2c measures on the four systems.
  **M2c + M3 MEASURED 2026-10-31 (760d871b, the `gridAudit` postSCF action: every term rebuilt on the other raster at the
  converged D, fixed-D, no SCF step).**  η = delivered / §1's bound S·x³e^{−x}, with x = −ln(the run's own predicted alias)
  per term; V_long against S_V = 2.35 Ha, both on max|ΔH_ij| (η_H) and on |ΔE| (η_E); Hartree / XC as η_E = (|ΔE|/|E_term|)/(S_ρ x³e^{−x}), S_ρ = 0.35.
  Runs in `~/Code/qchem6-runs/<Material>/C2M2c_*.log`; script `eta.py` in the lane-A scratch (the numbers are all here).

  | system (class) | term | x | max\|ΔH_ij\| (Ha) | ΔE (Ha) | η_H | η_E |
  |---|---|---|---|---|---|---|
  | Si k222 (covalent) | V_long | 27.0 | 1.1e-10 | 5.6e-11 | 1.2e-3 | 6.4e-4 |
  |  | Hartree | 27.0 | 3.5e-9 | 8.0e-11 | — | 1.0e-2 |
  |  | XC | 27.0 | 9.9e-3 | −1.36e-5 | — | **434** |
  | NaF Γ SR2 (ionic) | V_long | 6.56 | 1.4e-2 | 9.7e-9 | 1.5e-2 | 1.0e-8 |
  |  | Hartree | 4.20 | 9.2e-2 | 6.8e-5 | — | 1.6e-5 |
  |  | XC | 4.20 | 2.1e-1 | 8.6e-4 | — | 4.5e-4 |
  | MnO AFM-II free (TMO, magnetic) | V_long | 5.30 | 2.6e-1 | −4.8e-6 | 0.15 | 2.7e-6 |
  |  | Hartree | 4.66 | 2.0e-1 | −2.4e-5 | — | 5.3e-6 |
  |  | XC | Becke mesh | 0 | 0 | — | — |
  | β-MnO₂ A1AFM U0, exact resume (TMO, magnetic; the run is AliasFree, reference BallOnly: the Δ is BallOnly's error) | V_long | n/m | 1.96 | −2.6e-6 | n/m | n/m |
  |  | Hartree | n/m | 1.34 | −5.5e-5 | — | n/m |
  |  | XC | Becke mesh | 0 | 0 | — | — |
  | Al fcc k222 kT 1e-2 (metal) | V_long | 13.65 | 6.8e-7 | −2.9e-10 | 9.7e-5 | 4.2e-8 |
  |  | Hartree | 3.41 | 6.5e-3 | −1.7e-6 | — | 6.5e-4 |
  |  | XC | 3.41 | 3.5e-2 | +8.1e-6 | — | 2.2e-5 |

  **Verdicts (rulings by the controller, 2026-10-31):**
  1. **η is NOT flat** (V_long η_H 1e-4 … 0.15; Hartree η_E 5e-6 … 1e-2: 3–4 decades, no class ordering) ⇒ η is
     **SELF-CALIBRATING**: a pilot `gridAudit` sets η per run per term.  No class taxonomy, no shipped default η per class; η = 1 stays the rigorous mode.
  2. **XC is not covered by the ρ-pair bound** (Si: the pairs alias-free at x = 27, XC still −13.6 μHa, 434× over).  v_xc(ρ)
     is nonlinear, so it aliases even when ρ̃ is band-limited.  The bound must never be presented as covering XC.  XC needs its own model: open
     design in C2; the lead is that the pilot MEASURES v_xc[ρ]'s spectral tail (an FFT on the reference raster) instead of
     predicting it from ρ's.  Until then the `[accuracy]` line labels XC "measured / heuristic", never rigorous.
  3. **For forces, η is judged on H, never on E.**  ΔE cancels inside Tr(D ΔH): MnO V_long max|ΔH_ij| 0.26 Ha against ΔE 5 μHa;
     NaF 1.4e-2 Ha against 1e-8.  An η calibrated on energies is NOT valid for forces (also a forces-design-note item, `doc/QueueCampaignPlan.md`).
  (β-MnO₂ η n/m: `gridAudit` prints the RUN's own predicted alias, not the reference's, and here the run is the alias-free one (1.1e-7 / 2.5e-9);
  the audit should report the reference raster's verdicts too -- a small follow-up.  Its V_long max|ΔH_ij| 1.96 Ha against ΔE 2.6 μHa is verdict 3's
  strongest case.)
  5. **Rigid-translation (egg-box) gate** (lane B's NaF F-site work, via the controller 2026-10-31).  Free NaF rocksalt k222, auto densityEcut
     80 Ha BallOnly (audit UNDER, collocation alias 1.5 %): translating the WHOLE crystal by 0.03 bohr moves Etot by −1.0e-5 Ha.  That is most
     of the "F-site curvature" signal (−1.27e-5 at δ = 0.03): the concave F site was largely our grid.  CP2K shows the same physics: its egg-box
     falls to 2e-8 only at ≥ 2000 Ry, and at 600 Ry it fakes a minimum.  RULE: **E invariant under a rigid translation to within ε is a NECESSARY
     gate for any grid used for geometry.**  It is cheap (one extra point), needs no oracle, and is the energy-side twin of verdict 3's max|ΔH| rule
     (the egg-box force IS the raster's ΔH seen through dE/dR).  Build it as a `gridAudit` companion (a postSCF `eggBox {shift}` re-evaluates E at
     the converged D on the shifted cell's basis, fixed-D like the audit, plus a self-consistent variant for the gate) and a test.  C2 sizes the
     defaults so a geometry deck passes it at its ε.
     **BUILT 2026-10-31 (17e089e6, `SolidCalculation::EggBox` + postSCF `eggBox {shift}`; fixed D, the same per-term table as `gridAudit`).**
     Every analytic term reads 0 to roundoff (Si Γ ≤ 6e-15, NaF ≤ 1e-13: a built-in check that the comparison is sound).  E(shifted) − E(run):
     | system (run data) | shift (bohr) | Hartree | V_long ΔE / max\|ΔH\| | XC | total |
     |---|---|---|---|---|---|
     | Si Γ densityEcut 20, uniform XC raster (test) | (0.11, 0.07, 0.03) | −8.6e-11 | −5.6e-10 / 1.1e-10 | **+4.5 μHa** / 9.5e-3 | +4.5 μHa |
     | NaF k222 free, auto densityEcut, BallOnly, Becke XC (`NaF/C2_eggbox_k222.log`) | (0.03, 0, 0) | **−8.25 μHa** | −8e-10 / **6.1e-3** | −5e-12 | −8.25 μHa |
     |  | (0.03, 0.02, 0.01) | **−12.3 μHa** | −2e-9 / 7.9e-3 | −5e-12 | −12.3 μHa |
     NaF reproduces lane B's SCF −1.0e-5 at 0.03 bohr (fixed-D −8.25e-6: second order apart, plus a different XC mesh), and LOCATES it: the
     density raster (Hartree, alias 1.5e-2 at x = 4.2), not V_long (its ΔE cancels: 6e-3 Ha in H, 1e-9 in E -- verdict 3 again) nor XC
     (Becke follows the atoms).  Si's is the uniform XC raster.  So the egg-box gate checks the same terms the audit flags UNDER, and a
     geometry deck's default must size the DENSITY raster (and, on a uniform-XC run, the XC raster) to pass it.
  4. The bound is safe (η ≤ 0.15) on every LINEAR term measured; delivered errors sit 10–10⁵× under it, which is the room the
     self-calibrated η recovers.
- **M3 — the XC raster N vs 2N** on the same four, as §3's audit would print it.
- **M4 — V_loc long route (a)**: delivered error and wall at the C2 raster vs the present term.
- **M5 — R1.0r**: list every banked test on an anisotropic mesh (the re-bank set) before flipping.

## 6. The default ε and the re-bank

- **The default ε is DECIDED AFTER M1–M4** (user ruling 2026-10-30).  Until then -- and for the existing tests after it -- every unit
  and integration test runs at **ε = 1e-5, which is today's κ = 30 Ha** (L = 15 on the local-PP sweeps), so choosing a default moves
  no anchor by itself.  A test that needs tighter states its own ε (e.g. the geometry-optimization integration tests, whose forces
  need the energy converged below the force tolerance).  Estimate for the eventual default, to be replaced by M1/M2's numbers: with
  the pessimistic S ≈ 20 Ha, ε = 1e-6 gives L ≈ 17, i.e. G_N² ≥ 17p on the finest level, √(17/8) ≈ 1.46× AliasFree's N per dimension
  (≈ 3× its points) if the bound were tight; the measured delivered errors (≤ 1e-5 Ha already at G_N² = 11p) put the effective L
  several units lower, near AliasFree's 8–11p on the finest level and below it on the coarser ones.
- **One re-bank window** (pin 10: every moved anchor re-judged against an INDEPENDENT route, never refreshed): the raster N, the
  routing rule, the ball/vlocEps, κ, R1.0r and R1.5 move together, so every anchor moves at most ONCE.  Acceptance is the grid
  convergence of ρ (ρ_lost/N; ρ under N → 2N) and of a property — never ΔE (pin 3).  Independent routes: CP2K at matched grids
  (`doc/Benchmark.md` §5a), self-convergence at 2N, and the closed-form oracles.  `scripts/anchordiff` lists the moves.
- Vendored states (`decks/states/`) are refreshed by their cold gates after the re-judgement, by the explicit `cp` they print.

## 7. Steps (each one commit or a short series)

1. M1–M5 (probes only; numbers into §5/§6).
2. ✅ BUILT 2026-10-31 (17e089e6): `solid.accuracy` ε + the expert calibration η, the rigorous `AliasTail` math isolated (GPWTolerances), the banner's `accuracy:` line (rigorous / calibrated / used; XC MEASURED), gridAudit's bound + measured η per term, the egg-box gate.  Original text: `solid.accuracy` ε + the derivation table + the banner line; at the tests' ε = 1e-5 every derived value still EQUAL to today's (a pure re-plumbing:
   bit-identical).  Unit test: the table's formulas against hand values.
3. Raster N per level from ε; `RasterPolicy` deleted (one default, the two ctor defaults gone); the audit asserts.  Routing rule from
   the same condition.  ← the anchor-moving step.
   **BUILT 2026-10-31 (lane A; the re-bank table goes to the controller before the landing sweep):**
   - **The raster is a demand, not a name.**  `GPWTolerances::{rasterX, vlocRasterX, xcRasterX}` are alias exponents: every level's
     raster reaches G_N² ≥ x·p_cap(l), p_cap = the tightest effective pair exponent the pair→level rule can route there
     (`GPW_Evaluator::PairCap`: the finest level 2α_max, a coarser one the inverse of the relative rule) -- so the routing and the raster
     are ONE condition.  x = 0 is the ball's own raster, **bit for bit the old BallOnly** (gated, `GPW.RasterDemandIsDelivered`).  The PW
     library states it as `NyquistRaster{G_N}`; `RasterPolicy` survives only on the plane-wave orbital path.  `GPWParams::raster`,
     `SolidCalcOptions::raster` and the deck's `solid.raster` are gone (the deck refuses the key, saying why).  The audit's verdict is the
     demand's and THROWS when a level is built short of it.
   - **Three modes, one resolve** (`SolidSetup::Resolve`): `accuracy` 0 = no target, the rasters as stated (`tolerances.rasterX ...`;
     CP2K's ball raster is x = 0); a stated `accuracyCalibration` derives x = L(ε/η) per term (η = 1 rigorous); otherwise (the default)
     the setup is unresolved and the run's **pilot** sizes it (`SolidCalculation::CalibrateRasters`, after the final stage): audit every
     raster term against references 8 units of x finer (XC: at least 2×), fit each over-target term's scale under the rigorous shape
     (`GridAccuracy::MeasuredDemand`, C x³e^{−x} = target), rebuild on the raised rasters, re-converge from the audited D (a restored stage:
     the Ladder on its GDM rung at kT 0), repeat (≤ 4 rounds).  A `SolidSession` FREEZES its first point's settled rasters for the rest of
     the surface (`SizedSetup`).  Fingerprint: the demands + ε (a BallOnly state is x = 0, so old BallOnly states still resume exactly).
   - **XC is measured, and it is NOT the pair shape** (verdict 2 confirmed in the pilot): Si Γ densityEcut 20, uniform XC, H-graded --
     max|ΔH_xc| 1.0e-2 / 6.5e-3 / 2.2e-4 / 3.1e-5 at x 21 / 30 / 68 / 150, a POWER of x, never reaching 1e-5 per element at any affordable
     raster (x 137 is ~16× L0's points).  The pilot fits XC as a two-point power law (`MeasuredPowerDemand`), caps it at 2× x per round,
     and on a geometry run that does not settle says so with the advice: the Becke mesh follows the atoms.  So a **geometry run wants the
     Becke XC mesh** (the NaF path test already runs on it); a uniform-XC geometry run carries a measured egg-box.
   - **Three corrections found by the first measurement sweep (2026-10-31):**
     (1) a stated raster is ISOTROPIC (N = max over axes, as the policy rasters were): per-axis N broke the imposed raster star-average on Si's
     2×1×1 supercell (`GPW_Si.Γ_Imp_eqK211`: the reference read XC 2.3 Ha off and the re-sized run converged to −28 Ha);
     (2) the uniform XC's error is the DENSITY raster's below the density's rigorous x (L(ε) = 19.35): Na k222s at density x 6 read
     |ΔE_xc| 1.2e-4 with its XC raster at x 28 and at x 60 alike, and 5e-8 once the density raster went to x 16 -- so an over-target XC term
     raises the density raster there, and its own XC raster (the power law) only above it;
     (3) a restart onto an unresolved setup ADOPTS the saved state's rasters when the state was sized for the same ε (the pilot already ran):
     without it an exact resume restarted on the ball's rasters and re-sized (`GPW_NaF.Γ_Imp_eqExactResume_Panel` 7.7e-5 off).
   - Not in this step (follow-ups): the thread-count invariance of the collocation scatter (named above: it rides this re-bank only if
     it lands in the same window); the ball/vlocEps/κ derivations (step 4); the XC raster model proper (step 5, lane B's area).
   - **M5 (folded in, 2026-10-31): the R1.0r re-bank set** -- the banked tests on an anisotropic mesh: `GPW_Si.k211_Imp_Anchor`,
     `k311_eqCplx`, `k311_Uni_MOM_eqCplx`, `k211_Uni_Cplx_ChiRealSpaceEqSupercell_Long`, and the k211 lattices of `src/Response/tests/Kernel.C`
     (their `ResponseKernel` cases).  They move in step 6, not here.
   **Also here (controller, 2026-10-30): thread-COUNT invariance of the collocation scatter** -- a FIXED accumulator count
   independent of `QCHEM_OPENMP_THREADS` (5ec26b06 made it reproducible at a fixed count only; OpenWork §3 "The short tier is NOT
   green at QCHEM_OPENMP_THREADS=4").  It moves the serial bits once, so it rides this window's re-judgement rather than its own.
3b. **The XC route follows the consumer (user ruling 2026-11-01) -- BUILT 2026-11-01 (lane A, 6713e13a).**  A geometry point
   (`SolidCalculation::AtGeometry`, `ClientConsumes::MatrixElements`) whose XC is Auto takes the Becke mesh (pin 30: H_R = {e}); a stated
   mesh and `policy.beckeXC` are obeyed.  The pilot then SIZES it: `BeckeXCParams()` (40, Lebedev 29) is calibrated on the energy, so on an
   H consumer each round audits the XC at fixed D against the next rung of a NESTED ladder (`qcMesh::BeckeLadderAbove`: (40,29) (60,35)
   (80,41 GL) (100,47 GL) (120,53 GL); `SolidCalculation::XCMeshAudit`) and climbs while the geometric tail estimate d/(1 − r) of the rung
   differences exceeds ε.  The egg-box cannot gate it (the mesh moves with the atoms); the ladder is the gate.  MEASURED (max|ΔH_xc|, Ha):

   | system | 40/29 vs 60/35 | 60/35 vs 80/41 | 80/41 vs 100/47 | settled | audit vs 120/53 | pilot cost |
   |---|---|---|---|---|---|---|
   | Si Γ, auto floor | 1.8e-5 | 8.9e-6 (tail 1.8e-5) | 4.7e-7 | 80/41 GL | ≤ 1e-5 (gated) | 3 rounds, 9 its |
   | NaF Γ, SR2 | **1.0e-4** | 1.8e-5 (tail 2.2e-5) | 2.3e-6 | 80/41 GL | 2.7e-6 | 3 rounds, 5 its |

   The ladder does not halve per rung: it drops ~20× past 60/35, so the tail estimate over-sizes by one rung on Si (60/35's true error
   ~9.3e-6) -- the safe side, at a GL-41 mesh's cost.  A STATED Becke mesh is audited and reported, never raised; a bare `cellKind Becke`
   is the MeshParams struct default (30, Lebedev 5: 1.8e-2 per element on Si) -- OpenWork §3.  Gates: `SolidSession.GeometryRunPicksAndSizesTheBeckeMesh_Long`,
   `GPW_NaF.Γ_GeometryRunSizesItsBeckeMesh_Long` (~520 s each, serial).  Not done: energy runs keep their calibrated mesh (unaudited); a
   radial-vs-angular split of the ladder (one combined rung today).
4. Balls + vlocEps + κ + screenEps + densityEps from ε (whichever M2 shows bind).
5. §K: XC ball = density ball as a stated rule; the N-vs-2N XC audit.
6. R1.0r intersection.  7. R1.5 mixer G-space face.  8. V_loc long route (a) (ruled 2026-10-30; M4 sizes it).
9. The re-bank: one sweep, every moved anchor re-judged, vendored states refreshed.  Then a pin: "one accuracy target; each
   mechanism's tolerance derived from it and judged by its delivered error".
