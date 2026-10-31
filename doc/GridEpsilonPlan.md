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

**The derived tolerances.**  **L(ε) is the root of S_term·L³·e^{−L} = ε** (MEASURED, M2b: the alias tail of an s+p+d pair carries a
Hermite prefactor ~ x³; S_Vlong = 2.35 Ha, S_ρ = 0.35 ⇒ L(1e-5) ≈ 21.6 / 19.4).  ⚠ NOT the bare ln(S_term/ε) (≈ 14.5 at 1e-5), which
under-sizes every raster by ~7 units of x; the bare form survives only where a row says so (the balls, row 3, are pure Gaussians in G):

1. **Raster N, per multigrid level** (the one that failed).  A pair p sampled on a raster of inscribed Nyquist radius G_N aliases
   with amplitude e^{−G_N²/p} (GridSharpTerms §2).  So level l's raster is the smallest FFT-friendly N with
   **G_N,l² ≥ p_eff,max(l) · L(ε)**, p_eff = max(p, β_field) of the tightest pair ROUTED to l.  Today, on the finest level at C = 2 (its
   tightest pair p = 2α_max = Ecut): BallOnly G_N² = 2Ecut = 2p (alias e^{−2}: the ±2 Ha of GridSharpTerms §3), AliasFree 8p (e^{−8}).
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
  | ρ̃ at low G (relative) | S = 0.35 | — | 16.6 / 19.4 / 22.1 / 24.7 | x 16.5–19.7: 4e-5 … 2e-6 |

  Gate (x ≥ 4, below it the image is not a small tail): err ≤ S·x³·e^{−x}, S_Vlong = 2.5 Ha, S_ρ = 0.4.  **Consequences for §1 row 1:**
  (i) L(ε) is the root of S·L³e^{−L} = ε, not ln(S/ε) -- at ε = 1e-5, L ≈ 21.6 (V_long) vs ln(20/1e-5) = 14.5 from the bare form;
  (ii) AliasFree (x ≈ 16–20 on the finest raster) does NOT meet 1e-5 Ha on a tight pair's V_long -- a reason, beside the cost (M1),
  for route (a) taking V_long off the raster; (iii) these are WORST-CASE single-pair numbers -- a production H_ij weights the tight
  primitive by its tiny density share, which is what M2c measures on the four systems.
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
2. `solid.accuracy` ε + the derivation table + the banner line; at the tests' ε = 1e-5 every derived value still EQUAL to today's (a pure re-plumbing:
   bit-identical).  Unit test: the table's formulas against hand values.
3. Raster N per level from ε; `RasterPolicy` deleted (one default, the two ctor defaults gone); the audit asserts.  Routing rule from
   the same condition.  ← the anchor-moving step.
   **Also here (controller, 2026-10-30): thread-COUNT invariance of the collocation scatter** -- a FIXED accumulator count
   independent of `QCHEM_OPENMP_THREADS` (5ec26b06 made it reproducible at a fixed count only; OpenWork §3 "The short tier is NOT
   green at QCHEM_OPENMP_THREADS=4").  It moves the serial bits once, so it rides this window's re-judgement rather than its own.
4. Balls + vlocEps + κ + screenEps + densityEps from ε (whichever M2 shows bind).
5. §K: XC ball = density ball as a stated rule; the N-vs-2N XC audit.
6. R1.0r intersection.  7. R1.5 mixer G-space face.  8. V_loc long route (a) (ruled 2026-10-30; M4 sizes it).
9. The re-bank: one sweep, every moved anchor re-judged, vendored states refreshed.  Then a pin: "one accuracy target; each
   mechanism's tolerance derived from it and judged by its delivered error".
