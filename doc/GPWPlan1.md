# GPWPlan1 — the post-raw-XC execution plan (2026-07-23; restructured 2026-07-26)

The forward queue, superseding `doc/GPWPlan.md`'s TODO section (that file remains the authoritative RECORD of
the 2026-07 campaign, and `doc/GPWHistory.md` holds the deep archive).  **Read GPWPlan.md's "Durable pins /
invariants" before working here:** THERE IS NO CUT; no grad-student knobs (policy enums, not numeric dials);
spin-native is the formulation; correct > efficient > end-user > dev > readable.

Layout: **DONE** (condensed highlights — full detail in the cited commit messages), then the **forward roadmap**
(the agreed order of attack), then **pending** items from the original sequence, then **future considerations**.

---

# DONE (highlights — detail lives in the commit messages)

- **Platform**: raw-XC feed everywhere (XC-collapse basin removed, ρ_DM≥0 for aufbau); 0h occupation guards
  (persistent-hole → MOM release, never-silent non-aufbau warning); C=2 default; RasterPolicy BallOnly default
  on the GPW surface (≈1 mHa, CP2K's bet); BallOnly+OMP ≈ CP2K runtime parity.
- **GDM/variationality RESOLVED**: the GPW collocation E[ρ] IS variational; GDM converges it (NaF/Γ → −24.478,
  every geodesic step descends).  The "non-variational" scare was GDM's `ComputeStep` gated `if itsEn≥EMax`
  degrading to an UNMIXED diagonalize; fix = `EMax 1e-2→1.0`.  (env instruments `GPW_XCROUTE`, `GPW_GDMTRACE`.)
- **Item 2 — PER-SYSTEM ITERATION OUTPUT (both increments)**: `virtual DisplayColumnHeaders/DisplayColumns` on
  `tSCFIterator<T>`, `Molecular`/`Solid` subclasses; compact `tDensityMixer::Tag()` (Lin/Ker/Pul) + accelerator
  `Tag()`/`Count()` (Null/DIIS:N/GDM); `cfg` 1-char `*` on config change; Solid gap column (folds
  `ReportBandGap()`), virial dropped there; `ρ_lost/N` column (signed grid-charge leak per electron, via
  `EnergyBreakdown::GridChargeLost`).  Direct-min now provably disables the density mixer.
- **4b FERMI SMEARING — Increment 1 (per-block μ), `efeb251e`**: `SCFParams::SmearingkT` (0=off, zero
  regression).  `TOrbitals::TakeElectronsFermi(ne,kT)` — μ by bisection on Σgᵢfᵢ=nₑ, occ=gᵢfᵢ, Mermin
  −TS≤0.  Free energy plumbed via `EnergyBreakdown::MinusTS` + the iterator's `TotalEnergy()` helper, so the
  E-flat gate / facade / display all read A=E−TS.  (The "just a Hamiltonian term" sketch didn't fit — a term
  never sees occupations — so −TS is a WF-side scalar `GetEntropyTerm()` the iterator stamps in.)  Gates:
  `SmearingInertOnGap` (Si/Γ kT=1e-3 == −7.11506, MinusTS≈0), `SmearingConvergesDegenerateShell` (Si pseudo-
  atom degenerate 3p: aufbau never converges Δρ → smeared CONVERGED).  **Finding: kT must EXCEED the frontier
  splitting** (NaF: 1e-2 converges, 1e-3 slosh-rotates).
- **MOM-masked Fermi (`00f3b27b`)**: compose character-selection (MOM) with energy-smearing (Fermi) — Fermi on
  effective energies εᵢ+Λ(1−sᵢ)² so a diving diffuse ghost stays empty by CHARACTER while the frontier smears
  by energy.  `SCFParams::MOMSmearPenalty`.
- **Field-sharpness grid registration (`b3dad5be`) — the NaF ghost CURE** (retires the "near-gapless flapping"
  pathology): `PairLevel` resolves `max(αᵢ+αⱼ, β_field)`, β_field = 1/(2r_loc²) for the local-PP + (2/3)·α_max
  for the density/XC (V_xc~ρ^⅓).  Root cause was `63f20bd1`'s 3³ grid + the density relative rule routing
  diffuse pairs onto 27 points → over-attractive V_xc → a diffuse eigenvalue dives.  Cures NaF on DEFAULTS
  (−24.4311, 44 iters, no spike, no smearing); ρ_lost/N drops 4-5 orders; Si untouched.  `GPW_FIELDSHARP`
  overrides the 2/3.
- **Pivoted-Cholesky detector + vetting FOUNDATIONS (`6727b5af`, `9b546bc1`)** — the §4a DETECTOR: `Ortho::
  CholeskyPivoted` (rank-revealing, sparse V — SELECTS independent AOs, drops redundant); `qchem::
  PivotedCholeskyDrops(S)` (basis-neutral — eigen GAP gives the drop COUNT, pivoted Cholesky the ORDER);
  the molecular subset ctor.  Detector done; the ACTUATOR (prune/report) is open — see pending §4a.

New set of numbers.
**0. REPORTING FEATURE** — build `doc/RunReportPlan.md`.  `qchem.Reporting`: the report IS json, a generic
console renderer (layout INFERRED from json structure — table vs tree; NO per-section renderers), a global
sink (`GlobalReport` keyed + key-free `CurrentRunReport`), incremental section-by-section rendering, detail
level console-only.  Consolidates the scattered `[GPW grid]`/`[overlap S]`/cache-RAM/SCF-settings prints.
**✅ DONE (2026-07-26) — migration COMPLETE (all 4 sections)** — module + generic renderer + sink; provider
facet is `{CurrentRunReport, EmitSection(name, json)}` plus a **section CURSOR** (`Set` + RAII `Section`/`Row`)
that lets a five-layers-down provider (the `LASolver` conditioning) write context-free — no `report&`/irrep
threaded through the physics.  `scf` → `Calculation::Converge`; `basis` → `perIrrep` (conditioning) + `removed`
(`PivotedCholeskyDrops`) via the cursor; `grids` → `grids.ladder` from `GPW_Evaluator::EmitGridsReport` (`RunGPW`
brackets the GPW run, which also gave GPW `basis.perIrrep` free); `cache` → per-run snapshot of the singleton
integrals cache (`IntegralsCache::EmitReport`, cumulative, never cleared).  Rendered by `scfrun`/`RunGPW` verbose
(`SetConsole`); `ctest -j16` green (606/606).  Remaining: SOLID refactor (deferred), `basis.removed` naming, GPW
stream-cache section, disk/rolling-log SINK (deferred; see RunReportPlan "Renderer vs SINK").

**items 1–3 + IBZ — DONE (2026-07-28/29)** (detail in the commit messages):
- **1. VALGEN CLI + Al basis** — `CLIapps/valgen.C` shipped; Al q3 basis in `valence_lowq.bsd`.  (F⁻/anion
  DIIS limit-cycle + atomic-op(r) radial-only caveats recorded then; not on the metal path.)
- **2. DEGENERATE-SHELL METAL (FCC Al @ Γ)** — gates `AlFCCDegenerateShellAufbauStalls` (aufbau floors Δρ on
  the degenerate 3p) + `AlFCCMetalGlobalMu`/`AlFCCAnnealedMetal` (Fermi smearing + a kT-anneal driver
  `RunGpwAnnealed` converge it).  FINDING: **GDM ⊥ Fermi smearing** — its geodesic uses the fixed-occupation
  `[F,D]`, not the free-energy gradient, so it diverges under smearing (use DIIS/Kerker for smeared runs).
- **3. GLOBAL μ ACROSS k-BLOCKS + IBZ/k-star — done, EXACT.**  `Crystal_EC` global-μ metal mode (one μ across
  the BZ, `FillOrbitalsGlobalFermi`), and IBZ folding is EXACT (`AlFCCMetalIBZExact`: 8→3 k-points, matches the
  full mesh to ~6e-8, half the wall-clock).  Density star-average split by tolerance: Hartree on ρ̃ (G-space
  `SymmetrizeGMap`), XC on the non-negative ρ_DM raster (real-space `cFIT_SF_ABS::SymmetrizeRaster`),
  one-electron auto-projects.  Ops single-sourced from the basis + ctor-injected (no setters); one new abstract
  method + one basis accessor, nothing else touched.  Reciprocal↔direct + time-reversal knowledge lives in
  `SpaceGroup` (`ReciprocalPointOps` TR-on / `DirectPointOps` TR-off).  NON-symmorphic (diamond Fd-3m) now EXACT
  too (item 5 below, `SiDiamondIBZ_NonSymmorphic` -7.77847 vs full-mesh -7.77846): the glide τ enters the density
  star-average via `SpaceGroup::ReciprocalOps/DirectOps` (full {W|τ} point group, no guard) — G-space Hartree
  gets the e^{+2πi(Um)·τ} phase (`SymmetrizeGMap`), the real-space XC raster the exact FFT fractional shift
  ρ(W·x+τ) (`SymmetrizeRaster`); the k-fold still uses the linear TR-on ops.
- **4. MULTI-K Na metal — FIRST LIGHT** (gate `NaFCCMetalGlobalMu`): the honest half-filled-band Fermi surface
  (shifted MP + global μ + smearing, `VALENCE_LOWQ_SR2` at the real density).  Machinery validated; the
  converged-mesh / cohesive-energy "done properly" pass is what remains (below).

---

#roadmap (agreed 2026-07-26 — toward an honest metal + Fermi-smearing test)

*(Items 1–3 + IBZ/k-star are DONE — summarised in the DONE section above.  Remaining forward work:)*

**4. MULTI-K Na — the honest metal test, DONE PROPERLY** — first light is in (`NaFCCMetalGlobalMu`, above);
the remaining pass is the CONVERGED metal: a denser mesh (needs IBZ — now available), a metallic Na basis via
valgen (the SR2 gate validates the machinery, not a cohesive energy), and a physical did-E-move / bulk anchor.
Optional BCC once a `BCCUnitCell` exists (only FCC today).

**5. NON-SYMMORPHIC IBZ support — DONE** (`SiDiamondIBZ_NonSymmorphic`, -7.77847 vs full mesh -7.77846).  The
glide/screw τ-phase density symmetrization for diamond-type crystals: `SpaceGroup::ReciprocalOps()`/`DirectOps()`
expose the FULL {W|τ} point group (no symmorphic guard, no TR); the G-space Hartree star-average carries the
e^{+2πi(Um)·τ} phase (`SymmetrizeGMap`, scatter form U=Wᵀ) and the real-space XC raster the EXACT FFT fractional
shift ρ(W·x+τ) (`SymmetrizeRaster`, one FFT round-trip per distinct τ — exact even on the τ-incommensurate n=15
grid).  Unit-pinned by `GMapUT` (diamond structure-factor invariance + a W-only-corrupts control).  The k-fold
is unchanged (linear TR-on ops).  REMAINING here: the spin-polarized global μ (both spin channels, one μ).

**FIRST LIGHT DONE 2026-07-28 — committed gate `GPW_SCF.NaFCCMetalGlobalMu`.**  FCC Na, Zion=1 (3s¹) → ONE
electron/cell → half-filled band → μ cuts THROUGH it = a genuine Fermi surface.  Shifted Monkhorst-Pack 2×2×2
(`kShift=½`, k at ±¼) + global μ + smearing.  MEASURED (the textbook smeared Fermi surface): μ mid-band
(≈-0.014); the 2 k-points inside the surface fill (n_k=2.0, ε=-0.077), the 6 on it smear FRACTIONALLY
(n_k=0.67, ε≈μ, f=1/(1+e^{0.7})=0.33/spin); `Σ_k w_k n_k = 1.0000` exactly; converged ~26 iters, A=0.0455.
**Basis = `VALENCE_LOWQ_SR2` at the REAL FCC-Na density (a=10 au, matched to Na's atomic volume)** — SR2 drops
the diffuse Na s 0.0857 + p 0.05, so the Bloch overlap is well-conditioned at the CORRECT lattice constant
(cond≈38).  (First cut used `SR` at a stretched a=12 au — stretching the cell was a workaround the user rightly
flagged; SR2 at the real density is the honest fix, and 9× faster: 1.3 s.)  **Annealing composes**
(`NA_ANNEAL` in `DISABLED_NaFCCMetalExperiment`): kT {0.02→0.01→0.005}, re-seeded; internal E converges toward
T→0 — NOT kT-flat, unlike Al's degenerate shell (a real Fermi surface genuinely depends on T, so annealing
recovers the T→0 answer).  **CAVEATS / next:** (a) SR2 is a MINIMAL 6-function Na basis (no diffuse 3s) — the
gate validates the machinery + Fermi surface, not a cohesive E (A>0 here); a fuller metallic Na basis via
valgen (step 1) is the accuracy follow-up; (b) mesh convergence (2×2×2 → denser) needs IBZ/k-star reduction to
stay tractable; (c) BCC (Na's real structure) once a `BCCUnitCell` exists (only FCC today).  The machinery
(shifted MP + global μ + smearing + annealing) is validated end-to-end.

---

# Pending (from the 2026-07-23 sequence — not on the immediate roadmap)

### 1. PARAM-STRUCT GRADUATION
Move settled knobs out of env into typed structs; grade Standard vs Advanced (Advanced = sensible defaults,
touched only for pathological cases; guards self-correct so users don't hand-tune).  `SCFParams` keeps the
Standard surface; nested Advanced structs (`MixerParams`, `OccupationParams` policy enum {Aufbau, MOM, Fermi},
`GridParams`).  DELETE the `NAF_*`/`GC_*` env recipes (hard-code per test); keep verification instruments +
ops valves as documented env.  Fold in the BallOnly default.  (Partially superseded by the flat `SmearingkT`/
`MOMSmearPenalty` fields already added — the graduation would nest them.)

### 3. CACHE2/3 BYTE-BUDGET LRU (GPWPlan §5)
Byte budget + LRU eviction, policy per scope BY THE ALGORITHM (lattice = scoped size-1; molecular = generous);
per-cache RAM in the run report (→ the Reporting `cache` section).  Retires the clear-based band-aids.  Wanted
before the diffuse-basis campaigns (they bloat the geometry caches).

### 4a. DIFFUSE-BASIS ROBUSTNESS — the ACTUATOR (detector is DONE, above)
Turn `PivotedCholeskyDrops` into action so `VALENCE_LOWQ_SR`/`_LOWQ` "just work".  Two open decisions
(doc-analysed 2026-07-26):
- **Auto-prune vs report-only.**  Auto-prune = a per-type `IrrepBasisSet::Prune(indices)` + `BasisSet::Prune`
  swapping the pruned IBS into `itsBasisSets`; the vetting driver computes S, gets drops, prunes, reports the
  new λ_min; RKB analyses the LARGE sector and prunes paired large+small.  Report-only = just LIST the
  redundant functions (named by exponent/angular/atom) and ask the user to fix the basis (the 80/20, since the
  user reruns either way).  **Report-only lands in the Reporting `basis.removed` section (step 0).**
- **β_field pin + SR cost.**  Pin the density-rule β_field=(2/3)·α_max by ρ_lost/N SATURATION, not wall-clock
  (the 2/3 default gives back most of `63f20bd1`'s 3.4× by routing diffuse density pairs fine).  SR cost
  decomposes: TIME = diffuse real-space offsets (screening lever), MEMORY = grid fineness.  Grid-registration
  and Auto-trim are orthogonal in PURPOSE but compound in COST (trimming drops orbital directions, not
  collocation DIMENSION); pruning the basis (not just the orbitals) is what actually cheapens the streams.

### 4c. LIBCINT LATTICE ENGINE
Teach `PG_LibCint` to realise `Molecule::LatticeSum1E` so the GPW 1E/KB/3C build can select the faster engine
(collocation streams stay MnD).  User OOD: `PG_LibCint`/`PG_MnD` inherit a `LatticeSum1E` interface the
molecular SCF never sees (ISP); the Solid code sees only that interface, never MnD/LibCint except an enum at
the very top (LSP).

### 5. B_ij(R) k-INDEPENDENT 1E MEMO (GPWPlan 0.5(d))
LAST — payoff only on multi-k; time with the IBZ track.  Cache B(R), never M(k) — "keep k out of the key".

---

# Future considerations — Fermi smearing (captured 2026-07-26, so we don't lose them)

- **Principled kT selection.**  Today kT is a hand-set knob (finding: it must EXCEED the frontier splitting —
  NaF 1e-2 converges, 1e-3 sloshes).  Tie it instead to a physical quantity: the gap, the DOS at the Fermi
  level, or a target electronic entropy — so the software picks a sensible kT (the expert-system direction).
- **Smearing flavors beyond Fermi–Dirac.**  Gaussian / Methfessel–Paxton / cold smearing give a SMALLER
  entropy correction and better T→0 behaviour than plain Fermi–Dirac.  Decide per battery-materials need
  (Fermi–Dirac IS the true finite-T physics; MP/cold are numerically-nicer approximations for a T→0 answer).
- **T→0 extrapolation.**  The SCF minimises the free energy A=E−TS; to recover the physical internal energy E
  at zero broadening use the ½(E+A) extrapolation (or a MP-order-consistent one).  Report it alongside E, −TS, A.
- **Annealing as a general capability.**  The kT-schedule (ramp T→0 in steps, re-seed each stage from the last
  — the density-continuation machinery exists) is used in roadmap step 2, but it's a general convergence tool
  for any hard/near-gapless landscape, worth exposing as such.  (Built as the `RunGpwAnnealed` test-driver in
  item 2; promoting it to a library/facade capability is the open piece.)
- **GDM needs a smearing-aware search direction.**  Measured in item 2: the GDM direct-minimiser DIVERGES under
  Fermi smearing because its geodesic direction is the fixed-occupation `[F,D]`, not the free-energy gradient
  dA/d(rotation) (which carries an occupation-response term).  To let GDM tail-polish a smeared solid, build the
  step from the smeared gradient; until then, smeared runs use DIIS/Kerker/Pulay and GDM stays fixed-occupation.

---

# Future consideration — Becke XC grid (the last term without a near-ideal grid)

**STATUS 2026-07-30 — BUILT + GATED (uncommitted).**  The pieces, in landing order:
- **Periodic Becke mesh**: `qcMesh::UnitCellKind {Uniform,Becke}` in `MeshParams` (folded into `ID()`);
  `UnitCell::CreateIntegrationMesh` grows the Becke branch (`src/Structure/Imp/UnitCell.C`) — per-atom
  radial×angular grids, Becke cell weights with the competitor set over ALL periodic images, gathered in
  expanding Chebyshev cell shells until BOTH tails are provably < ε by the s(μ) collinear bounds (magnitude
  screen, NO radius anywhere — the pin holds).  ε=1e-8 ON PURPOSE (bounds are worst-case-collinear;
  quadrature error ~1e-4 dominates; 1e-10 measured several-fold costlier for zero gain).  Points wrapped to
  the home cell; far radial tails get ~0 weight and are dropped — a diffuse function's reach costs nothing.
  Acceptance: `LatticeMesh.Becke*` (UTStructure) — ∫1=Ω (−2.9e-4), lattice-Gaussian 1e-4, sharp α=100
  core field 1e-4 where the same-point-count uniform grid misses by >1e-2.
- **ANGULAR-RULE AUDIT** (`Mesh_Angular.*`, UTMesh): the 40-year-old hand-coded Lebedev tables are ALL
  CORRECT to their claimed degrees (6/8/12/24/30/50 → L=1/3/5/7/8/11, ≤1e-6 = the 7-figure constants);
  GaussLegendre is machine-exact at ANY L; EulerMaclaren has NO algebraic degree (degree-3 moments ~1e-3
  at L=29 m=2).  ⇒ the Becke shell (angular-limited) uses **GaussLegendre L=29** (Lebedev stops at L=11).
- **`PW_XC_Becke`** (`src/Hamiltonian/Internal/PWTerms.C`): the real-space sibling of PW_XC (one LDA
  functional per instance, Dirac+VWN pair) — ρ(r) evaluated ANALYTICALLY per mesh point off the density's
  ScalarFunction face (no FFT/collocation/fit; ρ_DM≥0 so the ρ>0 guard is inert), `⟨i|v_xc|j⟩` =
  `qcMesh::WeightedOverlap` (new tabulated-V overload) over the Bloch-summed orbital basis.  The MESH is
  built once by the caller and SHARED by the pair (as PW_XC shares its fit basis).  Hartree untouched.
- **GATE — Si Γ (`GPW_SCF.BeckeXCMatchesUniformXC_SiGamma`, enabled): PASSES.**  Same converged density,
  uniform (raw-adjoint, Ecut=60 reference) vs Becke (nr=40, GL-29): dE_xc=1.1e-4 Ha, ρ-lost=−8.9e-4,
  max|ΔV_xc|=3.5e-4.  LESSON: at the SCF-sufficient Ecut=20 the max|ΔV_xc| read 1.2e-2 — the RASTER's
  element error, not Becke's (dE_xc was already 1e-4).  The uniform reference must be matrix-grade.
- **GATE — NaF SR2 sharp-F stress (`DISABLED_BeckeXCMatchesUniformXC_NaFSR2`)**: dE_xc tiny throughout
  (1.6e-4 BallOnly / 1.3e-5 AliasFree; SCF at Ecut=160 hits the CP2K truth −24.4312 in 17-18 iters), but
  element-wise the UNIFORM raster is NOT converged on sharp-F pairs at practical Ecut: max|ΔV_xc| vs the
  same Becke matrix = 1.55e-1 (BallOnly/160) → 1.89e-2 (AliasFree/160, worst at a diagonal deep-F-core
  element) — which is this grid's reason to exist.  The committed gate is therefore Becke INTERNAL
  convergence (B40 vs 2×-refined B80) + energy-level agreement — **PASSES: max|B40−B80|=7.1e-5,
  dE_xc(B40,B80)=1.5e-5** — the Becke matrix is converged below 1e-4 INCLUDING the sharp F core, while
  the uniform raster's worst element wanders 0.16 with raster config (BallOnly@80 −1.077 / AliasFree@160
  −0.897 vs Becke −0.9154±7e-5; V_xc~ρ^⅓ is not band-limited at the core, so even exact-quadrature
  AliasFree cannot element-converge it).  Even E_xc drifts TOWARD Becke as the raster refines:
  −4.95438 (BallOnly@80) → −4.95483 (AliasFree@160) → −4.95524±1.5e-5 (Becke).
- **VISIBILITY + KNOBS (2026-07-30, user request)**: building a periodic Becke mesh now announces itself —
  a `[Becke grid]` console line (the `[GPW grid]` family: natom, radial/angular recipe, kept/dropped
  points) plus a `grids.becke` addendum in the run report (`report::EmitAt`; inert when no run is open).
  NO line ⇒ the uniform XC route is in use (e.g. `DISABLED_NaFRocksaltGamma` — Becke is gate-only today).
  Env instruments on the gate default (`BeckeXCParams`): `GPW_BECKE_L` (GaussLegendre order, any int),
  `GPW_BECKE_NR` (radial points), `GPW_BECKE_ALPHA` (MHL scale), `GPW_BECKE_ANG=lebedev` (the
  octahedral-orbit tables; `GPW_BECKE_L` then = direction count {6,8,12,24,30,50}) — explicit-args
  callers (the NaF B80 refinement probe) keep their pinned resolution.
- **ANGULAR-ORDER LADDER (Si gate, 2026-07-30 — the user's site-symmetry conjecture, measured).**  The
  T_d site symmetry makes the field content sparse (L∈{0,3,4,6,...}; pair products add ≤2·l_max), so the
  needed order is far below the partition-worst-case L=29: GL ladder dV_xc = 3.5e-4(floor)@L=17 ==
  L=29 → 9.5e-4@**L=11 (72 dirs, 4000 pts, sub-mHa dE_xc — the measured sweet spot)** → 5.6e-3@L=7 →
  1.4e-2@L=5.  Lebedev at the SAME degree: better on V_xc elements (Leb-24 3.1e-3 vs GL-7 5.6e-3 — the
  O_h-orbit cancellation is real) but 5-10× WORSE on the ρ-weighted integrals (Leb-50 dE_xc −7.0e-3 vs
  GL-11 7.8e-4): the (±1±1±1)/√3 orbit sits EXACTLY on the diamond bond axes — ρ's sharpest, most
  T_d-asymmetric angular feature.  ⇒ default stays GL-17 (safe) / GL-11 (fast); the principled endgame is
  a SITE-GROUP-ADAPTED grid whose orbits avoid special (bond) directions — upgrading the "symmetry-adapt
  this grid" icing item with concrete design guidance.
- **ROUTE SEAM WIRED (2026-07-30, user request)**: every `Ham_PW_DFT` ctor takes an optional
  `qcMesh::MeshParams xcMesh` (ctor-injected policy, no setter) — default `UnitCellKind::Uniform` is
  bit-identical to before; `::Becke` assembles the `PW_XC_Becke` pair on the geometry's periodic Becke
  mesh (built once, shared; no Vxc fit basis on that route).  BuildTerms ANNOUNCES the choice either way
  (`[XC quadrature] UNIFORM G-space raster` / `periodic BECKE atom-centred mesh` + `grids.xcQuadrature`).
  Exposed as `GpwOptions.xcMesh`; `DISABLED_NaFRocksaltGamma` carries the one-token switch
  (`o.xcMesh.cellKind = Uniform  // <-- ::Becke turns Becke ON`, recipe from `BeckeXCParams()`).
  SMOKE-VERIFIED in-SCF (NaF, 3 iters, L=11): route runs, charge conserved, ~15 s/iteration unoptimised.
- **OPEN (next increments)**: making Becke the DEFAULT for (diffuse) bases — after the perf item;
  per-SCF cost (share ρ across the pair via a provider, cache the Φ basis-values matrix per k → GEMM
  route) — today ~15 s/iteration on NaF at GL-11, fine for experiments, not production;
  spin-native (`PW_XC` itself is unpol today); symmetry-adapting the mesh points (the icing — now with
  the measured design pin: site-group orbits must AVOID bond directions); per-element radial scaling
  (one mhl_alpha serves all species today).

Framing (user, 2026-07-26): this is likely the **final step to give every Hamiltonian term a near-ideal grid
across every basis set** — diffuse included.  Per-term today:

| Term | Treatment | Ideal for diffuse? |
|---|---|---|
| Kinetic, Overlap | analytic (Gaussians) | ✅ no grid |
| Local-PP short, KB nonlocal | analytic (3-centre / separable) | ✅ no grid |
| Hartree | G-space Poisson `4πρ̃/G²` (FFT) | ✅ diffuse ρ is SMOOTH → cheap |
| Local-PP long | G-space form factor + collocation (field-sharpness) | ~ near-ideal (field-sharpness handles it) |
| **XC** | pointwise-nonlinear `V_xc~ρ^⅓`, SHARP at cores, uniform multigrid | ❌ diffuse → routed fine → **explodes** |

**XC is the one term left.** It is pointwise-nonlinear and sharp at the cores, and on the uniform multigrid a
diffuse pair × sharp field is a TWO-SCALE integrand the single-level-per-pair assignment can't split — so it is
routed to the fine grid over the function's whole (e.g. 21.5 au) reach to capture a small near-core correction
→ RAM/time explosion (see §4a).  Spatial adaptivity (fine near atoms, coarse elsewhere) has NO clean `{G}`-space
picture — it is inherently a REAL-SPACE construct — so the fix leaves the uniform/`{G}` framework for the XC.

**Design — Becke XC grid, NOT GAPW:**
- Right UnitCell creates a uniform grid.  We will need a parameter in src/Mesh/Mesh.C something like:
  enum class UnitCellKind  {Uniform,Becke}; and then add UnitCellKind to struct MeshParams.
- KEEP the uniform FFT grid for **Hartree** (the G-space Poisson is the reason GPW is fast; a Becke grid can't
  do it).  So this ADDS a second grid for the XC only; it does not replace anything.
- ADD an atom-centered **Becke/Voronoi** grid for the XC quadrature: collocate ρ on the Becke points, evaluate
  ε_xc pointwise, `⟨i|V_xc|j⟩ = Σ_g w_g χ_i χ_j V_xc(g)` — the STANDARD molecular-DFT quadrature.  Dense radial
  near each nucleus (resolves the sharp V_xc), sparse far out (a diffuse tail costs almost nothing) → the
  diffuse explosion simply does not happen (point count ~ atoms × radial × angular, NOT the diffuse reach).
- This is the scale-decomposition (fine near atoms) pictured cleanly.  It is **lighter than GAPW** (which
  augments BOTH Hartree and XC near cores + a compensation charge — deferred, out of first-pass scope): Becke
  is XC-quadrature ONLY, Hartree stays G-space.
- DISABLED_NaFRocksaltGamma might be a good test case, the F- ion tends to make sharp peaks in rho and Vxc.
- We now have full space group support, so as icing on the cake we should symmetry adapt this grid.  The other uniform grid shoudl already be symmetry adapted.

**Beyond LDA — the functional-ladder win (user, 2026-07-26).**  The Becke grid also unlocks the ADVANCED
functionals, and for the same reason.  A semi-local XC potential is a pointwise-nonlinear function of ρ AND its
derivatives — ∇ρ (GGA), τ / ∇²ρ (meta-GGA) — and on a Gaussian basis ρ, ∇ρ, τ are ALL ANALYTIC at each grid
point (from χ_i and its analytic derivatives).  On the UNIFORM grid, climbing the ladder demands FINER grids
to resolve those higher-derivative, sharper-at-the-core fields (the codebase already flags it: GGA needs
`MeshParams::relCutoff` bumped).  On a Becke grid the derivatives are exact per point and the mesh is already
dense at the cores, so **LDA → GGA → meta-GGA costs the same** — no relCutoff escalation.  Correlation
potentials (also semi-local) ride along identically.  For HYBRIDS: the EXACT-exchange fraction is the analytic
ERI machinery (`Vee`/`FittedVee`), NOT a grid — but the semi-local DFT portion (the GGA exchange + correlation
that is most of a PBE0/HSE) is exactly what Becke makes accurate + affordable.  So the Becke XC grid is the
enabler for the whole planned functional upgrade (GGA/meta-GGA/hybrid), not only a diffuse-basis fix — it makes
the grid near-ideal across all basis sets AND all functional rungs.

**Why it's bounded, not a rewrite:** the molecular Becke grid code ALREADY exists; adapting it for a unit cell
(periodic Voronoi partition + periodic BCs) lives in `qcStructure`.  Discussed and deferred weeks ago; the
diffuse-function cost is the NEW, concrete argument for it (beyond the old all-electron/GAPW motivation).

**Gate:** it lands with a cross-check — Becke XC `==` the uniform-multigrid XC (same E_xc, same V_xc matrix to
grid tolerance) on a CONDITIONED basis (Si / NaF SR2) — before it becomes the default for diffuse bases.

**Timing:** the interim (uniform multigrid + field-sharpness + vetting) keeps conditioned bases correct + fast;
the Becke XC grid is the principled endgame that makes GENUINELY-diffuse bases affordable — so it is the
natural thing to reach for when the metal/anion work (roadmap steps 3–4) makes diffuse bases central.

---

# Parked/background (unchanged from GPWPlan.md)
0.5(e) runtime folded into item 3.  0i analytic V_local LONG (fold the core charge into PW_Hartree's G-space
solve; after it, whether the completion rung still pays at BallOnly+C=2 is a measurable rung-gate question).
§2 low-q multi-species Si/NaF/CsI cross-validation.  §3 CP2K reference library growth.  §5 cleanups (Vxc-fit
ISP ctor with GGA; cMesh unification; DRY PP adapters; periodic external PP).

---

# valgen `--auto` : rule-based valence-window refinement (PLAN, not built — 2026-07-28)

The diagnostics valgen now emits (basis-usage heat map + per-shell advice + `--floor` gap + per-irrep cond +
pivoted-Cholesky `removed`) ARE a deterministic tuning algorithm — no ML needed.  `--auto` would close the loop:
iterate the even-tempered knobs (emin/emax/N, equivalently β) until the energy sits near the PP floor AND the
pop histogram is clean.  Captured so we don't lose the hand-tuning expertise built up in this session.

**Objective (per l-shell; s and p decouple at the atom level):** minimise `gap_mHa = E_window − E_floor` subject
to a CLEAN histogram (no alternating ±, diffuse end well-populated, cond below ~1e3).  Stop when `gap < tol`
(≈2 mHa) and clean, or at a max-iteration cap.

**Seed guess:** emin from the seed-density ⟨r⟩ (diffuse extent), emax from the PP `r_loc` (α_max ≈ 1/(2 r_loc²)),
a default N≈5–6.  (Or accept a user seed window and only refine.)

**Loop = the advice rules turned into ADJUSTMENTS (each step logs which rule fired):**
- diffuse end is the peak (pop₀ ≈ pmax, at boundary)  → LOWER emin (×0.5)         [density spills past]
- diffuse end barely populated (pop₀ ≪ pmax)          → RAISE emin (×2)           [wasted diffuse fn]
- pops ALTERNATE ±, or cond > 1e3, or β < 2           → raise β: N−1 (or widen range) [overcomplete]
- gap still large & sharp end still climbing          → extend emax (or N+1)       [incomplete inner region]
- gap already small & a sharp low-pop shape tail      → trim emax/N                [leaner, costs a few mHa]

**Guards:** clamp emin ≥ ~0.02 (the LDA-XC-on-a-diffuse-density NaN floor — see step-1 notes); on an SCF throw,
RAISE emin and retry (the diffuse XC-NaN or a conditioning failure).  Cap iterations; never oscillate emin (damp).

**Why deterministic, not an optimiser over a loss surface:** the floor is the stopping test; β/cond are the
conditioning guards; each move is one of the 5 rules above — explainable, reproducible, cheap (one atom SCF +
one pool-floor SCF per step).  Bisection-like on emin, integer nudges on N.

**Honest caveat (the pin):** the atom is only a PROXY.  The real objective is grid-convergence of ρ / a bulk
orbital-coefficient heat map IN THE SOLID (GPWPlan §1, NaF pin).  `--auto` produces a good STARTING window fast;
the last word is a bulk vetting pass once the basis is in GPW.  So build `--auto` for ergonomics, but don't let
its atomic `gap_mHa` become the objective function (that's the oracle-matching trap we already rejected).
