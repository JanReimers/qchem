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

**BECKE XC GRID — BUILT + GATED + PERF-CLOSED (2026-07-30/31)** (commits e601cada, 2b08f98c, 70b32905,
c800ec48, 87acb699, 719d6b6a, 33a94ce4, 2fd152c8; remaining increments in their own section below):
- **Periodic Becke mesh** — `qcMesh::UnitCellKind{Uniform,Becke}` in `MeshParams`; `UnitCell::
  CreateIntegrationMesh` Becke branch: per-atom radial×angular grids, competitor images gathered in
  expanding Chebyshev shells until both tails are provably < ε by the s(μ) collinear bounds (magnitude
  screen only, NO radius — the pin holds; ε=1e-8 by measurement); points wrapped home; far diffuse tails
  drop free.  Acceptance `LatticeMesh.Becke*`.  **ANGULAR AUDIT**: the 40-year-old hand-coded Lebedev
  tables are ALL CORRECT to their claimed degrees (permanent `Mesh_Angular.*` gates); GaussLegendre is
  machine-exact at ANY L; EulerMaclaren has NO algebraic degree.
- **`PW_XC_Becke` + `BeckeXC_Engine`** (pair-shared mesh + geometry-fixed Φ tables keyed `BasisSetID` +
  per-serial ρ; `tDM_CD::DM_RhoAtPoints` = the FIELD dual of `DM_ContractBlocks` — the density GEMMs the
  caller's tables against its private D; `tChargeDensity::EvalBatch` factorised-phase batched inverse FT
  for the mixed/Kerker iterations where the Fock build sees the ρ̃ map, not a DM): per-iteration cost
  4.8 → 0.37 s; NaF Becke ≈37 s vs uniform 34 s, energies bit-identical.  ROUTE SEAM: `Ham_PW_DFT`
  optional `xcMesh` ctor param → `GpwOptions.xcMesh`, one-token `UnitCellKind` switch in
  `DISABLED_NaFRocksaltGamma`; every run announces `[XC quadrature]` / `[Becke grid]` / `[uniform grid]`
  + `grids.*` report entries.  Resolution knobs `GPW_BECKE_L/NR/ALPHA/ANG`.
- **GATES**: Si Γ (enabled): dE_xc=1.1e-4 Ha, max|ΔV_xc|=3.5e-4 vs a MATRIX-GRADE Ecut=60 uniform
  reference (lesson: SCF-grade Ecut reads as a false Becke error).  NaF SR2 sharp-F stress (DISABLED,
  long): the Becke matrix is INTERNALLY converged < 1e-4 including the F core (B40-vs-B80 = 7.1e-5)
  while the uniform raster's worst F-core element wanders 0.16 with raster config, and E_xc drifts
  TOWARD Becke as the raster refines — the grid's reason to exist, measured.
- **ANGULAR-ORDER LADDER (user's site-symmetry conjecture — confirmed)**: T_d-sparse content ⇒
  converged at GL-11 (72 dirs, sub-mHa), breakdown ≤ degree 7.  Lebedev at the same degree is better on
  V_xc ELEMENTS (the O_h-orbit cancellation is real) but 5-10× worse on ρ-weighted integrals — its
  (±1±1±1)/√3 orbit sits exactly ON the diamond bond axes ⇒ site-adapted grids must AVOID special
  directions.  Defaults: GL-17 safe / GL-11 fast.
- **PIN — HartreeOnly routing FALSIFIED as a default**: dropping the ⅔·α_max field-sharpness floor when
  the XC is off-raster destabilises the diffuse NaF at EVERY reduced β (β=0 → +904 Ha) — **the floor
  protects the DENSITY's low-G integrity (under-resolved diffuse-pair FFT fold-back corrupts exactly the
  charge-transfer mode), NOT only V_xc** (`b3dad5be`'s attribution was partially incidental).  Ships
  default-inert: the `GPWParams::RasterFields` seam + `GPW_ROUTING`/`GPW_RELFIELDSHARP` calibration valves.
- **DIAGNOSTICS**: `[lattice sums]` economy line (α_min/max, the ε's with env names, worst-pair reach →
  `CellsInSphere` count — the numbers that jump with diffuse functions) + the `qchem.Reporting` TIMING
  LEDGER (`Timed` RAII, EXCLUSIVE times — disjoint buckets sum to wall — sorted by cost at run end).
  First measurement NAMED the diffuse-NaF 285 s: local-PP LONG grid sweep 180 s + stream build 65 s,
  both driven by the 791-cell per-pair offset enumeration (explains the κ-null + the small screen-ε
  wins); local-PP SHORT analytic = 0.37 s over the SAME enumeration (~500× per-term) ⇒ 0i's G-space
  fold is the structural kill of the biggest bucket.
- **OPENMP FOR REAL**: clang's `-fopenmp=libgomp` is LINK-COMPAT ONLY — no codegen; every "OpenMP" run
  2026-07-19→07-31 was silently serial (the old load-imbalance / 0-speedup verdicts are VOID).  C++20
  modules force OpenMP GLOBAL (the language mode is baked into each PCM; mixed imports refuse both ways)
  ⇒ `CMAKE_CXX_FLAGS` + global define/link, with the apt-libomp bridge documented beside
  `find_package(OpenMP)`.  First REAL threaded diffuse-NaF run: 4:09 → 1:08 — which exposed the
  pair-loop thread-safety bug (resolved same day; next bullet).
- **THREAD SAFETY + PARALLEL STREAM BUILD (2026-07-31→08-01, commits 1a98a56e, e9bb964e)**: the
  intermittent threaded abort was BLAZE, not our code — global `-fopenmp` defines `_OPENMP` → Blaze SMP
  auto-activates → every `smpAssign` (even tiny vector copies) opens `BLAZE_PARALLEL_SECTION`, whose
  guard is a GLOBAL non-atomic bool → two pair workers race → throw → terminate.  **PIN: any project
  compiling Blaze under global -fopenmp needs `BLAZE_USE_SHARED_MEMORY_PARALLELIZATION=0`.**  Plus throw
  CONTAINMENT in both pair loops (worker exception rethrown serially, never terminate); TSan+Archer
  certificate (`build/TSan`) = zero project races; threaded ≡ serial trajectory to ~1e-10 over 45 iters.
  The stream build is now a parallel TWO-PASS build (parallel count → serial budget walk on the counts →
  parallel build directly in final tier form): 64.9 → 20.1 s, readouts + replay BIT-IDENTICAL to serial,
  fp32 demotion transient eliminated.  Diffuse-NaF time-to-first-iteration ≈ 25 s threaded (was ~250 s
  serial-everything).
- **0i V_loc-long CUSTOM G-BALL — DONE (2026-08-01, commit 6e9b03b4; the user's design)**: a STATIC
  EXTERNAL field is entitled to its own cutoff — vlocEcut ~ 2αmax+β sharpened to the exact harmonic
  per-pair rule req(p)=2ln(1/ε)·pβ/(p+β) (`StaticFieldPairLevels` + explicit `pairLevels` on
  `IntegratePotential`), one custom top level (NaF 214 Ha vs the 160 rung).  NaF custom-ball vs κ-sweep
  = **5.5 μHa**; setup bucket **180 s serial / 34 s threaded → 1.7 s**.  `GPW_LONG_SWEEP=1` = A/B
  instrument + non-Gaussian-model fallback; `GPW_VLOC_EPS` = self-convergence instrument.
  `AlFCCMetalIBZExact` re-anchored −2.116812 → −2.1169707 (full-mesh and IBZ shift TOGETHER; folding
  stays exact).  The load-bearing insight (measured via the retracted smooth-fold experiment, 4.6 mHa):
  per-level band-limits are exact for per-iteration KS fields (adjoint-consistent with collocation) but
  a plain truncation error for external static fields.  The σ-split alternative is parked in git history.
- **Na METAL FIRST LIGHT (2026-07-28, gate `NaFCCMetalGlobalMu`)**: shifted-MP 2×2×2 + global μ +
  smearing on FCC Na (`VALENCE_LOWQ_SR2` at the real a=10 density, cond≈38) = the textbook smeared Fermi
  surface (μ mid-band, fractional n_k=0.67 on-surface, Σ_k w_k n_k = 1.0000 exact, 1.3 s); annealing
  composes (kT-dependence is real for a true Fermi surface; T→0 recovered).  Caveats live in roadmap 4.

---

#roadmap (agreed 2026-07-26 — toward an honest metal + Fermi-smearing test)

*(Items 1–3 + IBZ/k-star are DONE — summarised in the DONE section above.  Remaining forward work:)*

**4. MULTI-K Na — the honest metal test, DONE PROPERLY** — first light is in (`NaFCCMetalGlobalMu`, above);
the remaining pass is the CONVERGED metal: a denser mesh (needs IBZ — now available), a metallic Na basis via
valgen (the SR2 gate validates the machinery, not a cohesive energy), and a physical did-E-move / bulk anchor.
Optional BCC once a `BCCUnitCell` exists (only FCC today).

**5. NON-SYMMORPHIC IBZ support — DONE** (moved to DONE-adjacent summary in the items 1–3+IBZ block above;
gates `SiDiamondIBZ_NonSymmorphic` −7.77847 vs full mesh −7.77846 + `GMapUT` structure-factor pins; the τ-phase
machinery = `SymmetrizeGMap` G-space phase + `SymmetrizeRaster` exact FFT fractional shift).
**REMAINING here: the spin-polarized global μ (both spin channels, one μ).**

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

# Becke XC grid — remaining increments (the build/gate/perf record is in DONE above)

*(Original motivation, design, and gate spec served and condensed into DONE; full narrative in git history.
The framing that stands: the Becke grid gives the XC — the one pointwise-nonlinear, sharp-at-core term —
its near-ideal grid across every basis set, diffuse included, while Hartree stays on the G-space Poisson.)*

*(The 2026-07-31→08-01 thread-safety + parallel-stream-build + 0i custom-G-ball wave is condensed into
DONE above.  The Becke-ON anchor failure surfaced during it — a −50 Ha ghost dive, identical serial and
threaded — was RESOLVED by the user's recipe retune (relax 0.45→0.25, delayed DIIS): the test now
converges and passes the −24.4304 anchor (measured −24.4317, 2026-08-01).  With thread-safety and 0i
both settled, the "Becke as default" decision below is now UNBLOCKED.)*

- **SPACE-GROUP STREAM/COLLOCATION REDUCTION (PLAN, 2026-08-01 — feeds the symmetry code review; the
  SpaceGroup SRP encapsulation is the home for all of it).**  Today the streams use ONLY the Hermitian
  j≥i fold (×2); the space-group redundancy is untouched.  Where it lives (for a 1-atom-per-species
  primitive cell the ops FIX the atom assignment, so NOT mainly in pairs): (1) **cross-cell offsets** —
  a (pair, R) stream maps under a site-symmetry op to (pair′, W·R) + a raster-index permutation; the
  diffuse pairs' ~100s-of-offsets lists (most of the ~513M cached points) collapse by the pair's
  site-symmetry order (up to 48 in O_h) — the 5–20× lever; (2) **polarization components** — Cartesian
  monomials under cubic ops are signed axis permutations: pair→±pair one-to-one, no shell mixing.
  TWO ROUTES: (a) orbit-expansion replay (store irreducible streams, replay members through per-op
  raster permutations — exact, no physics constraint, cache-hostile gather); (b) collocate IRREDUCIBLE
  representatives with orbit weights + star-average ρ once per iteration (the IBZ `SymmetrizeGMap`/
  `SymmetrizeRaster` machinery EXISTS), mirror on the integrate-back (representatives + representation
  transform of h) — also cuts the PER-ITERATION scatter/gather by the same factor.  Route (b)'s string:
  it is exact only for D symmetric under the group — it IMPOSES the symmetry.
  **THE IMPOSED-GROUP POLICY (the spontaneous-symmetry-breaking question — LiMn2O4 charge/spin
  ordering):** the group used for ANY reduction (IBZ, star-average, streams) must be a RUN-LEVEL POLICY
  — "the group to impose", a subgroup of the geometric space group — NOT hard-wired to the lattice's
  full group.  Both ordering workflows need exactly this machinery: (1) SSB-search runs impose a
  subgroup (or trivial group) in an ordering-commensurate supercell, with a symmetry-broken SEED
  (SeedStrategy) and +U (delocalization error suppresses disproportionation at LDA/GGA); (2) the
  cluster-expansion/MC route ([[battery north-star]]) imposes each CANDIDATE ordering's stabilizer
  subgroup and compares energies — same code, enumerated groups.  Spin orderings want the MAGNETIC
  (Shubnikov) group axis: ops composed with spin flip / time reversal.  Supporting pieces: converge-
  symmetric-then-RELEASE staged workflow (density-continuation machinery exists); an ORDER-PARAMETER
  DIAGNOSTIC (monitor D's symmetry under the FULL group while imposing the subgroup — a run REPORTS the
  ordering it found, never silent).  SpaceGroup's SRP surface then owns: group/subgroup computation +
  selection, (pair, offset) orbit generation, monomial representation signs, raster-commensurability
  checks (below), star-averaging.
  **RASTER COMMENSURABILITY (pointwise ops only):** a POINTWISE index permutation exists iff every op
  maps the integer torus to itself: axes mixed by W need EQUAL N, and each translation component needs
  N_i·τ_i ∈ ℤ — diamond's τ=(¼,¼,¼) ⇒ 4 | N per axis ON EVERY LADDER LEVEL.  The 5-smooth FFT padding
  does NOT guarantee this (45, 27, 15, 3... are 5-smooth; NaF's measured ladder had N=45 and N=3
  levels), so enabling stream-level symmetry requires the raster menu to round up to 5-smooth MULTIPLES
  of lcm(τ denominators) (45→48 etc. — modest).  Diamond Si already FORCES this path in the suite
  (`SiDiamondIBZ_NonSymmorphic`: AutoGrid yields odd N, 5-smooth padding keeps them odd — 4∤N is the
  generic case); a trigonal/hexagonal screw material (quartz P3₁21, τ=c/3; P6₁-family, τ=c/6) is the
  extension gate for non-power-of-2 denominators once the feature lands.  The EXISTING IBZ
  symmetrization is exempt: its FFT fractional shift e^{iG·τ} is exact for the band-limited interpolant
  at ANY τ — the constraint is specifically for sparse compact streams, which cannot take an FFT shift
  without densifying.  **k≠Γ:** the offset lists also carry Bloch phases e^{ik·R_n}, and ops act on k as
  well — stream reduction at general k must stay consistent with the IBZ k-star folding (the same star
  bookkeeping; one more thing `SpaceGroup` should own rather than the evaluator).
- **Coarse-end routing calibration** (`kMinLevelN=3` "very coarse for a diffuse lobe" — the standing
  suspect after the HartreeOnly falsification; valves `GPW_ROUTING`/`GPW_RELFIELDSHARP` are in place).
- **Becke as the DEFAULT — ✅ FLIPPED 2026-08-01.**  `GpwOptions.xcMesh.cellKind` defaults to the new
  `UnitCellKind::Auto`, resolved by the driver (`ResolveXCMesh`): **Becke** (the calibrated
  `BeckeXCParams()` recipe) — EXCEPT under `reduceBZ`, where Auto falls back to the uniform raster WITH
  a console note (the Becke-route density star-average is unverified — the "Becke+IBZ" item, owned by
  the symmetry review; explicit Becke+reduceBZ is honored, which is how the verification run will be
  done).  Non-GPW consumers treat Auto as the historical Uniform (zero blast radius).  Gate sweep: 11/12
  GPW_SCF pass on the flip (metals, smearing, both IBZ tests via the carve-out); two gates pin Uniform
  explicitly BY PURPOSE (`BeckeXCMatchesUniformXC`'s uniform arm; the box gate below).
  **SURFACED by the flip — Becke × DEGENERATE OPEN SHELL (open):** the aufbau pseudo-atom-in-box gate
  (half-filled 3p, density rotates freely in the degenerate shell) OSCILLATES ~Ha-scale under Becke
  where uniform was grid-stable: V_xc is pointwise-nonlinear, so an anisotropic ρ's quadrature error on
  the FIXED-AXIS angular grid rotates with the density — the energy-neutral zero mode becomes an
  energy-oscillating mode.  The SMEARED sibling gate converges fine under Becke (fractional occupation
  restores symmetric ρ) — consistent with the mechanism.  Fixes to weigh: smear/GDM for degenerate
  open shells (the existing recipe), or the site-group-adapted angular grid below (an orientation-robust
  quadrature).  Recorded here so Becke's open-shell/atom-in-box claim waits on it.
- **Spin-native** (`PW_XC`/`PW_XC_Becke` are unpol today; spin-native is the formulation pin).  At
  minimum settle the spin-native INTERFACE shape before the symmetry review: magnetic (Shubnikov) group
  support must not be designed around the unpol special case.  **TEST-MATERIAL LADDER (user, 2026-08-01
  — the interface needs a magnetic material to test against):** (a) **Na pseudo-atom in a box, doublet**
  (S=½, the EXISTING Na q1 PP — the minimal end-to-end spin-native GPW run: two channels through
  collocation + XC, moment=1); (b) **O₂ in a box, triplet** (O q6 GTH) — CROSS-ANCHORED against the
  molecular facade's spin-native triplet-O₂ gate, the spin sibling of `SiPseudoAtomInBoxMatchesFinite`;
  (c) **MnO rocksalt AFM-II** (2-f.u. cell, moments along [111]) — the first REAL d-electron magnet:
  exercises the magnetic-group/imposed-subgroup policy (AFM breaks the full group), wants +U for honesty,
  and is the direct rehearsal for LiMn2O4 charge/spin ordering; (d) LiMn2O4 itself.  (a)+(b) need no new
  physics beyond spin-native XC — they are the interface-shape tests; (c) is where the symmetry-review
  machinery (Shubnikov ops, imposed group, order-parameter diagnostic) gets its first real workout.
- **Site-group-adapted angular grid** — with the measured design pin: orbits must AVOID special (bond)
  directions (the Lebedev cube-corner lesson); the exact required degree is deducible from the site group.
  **SCOPE WIDENED (user, 2026-08-01): this is another SYMMETRY-REVIEW item, and it SUBSUMES the GL-29 vs
  GL-17 default question** — the site point group should determine BOTH how many and WHICH angular
  directions each atom needs (a group-invariant orbit set of the required degree), rather than a global
  hand-set L.  Lebedev's octahedral-orbit grids can beat GL per point at O_h-compatible sites (the
  measured V_xc-element win) PROVIDED the orbits avoid the site's special directions; GL stays the
  generic-site fallback.  Also the candidate cure for the Becke × degenerate-open-shell oscillation
  above (an orientation-robust / site-symmetric quadrature removes the rotating-error channel).  Until
  then GL-29 stays the default (anchors verified against it); GL-17 is the interim suite-time lever.
- **Per-element radial scaling** (one mhl_alpha serves all species today; NaF's F core vs Na wants per-Z).
- **Becke+IBZ**: the real-space star-average of ρ at mesh points is untested on this route (k-fold weights
  ride in D; non-Γ star symmetrisation needs a check or a mesh-point star-average).
- **Mixer Becke-mesh shadow** (the raw-through-dynamics analogue of 0.5(f2)) if ball-Gibbs ever bites at
  the mesh points — on mixed iterations the Becke ρ feed is the ρ̃-BALL field (FourierMixCD's map).

**Beyond LDA — the functional-ladder win (user, 2026-07-26; the strategic payoff, still ahead).**  A
semi-local XC potential is a pointwise-nonlinear function of ρ AND its derivatives — ∇ρ (GGA), τ/∇²ρ
(meta-GGA) — and on a Gaussian basis those are ALL ANALYTIC at each mesh point.  On the uniform grid the
ladder demands ever-finer grids (`relCutoff` escalation); on the Becke grid the derivatives are exact per
point and the mesh is already dense at the cores, so **LDA → GGA → meta-GGA costs the same**.  Hybrids:
exact exchange stays the analytic ERI machinery; the semi-local portion (most of a PBE0/HSE) is exactly
what Becke makes affordable.  The Becke grid is the enabler for the whole functional upgrade, not only a
diffuse-basis fix.

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
