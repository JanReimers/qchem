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
  `BeckeXCParams()` recipe) — originally EXCEPT under `reduceBZ` (the Becke-route star-average was
  unverified), a carve-out **RETIRED 2026-08-02** with SymmetryUpgradePlan §7 step 5 (W2c): Auto is
  Becke UNCONDITIONALLY; an imposed run (`imposeSymmetry`, the renamed `reduceBZ`) builds the
  mixed-rule site-adapted invariant Becke mesh (~2× points at the production recipe) and star-averages
  ρ on it (Al IBZ re-pinned route-matched −2.1174805).  Non-GPW consumers treat Auto as the historical
  Uniform (zero blast radius).  Two gates pin Uniform explicitly BY PURPOSE
  (`BeckeXCMatchesUniformXC`'s uniform arm; the box gate below).
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

# Stream-cache RAM/CPU trade — promote the knob + the fast-recompute campaign (USER, 2026-08-13)

Born from the spherical MnO arm-2 runs (doc/SphericalLatticePlan.md): the CP2K-span basis (n=152
Cartesian inner, diffuse-heavy) drives the stream cache to ~9 GB where the SR-trimmed runs took ~6 —
the diffuse pairs' collocation boxes, not the function count, are the cost.  CP2K's design point is
the OPPOSITE trade: no stored streams at all — every pair product recomputed per iteration by
ferociously optimized compact-box grid kernels, with DBCSR block sparsity + EPS_PGF_ORB screening
making the recompute cheap.  Both are exactly variational (adjointness is a property of the
collocate/integrate PAIR, not of caching).  Two items:

1. **Promote the budget to a first-class user knob.**  It EXISTS today as env-only
   (`GPW_STREAM_BUDGET_PTS` / `_F32`, two tiers, overflow → on-the-fly; PG_Cart_MnD/Evaluator.C)
   with compile-time defaults sized to ~8.6 GB peak.  Surface it through GPWParams (and the
   facade), document the tiers (fp64 bit-identical / fp32 ~6e-8-relative / recompute), and report
   the chosen point + coverage in the run header so a user SEES the trade their box made.

2. **The fast-recompute campaign.**  The fallback path is today's SLOW analytic sweep (~9 min/
   iteration on the magnetic MnO cell when the budget ran out, 2026-08-05) — it must become a
   CP2K-class kernel so low-RAM boxes can run diffuse spans at speed: per-pair compact-box
   vectorized Gaussian-product evaluation (the box enumeration already exists — it BUILDS the
   streams today), pair/offset magnitude screening tiers reused from the stream build, OMP over
   pairs (the GPW_OMP_THREADS seam), and — once an engine answers spherical natively (the
   SphericalLatticePlan kernel-hoist) — pair spaces of 5 not 6 components per d shell, which cuts
   BOTH modes ~20% on d-heavy bases.  With (2) fast enough, (1)'s default could even flip toward
   recompute on small-RAM boxes.  **HEAD START (user, 2026-08-13): the per-pair ALGEBRA is already
   cached** — `struct Ω : Cacheable2` (PG_Cart_MnD/GaussianRF.C), content-addressed and SHARED with
   the spherical evaluator (same primitives ⇒ same registry indices) — so the recompute kernel
   replays Ω and pays only the vectorized point loop over the compact box; the pair algebra is
   never rebuilt.
   MEASURED MOTIVATION (runs 49 first attempt): the CP2K-span MnO cell (n=152 Cartesian inner,
   diffuse-heavy) exceeded the budgets so far that iteration 1 did not COMPLETE in 80 min on the
   on-the-fly path — the overnight retry runs at GPW_SCREEN_EPS=GPW_DENSITY_EPS=1e-8 (the
   sanctioned ε-truncation; sub-mHa vs the 45 mHa signal).

# ★ THE RUNTIME GAP, MEASURED — and where it actually lives (2026-08-15, doc/OpenWork.md item 1)

The charter above said "fast-recompute campaign".  **The ledger says otherwise for the runs we
actually pay for.**  First order of business was making the cost VISIBLE: `GPW_REPORT=1` forces the
run report (hence the timing ledger) on for ANY driver, and new buckets name the per-iteration work
that had none — `scf: collocate density (pair scatter)`, `scf: integrate-back (pair gather)`,
`scf: XC-mesh quadrature H_xc`, plus a split of the rho sampling into the Φ-GEMM path and the
MATRIX-FREE path.  Measured on the MnO AFM-II magnetic cell (Γ, n=122, 48k-point Becke mesh,
`GPW_OMP_THREADS=8`, `GPW_MNO_NMAX=2`; 171 s total):

| bucket | s | what it is |
|---|---|---|
| setup: XC-mesh Phi tables | 55.6 | the Bloch basis at every mesh point — SERIAL |
| setup: collocation stream build | 27.2 | (already threaded) |
| scf: XC-mesh rho sampling | 23.8 | Φ GEMM + matrix-free sampling — SERIAL |
| scf: iterate residue | 21.8 | mostly the H_xc quadrature GEMMs — SERIAL |
| setup: becke mesh build | 16.6 | (already threaded) |
| scf: collocate + integrate (THE PAIR LOOPS) | **7.6** | the charter's target |

**The pair loops are 4% of the run once the streams are cached.**  Every dominant bucket is dense
O(npts·n) / O(npts·n²) work on the ATOM-CENTRED XC MESH, and every one of them was serial while the
pair loops — the one place with OMP — were not.  The fast-recompute kernel remains the right answer
for the OVER-BUDGET case (item 2 above, unchanged: the CP2K-span cell that never finished an
iteration), but it is NOT what makes a fitting-in-budget magnetic cell slow.

**What landed (2026-08-15):**
- `qchem.Parallel::WorkerThreads()` — the ONE opt-in worker count (`GPW_OMP_THREADS`, serial by
  default), lifted out of `NR_Evaluator::PairThreads` so every hot site reads the same number.
- OMP at four new sites, each partitioned by OUTPUT so the result is **bit-identical at any thread
  count**: the Φ table rows (`XC_GridEngine::Phi`), the rho-sampling point blocks
  (`IrrepCD::DM_RhoAtPoints`), the H_xc quadrature's output column blocks (`XC_GridEngine::Matrix`),
  and the batched inverse FT over mesh points (`FourierMixCD::operator()(rvec3vec_t)` — the ρ̃-mixed
  density the Fock build sees on EVERY Kerker/Pulay iteration).
- `XC_GridEngine::Matrix` also went TRIANGULAR (H is Hermitian, `chmat_t` stores i≤j): half the GEMM
  flops.  The retired full-M build's `½(M+M†)` only symmetrised roundoff (w, v real ⇒ M Hermitian).
- SHELL HOIST in `PG_Cart::IrrepBasisSet::operator()`: a shell's components share one contracted
  radial (PGData stores them consecutively), so the exps are evaluated once per SHELL, not once per
  component — 6× fewer on a d shell, on THE pointwise sweep behind every Φ table.
- `SeedCD` gained a batched `operator()(rvec3vec_t)` (threaded; and the cell metric \f$(A^TA)^{-1}\f$
  lifted out of the point loop, where `op()` recomputed a 3×3 inverse PER POINT).

**Result, same cell, same recipe:** 171 s → 101 s (2-iteration benchmark), physics unmoved
(m_stag 0.3171 both).  Φ 55.6 → 7.0 s (7.9×); iterate residue 21.8 → 6.9; seed+ortho 11.2 → 4.2;
matrix-free sampling 26.9 → 8.3 (measured at 4 iterations).  Per-iteration ≈ 30 s → ≈ 9 s, so a
24-iteration production run of this cell goes from ~13 min to ~5 min — and the same sites carry
EVERY atom-centred-XC run (Becke is the GPW default; this is not an MnO-only win).

**What is left, in order (the ledger after the fix):** the collocation stream build (~30 s setup,
already threaded — the recompute campaign's real target), the Becke mesh build (~18 s; `BeckeCutoff`
alone is 14% of all CPU), the H_xc quadrature (still the largest per-iteration bucket after
halving), and the pair scatter/gather (~3.5 s/iteration).  Two structural levers NOT yet taken:
- **SCREEN the Φ table.**  It is stored dense (npts × n) but the true object is SPARSE — a function
  is numerically zero over most of an atom-centred mesh.  Batching the mesh (per atom / per radial
  shell) with a per-batch significant-function list turns every Φ-shaped cost (build, rho GEMM, H_xc
  GEMM) into O(npts·n_sig²) — this is what gives Gaussian codes their O(N) XC, and the win GROWS
  with cell size (MnO's 4 atoms understate it).  THE next increment on this path.
- **A real-Γ fast path.**  Γ-only runs (all the MnO/NaF/Si production) carry Φ and D as complex with
  exactly zero imaginary part: 4× the flops for nothing.

## Round 2 (same day): the BLAS pin, and Blaze's dispatch rule

**The pin came back.**  `openblas_set_num_threads(1)` had been lost in the 2026-08 machine migration
(the `#include <cblas.h>` in the test mains outlived the call).  It now lives ONCE, as
`qchem::PinBlasToOneThread()` next to `WorkerThreads()` in `qchem.Parallel`, called from all four
`main()`s (ITMain, the two shared UT mains, scfrun).  `qcCommon` links `openblas` by NAME on purpose:
a BLAS swap must fail LOUDLY at link, which is exactly what did not happen last time.  The reason it
is a correctness knob and not a tuning one: OpenBLAS auto-sizes its pool from machine LOAD, so the
reduction order inside a GEMM becomes load-dependent and the last ULP moves between runs of the same
binary — on top of oversubscribing our own OpenMP regions.

**`QCHEM_BLAZE_BLAS` (default ON, 2026-08-15).**  Blaze's own kernels build at the SSE2 baseline here
(`BUILD_MARCH_NATIVE=OFF`).  Measured on the XC-mesh shape (48033×122 complex, one thread):

| | GFlop/s |
|---|---|
| blaze native | 2.55 |
| OpenBLAS `zgemm` | **41.7** |

**★ THE DISPATCH RULE (the finding worth keeping): Blaze hands a product to BLAS only when the
DESTINATION is a plain matrix.**  Assign into a `submatrix` view and it silently falls back to its own
kernel — 1.87 vs 34.1 GFlop/s on the same operands.  So the hand-blocked, hand-threaded, triangular
form that round 1 added to `XC_GridEngine::Matrix` (halving flops, spreading over 8 threads) was
LOSING 13× to save 2× the moment BLAS mode came on.  Under BLAS mode the quadrature is now ONE
whole-matrix product, no blocking, no OpenMP — and it beats the 8-thread hand-rolled version.  The
same rule bit the ρ sampling from the other side: a `HermitianMatrix` ADAPTOR operand is not
BLAS-dispatchable either, so `DM_RhoAtPoints` copies D into a plain matrix (an n×n copy against an
npts·n² GEMM) before multiplying.  Both shapes are kept: `#ifdef BLAZE_BLAS_MODE` picks the
whole-matrix product, else the blocked threaded triangular one.

**Measured (MnO Γ, n=122, 48k mesh, 8 threads):** rho sampling 4.17 → **0.18 s**; H_xc quadrature
11.8 s/4 iters → 4.1 s/2 iters (2.95 → 2.07 s per iteration).  Per-iteration total ≈ 8.7 → ≈ 7.6 s.
`ctest -j8` 716/716 green, and the sweep itself dropped 615 → 540 s.

After this the per-iteration ledger is: pair scatter+gather ≈ 3.8 s, H_xc ≈ 2.1 s, matrix-free ρ̃
sampling ≈ 1.35 s, rho GEMM ≈ 0.09 s.  Setup is dominated by the stream build (~28 s) and the Becke
mesh build (~16 s).  The next levers are unchanged: Φ-table screening, the real-Γ/TRIM path (Γ AND
every k with 2k ≡ 0 — so a Γ-centred 2×2×2 mesh is entirely real), and the fast-recompute kernel.

## Round 3 (2026-08-19): the collocation streams — and why "make them REAL" was the wrong diagnosis

The real-TRIM flip (doc/RealComplexPlan.md) left one named increment behind it: *"the complex-internal
evaluator collocation streams, ~150 s of the MnO profile."*  The arithmetic checked out — 87.5 s
(collocate) + 41.1 s (integrate-back) + 23.5 s (stream build) = 152 s of the 402 s free-run ledger.
**The DIAGNOSIS did not.**

**THE MEASUREMENT** (`perf record --call-graph=fp -F199`, MnO AFM-II Γ, n=122, 8 threads,
`GPW_MNO_NMAX=3`; annotation of `CollocateDensity`'s pair-scatter lambda, which carries 10.7% of whole-run
self time):

| instruction | share of the lambda | what it is |
|---|---|---|
| `movss (%rsi,%r8,4)` | 25.4% | the fp32 VALUE stream |
| `movl  (%rax,%r8,4)` | 23.7% | the per-point INDEX stream |
| `addsd (%rcx,%r9,8)` + store | 8.5% | the scattered destination RMW |
| `mulsd %xmm2,%xmm0` | 0.5% | the actual multiply |
| `callq *0x18(%rdi)` (the Bloch-phase `std::function`) | — | did not register |
| `__muldc3` (the complex \f$D_{ij}\overline{e^{ikR}}\f$ contraction) | — | NaN path only |

**The replay is DRAM-BANDWIDTH bound, and the scalar type of the contraction is free.**  The phase and
the complex \f$D\f$ contract are PER (pair, offset) — amortised over a run of ~20–40 grid points each —
while the streams themselves were always real (`std::vector<double>`/`<float>`).  Re-typing the
`LatticeSum1E` collocation face to `double` would have bought approximately nothing here.  What the
bucket wanted was **fewer bytes per point**, and the same annotation says where they were: HALF the
traffic was the per-point raster index.

**WHAT LANDED INSTEAD — run-length stream geometry (`PairOffsetStream`).**  `ForPairBox` walks the box's
innermost axis in unit steps, so its wrapped indices arrive in long CONTIGUOUS runs (one per \f$(dx,dy)\f$
column, split only by the modulo wrap or the value screen).  A stream now stores `runBase`/`runLen` —
(first index, length) per run — and the values alone; both replays (scatter and gather) walk `val`
sequentially and the grid contiguously.  Same points, same order, same screens ⇒ **bit-identical replay**,
and the collocate/integrate adjoint stays machine-exact.  The readout (`[stream cache]` /
`grids.stream`) now prints `runs` and `meanRun` so the encoding's payoff is visible — and it must be
read, because the win scales with run length: a run head costs 8 B, so the encoding pays for
`meanRun > 2` and pays WELL above ~8.  **Measured `meanRun` = 10.16 on MnO** (647 M points in 63.7 M
runs) and 4.19 on the small Si `AnalyticCollocationCrystalChargeConservation` grid — so the tiers go
12 → 8.8 B/pt (fp64) and 8 → 4.8 B/pt (fp32) on MnO, i.e. **stream RAM 5.78 → 3.70 GB, a 2.1 GB cut
that lands directly on doc/OpenWork.md item 2 (the RAM gap)** at unchanged coverage.

**And the OTHER complex-bound bucket, which was genuinely complex: `FourierMixCD::operator()(rvec3vec_t)`**
— the ρ̃-mixed density sampled at every Becke mesh point on every Kerker/Pulay iteration (32k points ×
24k G × 2 channels), the single largest per-iteration bucket in the free-run ledger at 87.6 s.  Two exact
regroupings of the same series: (1) ρ is REAL, so the ±G terms pair as
\f$\mathrm{Re}[(c(G)+\overline{c(-G)})e^{iG\cdot r}]\f$ — fold onto a canonical half-space ONCE per batch
and the point loop runs over half the terms (an identity, not a symmetry assumption: an unpaired G lands
alone in the folded weight); (2) the map's key order runs z fastest, so the \f$(x,y)\f$ phase multiplies
the run's inner sum ONCE instead of once per G.  Together 2 complex multiplies per G over the full {G}
become 1 per G over half of it.

**MEASURED** (MnO AFM-II Γ, 3 iterations, 8 threads — the same run three times):

| bucket | before | after | |
|---|---|---|---|
| `scf: collocate density (pair scatter)` | 7.79 | **5.40** | 1.44× |
| `scf: integrate-back (pair gather)` | 3.40 | **2.42** | 1.41× |
| `scf: XC-mesh rho sampling (matrix-free density)` | 2.58 | **1.08** | 2.39× |
| per-iteration SCF total | 5.52 | **3.88** | **1.42×** |

Physics unmoved to every printed digit (`Efinal=-61.3944791877`, `lastΔρ=0.0012609100`, `m_stag=0.5977`,
`m_net=0.03998`, and the whole energy decomposition) — the stream change is bit-identical by construction
and the ρ̃ regrouping did not move a visible digit either.  Scaled onto the 402 s free-run ledger:
collocate 87.5 → ~61, integrate 41.1 → ~29, ρ̃ sampling 87.6 → ~37, i.e. ~90 s off a 402 s run.
`ctest -j8` 747 executed + 42 disabled = 789 registered, **0 failed** (the pre-change count, unmoved).

**What was left in the ~150 s, in order** *(item 1 is addressed by Round 4 below)*:
1. **`setup: collocation stream build` (~24–28 s), untouched — the SHELL-BLOCKED KERNEL is its fix**
   (charter item 2, unchanged): `ForPairBox` re-evaluates the contracted radial per CARTESIAN COMPONENT
   PAIR, so a d×d shell pair pays the same two `exp`s 36 times, on a basis whose cost is dominated by Mn
   d shells.  `exp` + `Polarization::operator()` + `uintpow` together are 21.7% of the whole run in the
   round-3 profile.  This is the SAME kernel that owns the over-budget regime (run 64's 4318 s), so it
   pays twice.
2. The scatter's remaining cost is now the value stream plus the scattered destination RMW — i.e. genuine
   work.  Below that lies only fewer POINTS (tighter `GPW_DENSITY_EPS`, or the fold).
3. `setup: becke mesh build` (~12–32 s): `BeckeCutoff` alone was 11.5% of the round-3 profile.

## Round 4 (2026-08-19): the SHELL-BLOCKED box walk — stage A, the stream build

Charter item 2's kernel, staged (user-chosen): **stage A blocks the stream BUILD only**, where the value
screen `epsEff` is uniform across a shell pair and blocking is therefore provably bit-identical; **stage B**
blocks the on-the-fly collocate/integrate branch (the over-budget regime) and is the one that changes
screening semantics.  Stage A landed here.

**THE OBSERVATION, bigger than the charter's.**  The charter names the redundant `exp`s: "`ForPairBox`
re-evaluates the contracted radial per CARTESIAN COMPONENT PAIR, so a d×d shell pair pays the same two
`exp`s 36 times."  True — but reading the kernel, *everything* it does before the value is a **shell**
property, because it reads `radials[i]` only: the pair→level assignment, the screened offset list
(`ForImageOffsets`), the product centre \f$P\f$, `reach`, the fractional half-widths, the ellipsoid
pre-screen, the incremental \f$r\f$ walk, `di`/`dj`/`ri2`/`rj2`, the modulo wrap and the raster index.  A
shell's components differ **only** in `pols[i]` and `ns[i]`.  And `PGData::Init` flattens the basis `Block`
by `Block` (one `GaussianRF` + its polarizations), so a shell is already a contiguous index run — the same
fact the pointwise sweep's shell hoist uses.

**WHAT LANDED.**  `ForShellPairBox(i0,nI,j0,nJ,…)` walks the box once per shell pair and hands each
surviving point the two per-component FACTOR arrays
\f$f^I_a=n_{i_0+a}\,\mathrm{pol}_{i_0+a}(d_I)\,\mathrm{rad}_I\f$ and \f$f^J_b\f$; the caller forms
\f$val=f^I_a f^J_b\f$ for the pairs it wants.  So a d×d shell pair evaluates **2 contracted radials + 12
polynomials per point instead of 72 + 72** — the polynomial hoist is the charter's missing half.  The
product is associated exactly as the unblocked `(ni*pols[i](di)*radI)*(nj*pols[j](dj)*radJ)` was, so values
are bit-identical.  `ForPairBox` is now the \f$1\times1\f$ case of it — one box walk in the codebase, not two.

**AND THE BUILD BECAME ONE PATH.**  The tiering ORDER is part of the result: each pair's fp64/fp32/drop tier
consumes the global budgets in ROW-MAJOR **component**-pair order, so a shell-major *single* pass would tier
differently.  Hence the two-pass shape is now the only shape — (1) count, (2) a serial row-major budget walk
on the counts, (3) build only the tiered pairs — which keeps the budget walk in its original order while
both evaluation passes go shell-major.  That **retired the serial one-pass path outright**, and with it the
transient-bound abort and the streaming fp32 demotion (a two-pass build knows every pair's size before it
materialises anything).  Threading partitions both passes over SHELL pairs, whose component-pair sets are
disjoint, so there is no reduction and serial/threaded stay bit-identical.

**MEASURED** (MnO AFM-II Γ, 3 iterations, 8 threads):

| | before | after | |
|---|---|---|---|
| `setup: collocation stream build` | 28.73 | **17.01** s | **1.69×** |
| `ctest -j8` sweep (the SERIAL build path) | 212 | **203** s | no regression |

Bit-identity verified at the strongest available point — the **stream cache readout is byte-for-byte
unchanged**: MnO `pts64 149999860`, `pts32 496931912`, `offsets 35028`, `runs 63666502`, `meanRun 10.1613`;
the small Si gate `fp64 528 (58084328 pts), offsets 39032, runs 13860311 (mean 4.19069)`.  Same points, same
tiering, same runs ⇒ same streams.  Physics unmoved to every printed digit.  All 24 `GPW.*` kernel gates
pass (replay-vs-on-the-fly adjoint, charge conservation, the budget/self-heal residency gate, the three
stream-fold gates through the new blocking); `ctest -j8` 747/747, 0 failed.

**WHY 1.69× and not 36×, honestly.**  Three ceilings: (a) MnO's shells are s(1)/p(3)/d(6), so the mean shell
pair is ~6–9 component pairs, not 36 — the 36× is the d×d best case, not the average; (b) the counting pass
is a full evaluation pass (the `|val| ≥ eps` screen needs the values), so the build pays two box walks; (c)
what remains per point is now genuinely per-component — the `fI[a]*fJ[b]` multiply, the screen, and the
`Append`.  The residual is the per-component work, which is where it should be.

### Stage B: the on-the-fly arms — and what the profile said about them

`CollocateDensity`'s and `IntegratePotential`'s uncached `else` arms, blocked the same way.  TWO PASSES over
the pair set: a pair with a CACHED stream replays in the original ROW-MAJOR order (so a fully-cached run —
every anchor — accumulates into `rho` exactly as before, bit-identical), and the uncached remainder goes to
a shell-blocked pass.  On the integrate side that pass also picks up the whole STATIC SHARP-FIELD sweep
(`MakeLocalPP`'s κ rule, the explicit-`pairLevels` V_loc ball), which holds no stream cache at all.

**The D-aware box, resolved across a shell pair** (the user's decision): the tolerance stays per component
pair (\f$\varepsilon_{ij}=\varepsilon/|c_{ij}|\f$), but one box must serve them all, so it is sized from the
UNION — the tightest ε present — and each component pair applies its own \f$|val|<\varepsilon_{ij}\f$
screen.  Only ever ADDS sub-eps terms, never drops one, so the no-cut discipline holds.

**★ THE ONE THING THAT ALMOST KILLED IT, AND THE FIX.**  A first cut measured **1.03–1.08×** on the SCF
on-the-fly path against **2.13×** on the local-PP sweep.  The difference is the tolerance: the local-PP
sweep's ε is uniform, the SCF path's is D-aware — and a blocked walk consulting only the *union* tolerance
silently WALKS terms that the unblocked per-component box killed WHOLE via the prefactor early-out
(\f$\mathrm{pf}\ge-\ln\varepsilon_{ij}\f$ ⇒ the entire box is sub-eps).  The kill is recoverable exactly,
because **`pfExp` is a SHELL property** (it depends only on the two \f$\alpha_{\min}\f$ and the centre
separation) while only the tolerance is per-component: `PairPrefactorExp` is now exposed and both shell
drivers pre-filter their component pairs on their OWN tolerance before walking — which also shrinks the
union box to the survivors.

**MEASURED** (over-budget forced with `GPW_STREAM_BUDGET_PTS=GPW_STREAM_BUDGET_PTS_F32=0`):

| | stage A | stage B | + `pfExp` filter | |
|---|---|---|---|---|
| `GPW_SCF.SiliconGammaConverges`, all on the fly | 65.9 | 61.3 | **56.5 s** | 1.17× |
| `GPW_SCF.MnAtomInBoxDChannel` (d shells), all on the fly | 39.4 | 38.1 | **36.8 s** | 1.07× |
| `GPW.LocalPPKappaSelfConverged` (static sharp field, uniform ε) | 15.9 | — | **7.5 s** | **2.13×** |

Energies identical to every printed digit in every arm, cached vs all-on-the-fly included (Si
`Etot=-7.11507`, same 12 iterations either way; Mn `Etot=-14.638`, 9 iterations).  Strictly the union box is
not *provably* bit-identical — a component pair may see a few extra points outside its own box, all of which
fail its magnitude screen — so this is "measured identical", not "identical by construction" as stage A was.
`ctest -j8` 747/747.

**★ AND THE HONEST FINDING: the over-budget regime is NOT geometry-bound, so the charter's kernel is worth
far less there than in the build.**  A `perf` profile of the all-on-the-fly Si run, AFTER blocking:

| | share |
|---|---|
| the `IntegratePotential` / `CollocateDensity` lambdas themselves | 61% |
| `Polarization::operator()` + `uintpow` (the polynomials) | 22% |
| `exp` (the contracted radials) | 10% |

Only the bottom ~32% is what shell blocking removes, and it removes it by roughly the shell width — small
for Si (s/p), still modest for Mn.  The 61% is the irreducible **per-component-pair emit**: the
\f$f^I_a f^J_b\f$ multiply, the magnitude screen, and the scattered `r[idx] += cw*v` / `b += v*V[idx]`.  No
blocking can remove that — it is one operation per (pair, point) by definition.  **So the next lever in the
over-budget regime is fewer POINTS, not a faster kernel:** a looser `GPW_DENSITY_EPS` (the sanctioned
ε-truncation), the T3 stream fold, or a smaller span.  The charter's framing — "make the recompute kernel
CP2K-class and the over-budget cell becomes affordable" — survives for the SETUP sweeps it also covers
(2.13× on the static sharp field, 1.69× on the build) but NOT for the per-iteration scatter/gather, where
the emit was always the floor.  Run 64's 4318 s was never going to fall by a large factor to this kernel;
that regime needs fewer terms, not cheaper ones.

### Step 0 of the real-arithmetic path: exact ±1 Bloch phases at TRIM (2026-08-16)

`BlochPhase(k,n)` now returns the PARITY \f$(-1)^{(2k)\cdot n}\f$ instead of `std::exp` whenever every
component of k is a half-integer — i.e. \f$2k\f$ is a reciprocal lattice vector (Γ and the zone-boundary
points; note a Γ-centred 2×2×2 Monkhorst-Pack mesh is k ∈ {0,½}³, TRIM THROUGHOUT).  One helper, both
phase sites (`BuildImages`'s image list and `CellPhase`'s per-offset closure).

Why it matters beyond tidiness: `std::exp(iπ)` leaves \f$\sin(\pm\pi)\approx1.2\times10^{-16}\f$ in the
imaginary part, so at every zone-boundary k the Bloch S, T, V_loc, KB blocks and the XC-mesh Φ tables
came out *nearly* real and were then symmetrised, having paid complex arithmetic throughout.  Exactly
±1 makes those blocks EXACTLY real, which is what turns "is this block real?" from a tolerance somebody
has to defend into a bitwise fact — the precondition for selecting a real-arithmetic path by TYPE
rather than by threshold.  Gate: `GPW.TRIM_BlochMatricesAreExactlyReal` asserts `MaxImag == 0.0`
(equality, not a bound) for S, T and χ^k(r) at k=(½,0,0) and (½,½,½); `GeneralK_PhaseIsLiveWithImages`
is the negative control at k=(¼,0,0).

The real-arithmetic path ITSELF is not built and its shape is undecided — see the design note below.

## Run 64 (2026-08-17): the OTHER regime, measured — and the conditioning hypothesis REFUTED

The complement of "THE RUNTIME GAP, MEASURED" above, on the same cell at the OVER-BUDGET point:
full-136 spherical span, `MNO_ORTHO_TOL=1e-3`, screens tightened to `GPW_SCREEN_EPS =
GPW_DENSITY_EPS = 1e-12`, 4 iterations, AFM only, 8 threads, under `scripts/memsafe -H 10G`.
Log: `doc/logs/mno_probe_run64_fullspan_tighteps.log`.  76 min total.

**Stream coverage: `dropped 6,948,869,544` points.**  The cache holds 1.0 B (150 M fp64 + 850 M
fp32); ~87% of the demand falls to on-the-fly evaluation EVERY iteration.

| bucket | s | share |
|---|---|---|
| scf: collocate density (pair scatter) | 2952 | 65% |
| scf: integrate-back (pair gather) | 1366 | 30% |
| setup: collocation stream build | 128 | 3% |
| everything else (local-PP, Phi tables, Becke mesh, H_xc, ρ sampling, 1E) | < 100 total | 2% |

**★ THE TWO REGIMES ARE DIFFERENT PROBLEMS, and both are now measured:**

| | streams CACHED (in budget) | streams DROPPED (over budget) |
|---|---|---|
| pair loops | **4%** | **95%** |
| XC mesh (Φ, ρ, H_xc) | dominant | ~0.5% |
| what fixes it | the threading/BLAS work of 2026-08-15/16 (per-iteration 30 s → 7.6 s) | the FAST-RECOMPUTE KERNEL (still unbuilt) |

So the campaign's two items were never competing: item 1's XC-mesh work owns the in-budget regime and
item 2's recompute kernel owns the over-budget one, and neither touches the other's bottleneck.  The
shell-blocked kernel sketched under item 2 aims at exactly the 4318 s above — `ForPairBox` re-evaluates
the contracted radial per CARTESIAN COMPONENT PAIR, so a d×d shell pair pays the same two exps 36
times, on a basis whose cost is dominated by Mn d shells.

**PERSPECTIVE — the recompute path is no longer a wall.**  Run 49 (looser 1e-8 screens, FEWER
functions) could not complete ONE iteration in 80 min.  Run 64 did setup + FOUR at 1e-12 in 76 min
(~18 min/iteration in the pair loops).  Slow, but finite and usable: the kernel is now an
optimisation, not an unblocker, and does not need to pre-empt the real-TRIM work.

**THE PHYSICS: the screen-discipline hypothesis (doc/OpenWork.md item 3) is REFUTED here.**
`iters=4 Efinal=-82.1867 Eamp(last4)=13.66 => UNSETTLED`, m_stag 0.290 → 0.034 → 0.533 → 0.452 (order
survived).  **−82.19 Ha sits 20.7 Ha BELOW CP2K's variational −61.47 on a SUBSET of CP2K's own span**
— unphysical by the same argument that convicted run 58's −67.28.  Tightening eps from 1e-10 to 1e-12
did not prevent the dive, so "CP2K's 1e-14 screens are what let it hold the full span" does not
survive contact.  Caveats kept honest: (a) `MNO_ORTHO_TOL=1e-3` dropped 4 AOs (indices 11, 9, 115, 1),
so this is 132 of 136 — though the near-null modes it was about are still present (λmin 1.96e-05,
cond 3.73e+05); (b) NMAX=4 stopped it mid-trajectory, but a healthy run of this cell is at −61.4 by
iteration 3–4 (run 30: E₃ = −61.394; run 61: −33.0 → −58.4 → −60.9 → −61.35, monotone from ABOVE) and
no variational path descends 20 Ha below the reference and returns.
**Instrumentation gap to fix on any rerun: set `GPW_MNO_VERBOSE=1`** — without it the log carries only
the fingerprint, not the per-iteration table, so the trajectory SHAPE (dive vs stall) is unavailable.
