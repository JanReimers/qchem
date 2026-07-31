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
  pair-loop thread-safety bug (RESOLVED same day: Blaze SMP's global guard flag, not our code — below).

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

# Becke XC grid — remaining increments (the build/gate/perf record is in DONE above)

*(Original motivation, design, and gate spec served and condensed into DONE; full narrative in git history.
The framing that stands: the Becke grid gives the XC — the one pointwise-nonlinear, sharp-at-core term —
its near-ideal grid across every basis set, diffuse included, while Hartree stays on the G-space Poisson.)*

- **GPW pair-loop THREAD SAFETY — ✅ RESOLVED (2026-07-31).**  The Cache2/3 suspicion was WRONG (the
  threaded pair path never touches them — audited).  gdb caught the real thrower: **Blaze** —
  `std::runtime_error("Nested parallel sections detected")`.  Mechanism: the global `-fopenmp` (2fd152c8)
  defines `_OPENMP` in every TU, which silently flipped ON Blaze's shared-memory parallelization
  (`BLAZE_USE_SHARED_MEMORY_PARALLELIZATION` defaults to 1); Blaze then wraps EVERY `smpAssign` — even
  tiny below-threshold vector copies like ForPairBox's `rvec_t` exponent/coeff copies — in
  `BLAZE_PARALLEL_SECTION`, whose guard flag is a **global non-atomic bool**; two pair workers race it,
  the loser throws, and an exception escaping an OpenMP region is `std::terminate` ("terminate called
  recursively").  FIX: (a) `BLAZE_USE_SHARED_MEMORY_PARALLELIZATION=0` global define (CMakeLists, beside
  the OpenMP block) — restores the serial-Blaze semantics every committed anchor was calibrated on
  (never-nested-OpenMP strategy); (b) throw CONTAINMENT in both pair loops (first worker exception
  captured, rethrown serially after the region — any future worker throw surfaces as a normal exception,
  never terminate).  VERIFIED: TSan+Archer (`build/TSan`, 4-thread NaF, setup+3 iters, exercising stream
  replay + on-the-fly fallback + memo recording) = ZERO project races (only libomp-internal false
  positives); a full threaded NaF run reproduces the serial trajectory to ~1e-10 through 45 iterations
  (the accepted cross-pair ULP reduction drift only); no aborts.
  **SURFACED while verifying (pre-existing, serial == threaded, NOT a threading artifact):**
  `DISABLED_NaFRocksaltGamma` as committed (Becke XC hardcoded ON, line ~1041) fails the −24.4304 anchor
  IDENTICALLY serial and threaded — an unoccupied ε dives to −50 Ha at iters 10–12 (guard releases, it
  re-dives at ~45, run ends non-aufbau at −19.594).  That is a Becke-on-this-anchor config question for
  the "Becke as default" decision below, not thread safety.
- **0i analytic V_loc-long** (the other measured diffuse-setup lever): fold the smooth Gaussian core
  charge into PW_Hartree's G-space Poisson solve — NO real-space sum at all; deletes the 180 s dominant
  ledger bucket.  (The measured short-part analytic-vs-grid ratio, 0.37 s vs 180 s over the same 791-cell
  enumeration, is the empirical case.  CP2K's numeric-long is free for THEM only because it rides their
  per-iteration total-KS-potential integrate-back, which our G-space assembly doesn't have.)
- **Coarse-end routing calibration** (`kMinLevelN=3` "very coarse for a diffuse lobe" — the standing
  suspect after the HartreeOnly falsification; valves `GPW_ROUTING`/`GPW_RELFIELDSHARP` are in place).
- **Becke as the DEFAULT for (diffuse) bases** — runtime-viable since the engine; decide after the
  thread-safety + 0i items settle the setup story.
- **Spin-native** (`PW_XC`/`PW_XC_Becke` are unpol today; spin-native is the formulation pin).
- **Site-group-adapted angular grid** — with the measured design pin: orbits must AVOID special (bond)
  directions (the Lebedev cube-corner lesson); the exact required degree is deducible from the site group.
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
