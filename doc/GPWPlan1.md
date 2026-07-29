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

---

# Forward roadmap (agreed 2026-07-26 — toward an honest metal + Fermi-smearing test)


**1. VALENCE-BASIS-GEN CLI + an Al basis** — extract `IntegrationTests/ValenceBasisGen_UT.C` into a standalone
`CLIapps/valgen.C` (mirror `CLIapps/scfrun.C`; thin arg-parser over the existing `qchem.ValenceBasisGen`
library — `GenerateValenceBasis`/`GenerateSeedDensity`/`AssembleBasisFile`).  Then generate a valence basis for
a metal (Al).  **The CLI uses the Reporting framework (step 0) for its console output — the first dogfood.**
DONE (2026-07-28): valgen shipped (`--q`/`--semicore`/`--electrons`/`--shell`/`--iterations`, valence-default PP
selection, exponent+population basis-usage heat map, virial gate OFF for PPs).  Al q3 basis in `valence_lowq.bsd`.
  - **F⁻ still does NOT converge** — and it is NOT the (now-disabled) virial gate.  The neutral-F/F⁻ valence SCF
    limit-cycles: `[F,D]` oscillates 0.2↔4.5, energy swings −18.9↔−21.0, hits `NMaxIter=20`.  A genuine DIIS
    instability on the diffuse anion (energy anchor ~−21 is still fine, so the emitted basis is OK).  Candidates
    to try: more iters, MOM, Kerker/relax tweak, or Fermi smearing — revisit with step 2's smearing/annealing.
  - **Open-shell atoms are ml-SPLIT (Hund), but atomic op(r) is RADIAL-ONLY** — Al 3p¹ reports irrep `p {-1}`
    (pol U): `Atom_EC` puts the single unpaired p e⁻ in a definite mₗ=−1 (Hund, lowest mₗ first).  BUT the atomic
    evaluator `Radial::operator()(rvec3_t)` = `gaussian(norm(r),l,es,ns)` uses ONLY |r| — no Ylm, direction
    discarded — so op(r)/ρ(r) is the RADIAL profile regardless of mₗ.  ⇒ the emitted seed ρ(r) is SPHERICAL (not
    direction-biased); the mₗ label only enters the angular ERI couplings (DirectAk/ExchangeAk) of the SCF energy.
    So the seed is the spherical radial density of the mₗ=−1 SYMMETRY-BROKEN solution, not of the true fractional
    (⅓-per-mₗ) ground state — usually fine for a seed, but the fractional/spherical version is what step 2's
    smearing gives.  Deeper caveat (known weakness): the atomic op(r) is "fake" (no angular dependence), so it
    cannot represent real-space angular structure — matters when seeding mixed-valence ions (e.g. LiMn₂O₄ Mn³·⁵⁺).

**2. DEGENERATE-SHELL METAL TEST (works with Increment 1, per-block μ)** — FCC Al @ Γ (3s²3p¹: the degenerate
3p triplet is partially filled, aufbau can't pick a p → won't converge without smearing), same physics as the
passing Si-3p test.  Add **annealing** (a kT schedule, ramp T→0 in steps, re-seed each stage).  NOTE: Γ-only, a
degenerate open shell — NOT yet a Fermi-surface metal (that's step 4).
- FYI Al basis set in BasisSetData/valence_lowq.bsd ... valgen CLI app all ready if you need to regenerate (problems with diffuse functions is a common scenario).
- **DONE (2026-07-28).**  Two IntegrationTests/GPW_SCF_UT.C tests (physical a=4.05 Å=7.653 au, 1-atom FCC
  primitive cell, Zion=3):
  - `AlFCCDegenerateShellAufbauStalls` — the MOTIVATION: aufbau settles the energy (|ΔE/E|~1e-13, gap
    column ~0 = metallic) but |Δρ| floors (~5e-4) and never converges (the degenerate-3p rotation).
  - `AlFCCAnnealedMetal` — the CURE: a descending kT schedule {0.02, 0.01, 0.005} Ha via a new
    `RunGpwAnnealed` driver (re-seeds each stage from the prior converged density; same basis/grids/Ham, only
    kT changes → direct DM transfer, bs outlives the stages).  Every stage converges; the INTERNAL energy
    E=A−(−TS) is kT-independent to ~1e-7 across all three (−1.92115) — the physical T→0 answer; free energy
    A=E−TS at kT=0.005 is −1.93466, −TS<0 (gate iii).  Annealing cheapens the cold end: kT=0.005 cold-start
    ~44 iters → re-seeded ~17.
  - **Basis: Al block added to `valence_lowq_sr.bsd`** (dropped the diffuse s 0.04 + p 0.05).  The full
    `valence_lowq` Al p 0.05 makes the primitive-cell Bloch overlap rank-deficient (VetBasis drops the p
    triplet, ABORTs) AND explodes the collocation grid (reach ~21 au) — the SIPP→SIPP_SR / NaF-SR lesson.
    SR Al @ a=7.653: λ_min=3e-3, cond=2e3, 17 functions.
  - **Annealing is now a driver (`RunGpwAnnealed`), not yet a first-class library capability** — the
    "Annealing as a general capability" future-consideration below still stands.
  - **FINDING — GDM tail rung is INCOMPATIBLE with Fermi smearing (measured on this test, 2026-07-28).**
    Ran the annealed Al on the DIIS→GDM `Ladder` (`ScheduleSignal::EnergyChange`).  DIIS descended cleanly to
    the physical −1.97521 ([F,D]~4e-3), handed to GDM, and GDM's FIRST geodesic step was a persistent FALLBACK
    (`GPW_GDMTRACE`: best(Et−Ecur)=+0.42, all 12 backtracks uphill) → occupations flipped (cfg `*`) and the run
    never recovered.  ROOT CAUSE (not grids / not non-variationality — DIIS converges the same E[ρ]): GDM builds
    its geodesic DIRECTION from the fixed-occupation electronic gradient `[F,D]`, but line-searches the FREE
    energy A=E−TS with occupations Fermi-refilled per trial (`DirectMinStep`).  Under fractional occupation A is
    stationary where the SMEARED gradient (not `[F,D]`) is zero, so `[F,D]≠0` at the DIIS fixed point and GDM's
    direction is not downhill for A.  FIX (deferred, a real change): give GDM the occupation-response term so its
    search direction is the free-energy gradient dA/d(rotation), not `[F,D]` — then GDM could tail-polish a
    smeared solid.  Until then: DIIS (or Kerker/Pulay) for smeared runs; GDM only at fixed integer occupation.

**3. INCREMENT 2 — GLOBAL μ ACROSS k-BLOCKS** (structural: today each Bloch block fills to a FIXED per-block
nₑ; a metal needs ONE μ across the BZ with charge sloshing between k-points).  The enabler for a true metal.

**Approach agreed 2026-07-28 (two scoping decisions):**
- **DECOUPLE from IBZ — full unreduced Γ-centred mesh FIRST** (weights = 1/N_k).  Global-μ physics needs no
  symmetry reduction; building it on a full mesh now buys multi-k experience + a working metal fill sooner, and
  IBZ rides in later (only the weights change).  (Was "time with the IBZ track"; deliberately un-timed.)
- **EC surface = a MODE on `Crystal_EC`**, paralleling the existing `UsesAufbau()` branch — a "global total-N"
  ctor + a virtual (e.g. `UsesGlobalFermi()`); NOT a new `Metal_EC` type, NOT (yet) an `OccupationPolicy` enum
  (that's the item-1 graduation).

**Design (where it lives — traced 2026-07-28):**  Today the composite WF's `FillOrbitals` loops k-blocks
independently (`w->FillOrbitals(ec->GetN(irrep))`); under smearing each block runs its OWN μ-bisection
(`TOrbitals::TakeElectronsFermi(ne_k,kT)`) to hold exactly `ne_k`.  Item 3 lifts that bisection from per-block
to ONE μ across the mesh.  The coupling is SMALL: Fock build / diagonalization / DIIS stay per-block; only the
OCCUPATION step goes global.  It reduces EXACTLY to today at Γ (one block, weight 1 → global μ ≡ per-block μ),
so the Al-Γ tests cannot regress.  Aligns with `doc/SCFStrategyPlan.md` §5 (occupation = first-class seam; the
μ-solver "per k-block → global" is the named upgrade; §5 already pre-flags the GDM×smearing pitfall item 2 hit).

**Increments — ALL DONE 2026-07-28 (full-mesh scope; commits d9a416e0, e474cf67, + inc 3):**
1. ✅ **Split the Fermi primitive.**  `TakeElectronsFermi(ne,kT)` = free `FermiLevel(e,g,target,kT)` μ-solver +
   virtual `SetFermiOccupationsAtMu(μ,kT,eShift)`.  The solver is weight-agnostic (g = degeneracy OR
   w_k·degeneracy), so ONE bisector serves per-block and cross-k.  Γ bit-identical (commit d9a416e0).
2. ✅ **`Crystal_EC` global mode.**  `globalFermi` flag (defaulted off) + `UsesGlobalFermi()`; Nval is the
   whole-mesh total.  Purely additive (commit e474cf67).
3. ✅ **Composite global fill** `FillOrbitalsGlobalFermi`: gathers `(ε_i,k, w_k·g_i)` across the channel via the
   SAME `GetQNs().sym->GetWeight()` the density uses, bisects one μ on `Σ_k w_k Σ_i g_i f = N_total`, then
   `tIrrepWF::FillOrbitalsAtMu(μ)` on each block.  `−TS` BZ-weighted in `tIrrepWF::GetEntropyTerm` (w_k·(−TS_k)),
   consistent with the BZ-weighted E.  Branch in `FillOrbitals` next to `UsesAufbau()`.
4. ✅ **Validated on Al 2×2×2** — committed test `GPW_SCF.AlFCCMetalGlobalMu` (8-point Γ-centred mesh):
   - **Γ invariant**: global μ ≡ per-block Fermi to 1e-12 at a single k (`DISABLED_AlGlobalMuExperiment`).
   - **Charge conserved**: `Σ_k w_k n_k = 3.00000000` exactly (the weight-consistency guard holds).
   - **Metal signature**: global μ CONVERGES to A=-2.1168 where per-block filling (AL_GLOBAL=0) forces 3 e⁻ at
     every k and lands non-converged garbage A≈-0.46 → charge MUST redistribute between k-points.
   - **Dispersion**: k-sampling binds (-1.92 Γ-only → -2.12 at 2×2×2).  Full ctest 611/611.

**Multi-k risks — resolved:** (a) weight consistency — the μ constraint reads `GetQNs().sym->GetWeight()`, the
SAME accessor `TOrbitalsImp::GetChargeDensity` scales D by; charge conserves to 1e-8, confirming it.  (b)
complex-k — occupations/eigenvalues real, so the global bisection is real arithmetic (2×2×2 has complex blocks
and works).  (c) spin-polarized metals share ONE μ across BOTH spin channels — STILL DEFERRED (this is
spin=None, one channel; a magnetic metal needs the μ-solve to span both `itsSpinWFs` channels).  (d) kT must
exceed the inter-k level spacing near E_F (the item-2 "smear wider than the splitting" finding).

**Two follow-up findings (2026-07-28, from inspecting the AlFCCMetalGlobalMu run):**
- **Eigen-table display was never tested on a metal — FIXED.**  `tUnPolarizedWF::DisplayEigen` broke the level
  list at `e>0.0` (a MOLECULAR idiom: bound states sit below the vacuum level at 0).  In a solid the energy
  zero is arbitrary (PP G=0/alignment) so the Fermi level can be POSITIVE — Al's μ=+0.28 — and the break hid
  every occupied level above 0, showing only the one negative-energy Γ level.  Fix: break on OCCUPATION
  (`occ<1e-6`), and show FRACTIONAL occupations (the metal now honestly displays 1.90/2, 0.10/2 at the Fermi
  surface).  Atom/molecule output byte-identical (integer occ, bound levels).  Instrument `GPW_METALTRACE`
  dumps μ + per-k n_k (the charge sloshing).
- **Small `Eee` on the metal is PHYSICAL, not a bug.**  GPW's Hartree is the G≠0 Poisson solve `4πρ̃/G²`; the
  G=0 self-energy is dropped into `Ealign` (the neutralising-background alignment).  So `Eee` measures only the
  density's NON-uniformity — and a metallic (near-uniform) density has small `Eee` (Γ-only 0.020 → 2×2×2 0.006
  as the density metallizes; a uniform gas → 0).

**Remaining for item 3:** IBZ integration (weights from the symmetry-reduced mesh — the qchem7 track; the fill
is weight-agnostic so only the weights change) and the spin-polarized global μ (both spin channels, one μ).

**4. MULTI-K Na — the honest metal test** — FCC and/or BCC Na (1 valence e, half-filled conduction band): a
real Fermi surface, one μ across the BZ, smearing + annealing.  Gate (iv), done properly.  (May want a more
diffuse metallic Na basis from step 1 than the ionic-NaF-tuned `valence_lowq` Na.)

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
- KEEP the uniform FFT grid for **Hartree** (the G-space Poisson is the reason GPW is fast; a Becke grid can't
  do it).  So this ADDS a second grid for the XC only; it does not replace anything.
- ADD an atom-centered **Becke/Voronoi** grid for the XC quadrature: collocate ρ on the Becke points, evaluate
  ε_xc pointwise, `⟨i|V_xc|j⟩ = Σ_g w_g χ_i χ_j V_xc(g)` — the STANDARD molecular-DFT quadrature.  Dense radial
  near each nucleus (resolves the sharp V_xc), sparse far out (a diffuse tail costs almost nothing) → the
  diffuse explosion simply does not happen (point count ~ atoms × radial × angular, NOT the diffuse reach).
- This is the scale-decomposition (fine near atoms) pictured cleanly.  It is **lighter than GAPW** (which
  augments BOTH Hartree and XC near cores + a compensation charge — deferred, out of first-pass scope): Becke
  is XC-quadrature ONLY, Hartree stays G-space.

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
