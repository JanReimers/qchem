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

---

# Forward roadmap (agreed 2026-07-26 — toward an honest metal + Fermi-smearing test)

**0. REPORTING FEATURE** — build `doc/RunReportPlan.md`.  `qchem.Reporting`: the report IS json, a generic
console renderer (layout INFERRED from json structure — table vs tree; NO per-section renderers), a global
sink (`GlobalReport` keyed + key-free `CurrentRunReport`), incremental section-by-section rendering, detail
level console-only.  Consolidates the scattered `[GPW grid]`/`[overlap S]`/cache-RAM/SCF-settings prints.

**1. VALENCE-BASIS-GEN CLI + an Al basis** — extract `IntegrationTests/ValenceBasisGen_UT.C` into a standalone
`CLIapps/valgen.C` (mirror `CLIapps/scfrun.C`; thin arg-parser over the existing `qchem.ValenceBasisGen`
library — `GenerateValenceBasis`/`GenerateSeedDensity`/`AssembleBasisFile`).  Then generate a valence basis for
a metal (Al).  **The CLI uses the Reporting framework (step 0) for its console output — the first dogfood.**

**2. DEGENERATE-SHELL METAL TEST (works with Increment 1, per-block μ)** — FCC Al @ Γ (3s²3p¹: the degenerate
3p triplet is partially filled, aufbau can't pick a p → won't converge without smearing), same physics as the
passing Si-3p test.  Add **annealing** (a kT schedule, ramp T→0 in steps, re-seed each stage).  NOTE: Γ-only, a
degenerate open shell — NOT yet a Fermi-surface metal (that's step 4).

**3. INCREMENT 2 — GLOBAL μ ACROSS k-BLOCKS** (structural: today each Bloch block fills to a FIXED per-block
nₑ; a metal needs ONE μ across the BZ with charge sloshing between k-points).  The enabler for a true metal.
Time with the IBZ/space-group track (qchem7 `lattice-3d-spacegroup`), like item 5.

**4. MULTI-K Na — the honest metal test** — FCC and/or BCC Na (1 valence e, half-filled conduction band): a
real Fermi surface, one μ across the BZ, smearing + annealing.  Gate (iv), done properly.  (May want a more
diffuse metallic Na basis from step 1 than the ionic-NaF-tuned `valence_lowq` Na.)

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
  for any hard/near-gapless landscape, worth exposing as such.

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
