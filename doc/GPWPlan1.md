# GPWPlan1 — the post-raw-XC execution plan (2026-07-23)

The forward queue, superseding `doc/GPWPlan.md`'s TODO section (that file remains the authoritative
RECORD of the 2026-07 campaign — §0.5(b)/(c)/(f1)/(f2), 0h, the C=2 flip, the raster-policy A/B — and
`doc/GPWHistory.md` holds the deep archive).  Read GPWPlan.md's "Durable pins / invariants" section
before working here: THERE IS NO CUT; no grad-student knobs (policy enums, not numeric dials);
spin-native is the formulation; correct > efficient > end-user > dev > readable.

## Where the code stands (the platform this plan builds on)
- **Raw-XC feed everywhere** (tensors, terms, mixer dynamics): the XC-collapse basin is REMOVED
  (negCharge ≡ 0), ρ_DM ≥ 0 by construction for aufbau densities.
- **Occupation guards** (0h): persistent-hole detection → MOM release/re-capture; NEVER-SILENT
  non-aufbau warning at run end.  Detection is complete; the REMEDY for genuinely near-degenerate
  frontiers is item 4b below.
- **C=2 default** (density resolved at its own product exponent); ladder top-rung gate decoupled
  (`kRungGateC=8`).
- **RasterPolicy A/B measured + DEFAULT FLIPPED to BallOnly on the GPW surface** (user 2026-07-23:
  ≈1 mHa at/above the C=2 floor, CP2K's bet confirmed; PW path keeps AliasFree pending its own A/B;
  three exact-quadrature kernel gates pinned {.raster=AliasFree}; NaF production anchor −24.4304 =
  0.8 mHa from CP2K).
- **Runtime vs CP2K (single-thread us, threaded them)**: Si 45 s AliasFree / 7.2 s BallOnly vs 3.5 s;
  NaF 63 s at auto vs 5.8 s.  The structural gap is closed; BallOnly + OMP ≈ parity.
- **Known pathology, one family left**: near-gapless occupation flapping at particular
  discretizations (NaF Ecut=160, three sightings; the guard WARNS but cannot fix).  Item 4b cures it.

## The sequence (agreed 2026-07-23; dependencies annotated)

### 1. PARAM-STRUCT GRADUATION (unblocks the user's experiments + the BallOnly default decision)
Move settled knobs out of env into typed structs; grade **Standard** vs **Advanced** (Advanced =
sensible defaults, touched only for pathological cases — and the long-term goal is that the guards
make pathological cases self-correcting: the software behaves like an expert system).
- Structs over JSON (user decision): compiler-checked, IDE-discoverable, documented at the
  declaration.  If config FILES are ever wanted, serialize the structs — the dependency runs one way.
- Shape sketch: `SCFParams` keeps the Standard surface; nested per-seam structs for Advanced
  (`MixerParams` α/G0/Pulay, `OccupationParams` policy enum {Aufbau, MOM(start), Fermi(kT)} + guard
  tunings (hole-persistence K, release cap), `GridParams` densityEcut/cutoffFactor/RasterPolicy).
- The `NAF_*`/`GC_*` env overrides DELETE (recipes hard-coded per test); verification instruments
  (`GPW_MGRID_ECUTS`, `GPW_RELCUTOFF`, `GPW_RASTER_POLICY` until promoted, `GPW_ILLCOND_*`) stay env
  but are documented as instruments; ops valves (`GPW_STREAM_BUDGET_*`, `GPW_OMP_THREADS`) stay env.
- Present everything clean and human-readable in the Si and NaF integration tests (the tests ARE the
  user documentation).
- **Fold in: the BallOnly default decision** (RasterPolicy lands on the struct; flipping the default =
  one field + the anchor re-pin wave, calibration already recorded in GPWPlan.md 0.5(a)).

### 2. PER-SYSTEM ITERATION OUTPUT (virtual DisplayColumns/DisplayColumnHeaders)
Rides the same presentation cleanup as item 1 — do while the SCF surface is open.
- Atoms/Molecules: {E, [F,D], Δ[F,D], Δρ, virial, ρ_mix, accel, config}.
  ρ_mix fixed # of chars to identify itself (Lin,Ker,Pul) and a number
  accel fixed # of chars to identify itself (Null,DIIS,GMD) and a number
  config 1 char,  put * if the config changed (full configs are way too long)
  open to suggestions on any of this
- Solids/PP: {E, [F,D]?, ΔE, Δρ, ρ_lost, ρ_mix, accel, config, gap, ?? 
  DROP the virial (2+V/K assumes Coulombic homogeneity; GTH local + KB projectors break
  it — the idealVirial fudge retires) and show the health instruments that already exist as ad-hoc
  lines: gap (hole-flagged, from the fixed HomoLumo)
- Absorbs the `ReportBandGap()`/`ReportGridCharge()` process-wide flags into the per-system display.

### 3. CACHE2/3 BYTE-BUDGET LRU (GPWPlan §5 design, user-approved)
Independent; land BEFORE item 4a because the diffuse-basis campaigns are exactly the workload that
bloats the geometry caches on image clones (the reason `ClearGeometryCaches()` + `Overlap3CStream`
exist).  Byte budget + LRU eviction, policy selected per scope BY THE ALGORITHM (lattice = scoped
size-1 budget preserving the const& contract; molecular = generous default); per-cache RAM in the
end-of-run report; the clear-based band-aids then RETIRE.

### 4. ROBUST PHYSICS DEFAULTS — the two remaining robustness fronts
**4a. Diffuse-basis robustness (GPWPlan §1, "a grad student can add diffuse functions at will"):**
rectangular V through the periodic stack (truncated ortho → per-k orbital dims through Crystal_EC /
cDM_CD / collocation), LASolver gap-detection auto-tol, NEVER SILENT drops, promote the working
recipe into facade defaults.  Gates: full-basis NaF == SR/SR2 ± the dropped-mode mHa; Si anchors
untouched.
**4b. FERMI SMEARING (the occupation seam's missing policy — added 2026-07-23, user-approved).**
NOT major surgery, by construction:
- The seam exists: `TOrbitals::TakeElectrons` / `tIrrepWF::FillOrbitals`.  New policy: occupations
  f_i = 1/(1+exp((ε_i−μ)/kT)), μ by bisection on Σf = nₑ.  Fractional f already flows through the
  density build (degenerate levels use it today).
- The entropy NEVER touches H: at fixed T all operators are unchanged (f enters only D).  The Mermin
  term is a SCALAR from the occupations: S = −k Σ w[f ln f + (1−f)ln(1−f)]; EnergyBreakdown reports
  E, −TS, A = E−TS (+ optionally the ½(E+A) T→0 extrapolation); the iterator's E-flat gate reads A.
- **Design (user 2026-07-23): −TS is just ANOTHER TERM** — a Hamiltonian-object term contributing to
  `EnergyBreakdown` — so `SCFIterator::Iterate` keeps gating on "Etotal" (which quietly IS the free
  energy A when smearing is on), and the virtual DisplayColumns/DisplayColumnHeaders (item 2) label
  it honestly per system.  No iterator surgery at all.
- Increment 1: per-block μ (covers Γ-only — the three-sighting NaF Ecut=160 repro).  Increment 2:
  GLOBAL μ across k-blocks (structural: today each block fills to a fixed per-block nₑ) — timed with
  the IBZ track, like item 5.
- Gates: (i) gapped regression — Si at kT ≪ gap reproduces −7.11501 tightly (smearing inert where it
  should be); (ii) the cure — NaF Ecut=160 converges to −24.431-class instead of −24.078 flapping;
  (iii) A decreases monotonically where E need not.
- With 4b landed, the 0h warning's advice ("check the occupation recipe") gains an actionable
  default: the expert-system loop closes for this pathology family.

### 4c. LIBCINT LATTICE ENGINE (added 2026-07-23, user): `BasisSet::Molecule::Engine::LibCint` is
much faster than the MnD kernels — teach `PG_LibCint` to realise `Molecule::LatticeSum1E` so the
periodic/GPW path can select it (GPW itself is unchanged either way — the engine is a molecular-side
switch; the GPW_UT header already anticipates exactly this).  Slot after 4a/4b or opportunistically —
it is orthogonal to both; the collocation STREAMS stay MnD (they are ours), so the win lands on the
1E/analytic-KB/3C build side.

### 5. B_ij(R) k-INDEPENDENT 1E MEMO (GPWPlan 0.5(d))
Deliberately LAST: its payoff only materializes on multi-k runs — time it with the IBZ/space-group
track's arrival (the qchem7 `lattice-3d-spacegroup` work), so it never carries untested-in-anger.
Design pin (user): cache B(R), never M(k) — "keep k out of the key".

## Parked/background (unchanged from GPWPlan.md)
0.5(e)'s runtime aspect folded into item 3.  0i analytic V_local LONG (robustness; fold the core
charge into PW_Hartree's G-space solve) — NOTE (user question 2026-07-23): the ladder's L>0 levels +
top rung currently serve TWO clients, the collocation completion (REL-rule) and the V_long κ-ruled
grid sweep; after 0i the sweep client disappears and whether the completion rung still pays at
BallOnly+C=2 becomes a measurable rung-gate question.  §2 low-q multi-species Si/NaF/CsI cross-validation.
§3 CP2K reference library growth.  §5 remaining cleanups (Vxc-fit ISP ctor lands with GGA; cMesh
unification; DRY PP adapters; periodic external PP).
