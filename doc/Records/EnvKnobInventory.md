# Environment-knob inventory (D-ENV step 1) — 2026-10-03

**RECORD, not a queue.**  The plan this feeds is `doc/CleanCode.md` **D-ENV** (design agreed 2026-10-03: four tiers,
facade-constructor slices, an option registry, a reproducible input deck).  This file is the inventory that plan's
Step 1 asked for: every environment variable the code reads, what it does, which tier it belongs to, whether it can
move a number, and a PROPOSED disposition.  **The user reviews the dispositions** (§5) before any code moves.

Method: `grep 'getenv("' src/` (52 call sites, 41 distinct names after de-duplication, excluding the five
thread/OpenMP names already settled by D-THREADS), plus the 9 names `RunPolicy` resolves, plus the two places
outside `src/` that read the environment as a poor-man's input deck (§4).  Names mentioned only in comments are not
counted (565 distinct `GPW_*`/`QCHEM_*` strings appear in `src/`, nearly all in comments and docs).

Tiers (D-ENV): **T1** calculation choice · **T2** numerical-method policy · **T3** resource (never changes a number) ·
**T4** diagnostic / A-B hatch.  "Moves numbers?": **no** = output/timing only; **~ulp** = bit-level or
sub-1e-12; **yes** = changes a converged value or the trajectory.

## 1. Tier 4 — diagnostics (print or record something; the answer is identical) — 22 names

| knob | read at | what it does | proposed disposition |
|---|---|---|---|
| `GPW_BECKE_ATOMS` | `Structure/Imp/UnitCell.C:300` | per-atom Becke share/kept/dropped breakdown | keep as diagnostic |
| `GPW_BECKE_COUNT` | `UnitCell.C:145` | census of Becke point costs | keep as diagnostic |
| `GPW_DM_RANK` | `ChargeDensity/Internal/Imp/IrrepCD.C:212,385` | rank/PSD census of D, IPR of the Cholesky orbitals | keep (feeds the pin-21 canary); candidate for `gpwprobe` |
| `GPW_EH_TRACE` | `Hamiltonian/Internal/Imp/PWTerms_Hartree.C:145` | E_H G-space pairing trace | keep as diagnostic |
| `GPW_FIELD_SPECTRUM` | `BasisSet/PlaneWave/Evaluators/Imp/Evaluator.C:94` | cumulative |c|² spectrum of a fitted field | keep as diagnostic |
| `GPW_GDMTRACE` | `SCFIterator/Imp/SCFIterator.C:544` | line-search energy trace | keep; should become `SCFParams::Verbose` level |
| `GPW_INTEGRATE_CENSUS` | `PG_Cart_MnD/Evaluator.C:2033` | NEW/REPEAT census of integrate-back calls | keep as diagnostic |
| `GPW_KERKER_SPECTRUM` | `ChargeDensity/Internal/Imp/FieldMixer.C:51` | residual power binned by G/G0 | keep as diagnostic |
| `GPW_LOCALPP_RELCUTOFF` | `GPW_Evaluator.C:90,1238` | times the local-PP build and prints | keep as diagnostic |
| `GPW_MESH_ORTHO` | `UnitCell.C:396` | plane-wave orthogonality error of the XC mesh, binned in \|ΔG\| | keep as diagnostic |
| `GPW_METALTRACE` | `WaveFunction/Internal/Imp/CompositeWF.C:565`, `IrrepWF.C:193` | shared-μ fill per-block trace | keep; `SCFParams::Verbose` level |
| `GPW_NL_PER_L` | `Hamiltonian/Internal/Imp/PWTerms_PP.C:121` | banks per-l projector blocks (I0 diagnostic; spherical-lattice plan) | keep as diagnostic |
| `GPW_PHI_SPARSITY` | `BasisSet/Imp/DeltaFit_IBS.C:52` | Φ-table block-sparsity ceiling | keep as diagnostic |
| `GPW_RHO_NEGATIVE` | `DensitySampler_Singles.C:195` | negative-ρ census per route | keep as diagnostic |
| `GPW_RSS_TRACE` | `PG_Cart_MnD/Evaluator.C:163` | resident-memory breadcrumbs in the lattice sum | keep as diagnostic |
| `GPW_XC_ALPHA` | `DensitySampler_Singles.C:161` | prints the XC mix's effective α each step | keep as diagnostic |
| `GPW_XCROUTE` | `DensitySampler_Pair.C:220` | prints which V_xc route fires | keep as diagnostic (R2.21 wants a FORCE, not a print) |
| `QCHEM_ANGMESH_DEBUG` | `Structure/Lattice_3D/SymmetrizeMesh.C:522` | NNLS site-adapted angular-mesh debug print | keep as diagnostic |
| `QCHEM_DUMP_H` | `WaveFunction/Internal/Imp/IrrepWF.C:52` | ‖F‖, trace, max imaginary part of each H | keep as diagnostic |
| `QCHEM_MOM_SCORES` | `IrrepWF.C:165` | sorted head of the MOM scores at each fill | keep as diagnostic |
| `QCHEM_SITE_MOMENTS` | `DensitySampler_Singles.C:440`, `SCFIterator.C:269` | integrated site moments each iteration | **candidate for promotion**: a measured *integrated observable* users want; belongs in the report, not an env flag |
| `QCHEM_U_TRACE` | `Hamiltonian/Internal/Imp/Hubbard.C:760` | per-refresh +U occupation line to stdout | keep; candidate for `SCFParams::Verbose` |

None of these changes a number.  Proposed: move them behind **one** documented diagnostics reader
(`qchem.Diagnostics`: `Enabled("dm_rank")`), registered with a one-line description, so a typo'd name is caught and the
list is greppable in one place; they stay OUT of the deck.  The `*TRACE`/verbose ones are really verbosity levels and
should fold into `SCFParams::Verbose`.

## 2. Tier 4 — A/B hatches and verification instruments (a different ROUTE or tolerance; can move numbers) — 12 names

| knob | read at | default | what it does | moves numbers? | proposed disposition |
|---|---|---|---|---|---|
| `GPW_CONTRACT_CUBE` | `Evaluator.C:1171` | on | `=0` selects the reference box walk | ~ulp (agree 1e-14, unit-tested) | keep as hatch; **now has the `ContractCubeOverride()` test hook** — the env read can go to the registry |
| `GPW_EXP_RECURRENCE` | `Evaluator.C:838` | on | exp recurrence vs direct | ~ulp (10 s.f. unchanged) | keep as hatch |
| `GPW_SPHERE_SCREEN` | `Evaluator.C:618` | on | spherical vs box log-cut in the collocation | yes (screen) | keep as hatch |
| `GPW_LONG_SWEEP` | `GPW_Evaluator.C:1282` | off | old long-range local-PP sweep (A/B for the custom G-ball) | yes (5.5 µHa) | keep as hatch; candidate to DELETE once the A/B is retired |
| `GPW_XCGRID_NOSELECT` | `Mesh/XCPolicy.C:379` | off | `=1` disarms the Becke/uniform selector (always Becke) | yes (mesh) | **duplicate of a typed option** (`xcMesh.cellKind` Becke/Uniform) → delete the env |
| `GPW_RASTER_POLICY` | `GPW_Evaluator.C:266` | AliasFree | `ball` forces BallOnly on every grid | yes | **duplicate of typed `SolidCalcOptions::raster`** → delete the env |
| `GPW_BECKE_ANG` | `Mesh/XCPolicy.C:271` | Lebedev ≥ degree 29 | forces the angular scheme | yes | **duplicate of typed `MeshParams::angular`** (D5: applies only to args passed `<0`) → fold into D5 |
| `GPW_MGRID_ECUTS` | `GPW_Evaluator.C:855` | — | explicit multigrid cutoff list (CP2K-matching: CUTOFF/3^i) | yes | instrument for CP2K-parity runs; promote to a typed `Grids::ladder` ONLY if parity runs stay routine |
| `GPW_RELCUTOFF` | `Evaluator.C:455` | 0 | forces the absolute pair→level rule at a given κ | yes | instrument (CP2K-matching); same disposition as above |
| `GPW_XC_DM_MIX` | `DensitySampler_Singles.C:146` | — | overrides the XC-from-D mix α | yes | instrument tied to `GPW_XC_DM_SOURCE`; candidate for deletion with N4 |
| `GPW_XC_DM_BOOST` | `DensitySampler_Singles.C:151` | 1 | scales the XC-from-D mix α | yes | same |
| `QCHEM_SPINBLIND_KERKER` | `ChargeDensity/Imp/DensityMixer.C:67` | off | single-map mixer on a polarized density (re-creates the MnO collapse; pinned by a test) | yes | **test-only valve**; delete the env, keep a test-only constructor argument |

Three of these shadow a typed option (`xcMesh.cellKind`, `raster`, `MeshParams::angular`): the typed field and the env var
can disagree, and the banner shows only one.  Deleting the env twin is the cheapest win in the whole inventory.

## 3. Tier 2 — numerical-method policy that deserves a typed field (a real, documented choice) — 6 names + RunPolicy

| knob | read at | default | what it sets | moves numbers? | proposed typed home |
|---|---|---|---|---|---|
| `GPW_DENSITY_EPS` | `Gaussian/Lattice/LatticeScreener.C:45` | 1e-10 | collocated-density screening tolerance (`CollocationEps()`) | yes (A/B 1.2e-10 at 1e-10 → 1.4e-12 at 1e-12) | `Screening::densityEps` |
| `GPW_SCREEN_EPS` | `Evaluator.C:107` | 1e-10 | analytic 1E lattice-sum magnitude screen (CP2K EPS_PGF_ORB analog) | yes | `Screening::overlapEps` |
| `GPW_BECKE_EPS` | `UnitCell.C:127` | 1e-6 | Becke competitor-series tolerance (gate: equivalent-site shares equal to 1e-8) | yes (~1e-6 relative weights) | `Grids::beckeEps` |
| `GPW_VLOC_EPS` | `GPW_Evaluator.C:1287` | 1e-5 | local-PP long-range G-ball tolerance | yes | `Screening::vlocEps` |
| `GPW_FIELDSHARP` | `Evaluator.C:476` | 2/3 | κ: KS-field core exponent as a fraction of α_max (pinned by ρ_lost/N, not wall-clock) | yes | `Grids::fieldSharp` |
| `GPW_RELFIELDSHARP` | `GPW_Evaluator.C:277` | 1/3 | the HartreeOnly floor fraction of α_max | yes (β=0 diverges, +904 Ha) | `Grids::hartreeFieldSharp` |
| **`RunPolicy` ×9** | `Common/Imp/RunPolicy.C:33-57` | see file | `CP2K_COMPAT` (umbrella) and 8 route switches: `QCHEM_DM_LOWRANK`, `GPW_STREAM_FOLD`, `QCHEM_MIX_RHO_M`, `GPW_XC_DM_SOURCE`, `QCHEM_IMPOSE_SYMMETRY`, `QCHEM_BECKE_XC`, `GPW_DAWARE_SCREEN`, `QCHEM_U_EIGEN` | yes (each a declared deviation from CP2K) | **already the model**: typed `Deviation{knob, what, cp2kValue, value, stated}` with provenance and banner.  Becomes the registry's first client; rename the `GPW_*` names that are not about basis sets (`GPW_STREAM_FOLD`, `GPW_XC_DM_SOURCE`, `GPW_DAWARE_SCREEN`) with an alias period |

The six tolerances are the "very technical advanced options": sensible defaults, almost never touched, but a CP2K-parity
run or a convergence study sets them, so they must be reproducible from a deck, not from a shell.  `GPW_FIELDSHARP` and
`GPW_RELFIELDSHARP` already have typed twins (`RasterFields` / `GPWParams::relFieldSharp` arguments flow through the
evaluator) — they are read a second time from the environment as an override.

## 4. Tier 3 — resources, and the two places the environment is already a deck

**Resources** (never change a number; settled by D-THREADS): `QCHEM_OPENMP_THREADS` (alias `GPW_OMP_THREADS`),
`QCHEM_BLAS_THREADS`; read only for the banner: `OMP_NUM_THREADS`, `KMP_BLOCKTIME`.  **`GPW_COLLOC_MEMO`**
(`GPW_Evaluator.C:506`, density replay-cache depth, default 4) is also a resource/speed knob (memory vs time), no numbers —
it belongs here, not in §2.

**The environment as an input deck (the thing the deck concept replaces):**

* **`IntegrationTests/GPW/Harness.C` `EnvOverrides`** maps 14 environment variables onto typed fields before a test runs:
  `GPW_MEASURE`→`SCFParams::Δρmeasure`, `GPW_EPS`→`MinΔρ`, `GPW_NMAX`, `GPW_PULAY`, `GPW_PULAY_START`, `GPW_MOM`,
  `GPW_ACC`→accelerator, `GPW_KERKER_G0`, `GPW_IMPOSE`, `GPW_ORTHO`, `GPW_REAL`, `GPW_SEED`, `GPW_SMEAR`, `GPW_VERBOSE`.
  Every one of them is a T1/T2 field of `SolidCalcOptions`/`SCFParams` — i.e. exactly the **table of fields the deck
  serializes**; `EnvOverrides` is the deck loader in embryo, written as an env reader.  (`scripts/retake5a` and the
  Benchmark protocol drive runs through it.)
* **`CLIapps/gpwprobe.C`** reads ~57 `<P>_*` variables (`MNO_*`, `NIO_*`, …; `Env`/`Envd` helpers): the geometry
  discriminators (`_SWAP_SUBLATTICE`, `_SWAP_ORDER`, `_SHIFT`, `_KMESH`), the recipe (`_ORTHO_TOL`, `_CUTOFF_FACTOR`, `_ECUT`,
  `_SHARED_MU`, `_MOM_SEED`, `_REAL`, `_IMPOSE`, `_XC_UNIFORM`, `_XC_ECUT`, `_VET`, `_NR`, `_L`, `_ALPHA`, `_KERKER_G0`,
  `_XC_CUSP`, `_PULAY`, `_PULAY_START`, `_MOM*`, `_KT`, `_EPS`, `_MEASURE`), +U (`_U`, `_ACBN0`, `_U_RADIAL`, …), the response
  runs (`_CHI0`, `_CHI`, `_CHI_PERTURB`, `_CHI_FD`, `_CHI_RESTART`, …), the schedule (`_ANNEAL`, `_ACC`, `_ANNEAL_PENALTY`),
  the arms (`_SKIP_AFM`, `_SKIP_FM`) and the state (`_SAVE`, `_RESTART`).  Almost all are T1/T2 fields; the discriminators are
  *structure edits* and belong in the deck's structure section.  The `.cmd` files kept beside some logs are people
  reconstructing this by hand.

## 5. Findings and the decisions for the user

1. **No T1 choice lives in the library's environment.**  Basis family, PP vs all-electron, functional, +U are already typed
   (`SolidCalcOptions`, `CalcOptions`).  The environment holds T2 tolerances, T4 hatches/diagnostics, and T3 resources —
   and the *tests/probes* use the environment as a T1/T2 deck.  So D-ENV is less "pull options out of the environment"
   than "give the harness's env-deck a real, versioned home".
2. **Cheapest wins (no design needed):** delete the three env twins of typed options (`GPW_XCGRID_NOSELECT`,
   `GPW_RASTER_POLICY`, `GPW_BECKE_ANG`); delete or test-constructor the test-only valve `QCHEM_SPINBLIND_KERKER`;
   fold the `*TRACE` flags into `SCFParams::Verbose` levels.
3. **Decisions I need from you:**
   a. Are the six §3 tolerances to be typed fields (`Screening`/`Grids`) now, or left env-only until a parity campaign needs
      them?  (They ARE what makes a CP2K-parity run reproducible.)
   b. `GPW_MGRID_ECUTS` / `GPW_RELCUTOFF` (CP2K-matching): keep as instruments, or promote to typed `Grids` fields?
   c. Rename the non-basis `GPW_*` names in `RunPolicy` (`GPW_STREAM_FOLD`, `GPW_XC_DM_SOURCE`, `GPW_DAWARE_SCREEN`) to
      `QCHEM_*` with an alias period, as was done for the thread knob?
   d. `QCHEM_SITE_MOMENTS`: promote to the run report (an integrated observable) and drop the env flag?
   e. Is the 22-name diagnostics list to be one `qchem.Diagnostics` reader, or left as-is?
4. **Not inventoried here:** the C++ `#ifdef`/`-D` compile-time options (`QCHEM_OPENMP`, `QCHEM_RELCHECKED`,
   `QCHEM_ARCH_EXPERIMENT`, `-DCP2K_USE_LIBXC`, data-path macros); those are build configuration, not run options, and are
   correctly in CMake.

## 6. DECISIONS (user, 2026-10-03) and the resulting work order

| | decision | consequence |
|---|---|---|
| a | The six §3 tolerances become **typed fields NOW** | `Screening{densityEps, overlapEps, vlocEps}` + `Grids{beckeEps, fieldSharp, hartreeFieldSharp}`, constructor-injected from the facade (they move numbers, so they are NOT global); defaults = today's values ⇒ no anchor moves; env reads become a registry override layer |
| b | CP2K/QE/VASP/ABINIT **parity runs stay routine for a long time**; the user does not know these options well | `GPW_MGRID_ECUTS` and `GPW_RELCUTOFF` are therefore **promoted to typed `Grids` fields** (`Grids::ladder`, `Grids::relCutoff`) with the CP2K meaning written in the doc string and a one-line "what this is for" in the registry; the parity recipe becomes a named deck/preset (`parity-cp2k`) rather than a bundle of env vars |
| c | Take the **GPW prefix off every non-GPW-specific setting** | rename to `QCHEM_*`, old name kept as a deprecated alias (as `GPW_OMP_THREADS` was): at least `GPW_STREAM_FOLD`, `GPW_XC_DM_SOURCE`, `GPW_DAWARE_SCREEN`, and every §3/§2 tolerance or hatch that is not about the Gaussian-plane-wave basis itself (`GPW_DENSITY_EPS`, `GPW_SCREEN_EPS`, `GPW_BECKE_EPS`, `GPW_VLOC_EPS`, `GPW_CONTRACT_CUBE`, `GPW_EXP_RECURRENCE`, `GPW_SPHERE_SCREEN`, ...); the registry carries the alias |
| d | **Promote `QCHEM_SITE_MOMENTS`** | integrated site moments become a standard section of the run report (an integrated observable, `doc/` rule), the env flag is retired |
| e | **One diagnostics module** | `qchem.Diagnostics`: a process-wide, read-only-after-startup registry (a global, deliberately -- diagnostics gate OUTPUT only, never a computation; it writes through `CurrentReport`), each id registered with a one-line description, typo warning, `--list`, scoped RAII override for tests; excluded from the deck, recorded in the output header; the `*TRACE`/verbose ones fold into `SCFParams::Verbose` levels |

Order (cheap and no-anchor first): (1) delete the three env twins + the test-only valve (§5.2); (2) `qchem.Diagnostics` and
migrate the 22 diagnostics; (3) site moments into the report; (4) the `GPW_` → `QCHEM_` rename with aliases; (5) the typed
`Screening`/`Grids` fields, facade injection, registry override layer; (6) the deck (D-ENV step 3).

## 7. Corrections found while doing steps 1-2 (2026-10-03)

* **The first grep missed knobs read through a wrapper.**  `Mesh/XCPolicy.C` `BeckeXCParams` reads `GPW_BECKE_NR`,
  `GPW_BECKE_ALPHA`, `GPW_BECKE_L` and `GPW_BECKE_ROT` through `envi`/`envd` lambdas (a `getenv(n)` on a variable).  They
  are **tier 2** grid parameters whose typed twin is the `BeckeXCParams(nRadial, mhlAlpha, angularDegree)` arguments (the
  `-1` sentinel means "read the env") and `MeshParams::angRot`: they belong with the `Grids` fields of step 5.  (RunPolicy's
  `IsSet`/`AsBool` wrappers were already counted.)  Inventory total: **45 distinct names** in the library + the 9 RunPolicy.
* **`GPW_LOCALPP_RELCUTOFF` is DUAL-USE**: `GPW_Evaluator.C:90` reads it as the NUMERIC local-PP pair→level κ (default 30 Ha;
  the self-convergence check sets 60) — a **tier-2** knob, `Grids::localPPRelCutoff` in step 5 — while `:1238` used the same
  variable as a timing switch.  Split: the numeric read stays an environment read until step 5; the timing switch is now the
  diagnostic `localpp_timing` with NO legacy alias (so setting κ no longer prints timings — the one deliberate behaviour change
  of step 2).
* **Step 1 done**: `GPW_XCGRID_NOSELECT`, `GPW_RASTER_POLICY`, `GPW_BECKE_ANG` deleted (each shadowed a typed option and could
  disagree with the banner); `QCHEM_SPINBLIND_KERKER` replaced by the test-only `KerkerParams::spinBlind` field.
* **Step 2 done**: `qchem.Diagnostics` (qcCommon): `QCHEM_DIAGNOSTICS=id[=value],...`, `=list`, typo warning, the old names as
  legacy aliases (ON unless `"0"`; the value is the argument, e.g. `GPW_MESH_ORTHO=4`), `Scoped` for tests; 22 diagnostics
  migrated (21 + `localpp_timing`).  `GPW_DM_RANK=1` and `QCHEM_DIAGNOSTICS=dm_rank` both verified on NaF (34 report lines; none
  without).  Caveat: sites that cache `static const bool on=Enabled(...)` read it once, so `Scoped` cannot toggle them mid-process.
