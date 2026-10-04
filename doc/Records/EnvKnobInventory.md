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
* **Step 3 done**: site moments are a standard end-of-run line on every polarized result (`SolidCalculation::Converge`;
  `Becke-partitioned Integral w (rho_up-rho_dn) d3r [e]: 0:+1.0000 1:+1.0000 net=+2.0000` on the O2 triplet), with an explicit
  DEFECT message for a polarized Becke run that has no site blocks and an `n/a` for a uniform mesh.  `QCHEM_SITE_MOMENTS`, its
  per-iteration console line and the sampler's "unavailable" print are gone (the json `scf.siteMoments` and the `m_site` column
  already carried it).  `site_moments` left the diagnostics registry (21 diagnostics + `localpp_timing` => 21).
* **Step 4 done, and a correction to §6(c)**: `qchem.Environment::Env(name, legacy)` (new name wins; old name works with a
  one-time stderr notice).  Applied by the user's criterion -- *not specific to the Gaussian-plane-wave basis*: the Becke XC mesh
  `GPW_BECKE_NR/ALPHA/L/ROT/EPS` → `QCHEM_BECKE_*` and the V_xc feed `GPW_XC_DM_SOURCE/MIX/BOOST` → `QCHEM_XC_DM_*`.  §6(c)'s first
  draft also listed `GPW_STREAM_FOLD` and `GPW_DAWARE_SCREEN`; they ARE collocation machinery of the GPW basis (orbit fold of the
  collocation pair streams; the collocation box tolerance), as are the tolerances (`GPW_DENSITY_EPS`, `GPW_SCREEN_EPS`,
  `GPW_VLOC_EPS`), the contraction/recurrence/sphere switches, `GPW_MGRID_ECUTS`, `GPW_RELCUTOFF`, `GPW_COLLOC_MEMO`, the field
  sharpness pair and `GPW_LOCALPP_RELCUTOFF` -- so they KEEP the `GPW_` prefix.  The diagnostics' old names stay as legacy aliases.

## 8. Step 5 findings (2026-10-03): WHERE each tolerance can be injected

**Done (5a):** `MeshParams::beckeEps` (default 1e-6, folded into `MeshParams::ID()` in scientific notation because
`to_string(1e-7)=="0.000000"`); `BeckeXCParams` applies the env override (`QCHEM_BECKE_EPS`, alias `GPW_BECKE_EPS`) as the
override layer; `MakePeriodicBeckeMesh` reads `mp.beckeEps`.  The default reproduces the NaF anchor exactly (−24.43040329);
`QCHEM_BECKE_EPS=1e-5` still changes the mesh (16752 → 16416 points).  The four `QCHEM_BECKE_NR/ALPHA/L/ROT` already had typed
twins (the `BeckeXCParams` arguments, `MeshParams::angRot`).

**The remaining tolerances split by WHO OWNS THE CODE THAT USES THEM, not by tier:**

| owner | knobs | how it is built today | injection point |
|---|---|---|---|
| `GPW_Evaluator` (`BasisSet/Gaussian/Lattice`) | `GPW_VLOC_EPS`, `GPW_LOCALPP_RELCUTOFF` (κ), `GPW_RELFIELDSHARP`, `GPW_MGRID_ECUTS` | one per k-block, built by `GPW_IBS` from `GPWParams` via `GPWFactory` (10 positional ctor args) | **straightforward**: a `GPWTolerances` value in `GPWParams`, forwarded through `GPWFactory` → `GPW_IBS` → `GPW_Evaluator` |
| `NR_Evaluator` (`PG_Cart_MnD`, the molecular pair-loop engine) | `GPW_SCREEN_EPS` (`kScreenEps()`, 13 sites), `GPW_FIELDSHARP` (`kFieldSharp`), `GPW_RELCUTOFF` (`kEnvRelCutoff`) | constructed INSIDE the molecular `Real_BS` that the CALLER builds with `Gaussian::Factory(BasisSetData, cell, engine, angular)` and hands to `GPW_Evaluator` as `itsMol` | **needs a decision** (below): the tolerance is fixed before any GPW/solid option exists |
| `LatticeScreener` (`CollocationEps()`) | `GPW_DENSITY_EPS` | the screener is built "at it"; 5 + 17 (`kDensityEps`) sites | the screener already has a seam (`GeometryOnlyScreener(eps)`, the `DAware` screener): the eps should be a CONSTRUCTOR argument there, supplied from `GPWTolerances` |

**The decision for the NR_Evaluator knobs** (the user's rule: *the facade hands each object the slice it needs at construction*):

* **A. Put the slice in `Gaussian::Factory`'s signature** (`Factory(data, &cell, engine, angular, tol)`), and have `SolidCalcOptions`
  carry the same struct so the CALLER passes it twice.  Honest about where the object is built, but the options exist in two
  places and can disagree.
* **B. `GPW_Evaluator` re-derives its pair-loop evaluator with the tolerances** (a `WithTolerances(tol)` clone of the orbital block,
  the molecular basis staying the caller's).  One source of truth (`SolidCalcOptions`), at the cost of a clone operation on the
  evaluator face and care that `PGData`'s lazy caches are not shared across tolerance settings.
* **C. Leave the three `PG_Cart_MnD` knobs as instruments** (A/B hatches, environment only) and type only what the GPW layer
  and the screener own.  Cheapest; the screening epsilon that CP2K-parity runs most want to set (`GPW_SCREEN_EPS`) stays in
  the environment, which is exactly what the deck (step 6) is meant to remove.

My recommendation: **B**, because the deck needs ONE place to write `screenEps`, and a parity run must be reproducible from the
deck alone.  It is the larger change (the evaluator face, the cache-sharing audit), so it should be decided before it is started.

**DECIDED (user, 2026-10-03): option B** -- `GPW_Evaluator` re-derives its pair-loop evaluator with the tolerances applied
(`SolidCalcOptions` is the one source of truth; the deck can write `screenEps` once).  Work order for step 5: (i) `GPWTolerances`
in `GPWParams` for the `GPW_Evaluator` knobs (`vlocEps`, `localPPRelCutoff`, `relFieldSharp`, `mgridEcuts`) + the screener's
constructor eps (`densityEps`); (ii) option B for `screenEps`, `fieldSharp`, `relCutoff` (the evaluator-face clone and the
`PGData` cache-sharing audit); (iii) the facade carries the struct in `SolidCalcOptions`; (iv) the env stays only as the
override layer applied where the options are built.

**Step 5 (i) DONE 2026-10-03**: `GPWTolerances` (`qchem.BasisSet.Gaussian.Lattice.GPWTolerances`: `vlocEps`, `localPPRelCutoff`,
`relFieldSharp`, `mgridEcuts`) on `SolidCalcOptions::tolerances` -> `GPWParams::tol` -> `GPW_IBS` -> `GPW_Evaluator::itsTol`; the
four environment reads inside the evaluator are gone; `ApplyEnvOverrides` is the ONE place the environment is applied (the facade,
at option resolution) and names what it changed on the run banner (`GPW tolerances (non-default): ...  [environment overrides: ...]`);
a non-default tolerance enters `IDFragment()` (a different basis for caching).  Verified on Si: default silent; `GPW_VLOC_EPS=1e-6`
moves Etot by 2e-8 Ha; ctest 1007/1007.  REMAINING for step 5: option B for `screenEps`/`fieldSharp`/`relCutoff` (the
`PG_Cart_MnD` evaluator) and `densityEps` (the lattice screener's constructor argument).

**Step 5 (ii) DONE 2026-10-03**: `screenEps`, `fieldSharp`, `relCutoff`, `densityEps` joined `GPWTolerances` (env overrides `GPW_SCREEN_EPS`/`GPW_FIELDSHARP`/`GPW_RELCUTOFF`/`GPW_DENSITY_EPS`, applied only in `ApplyEnvOverrides`). Option B as built: `GaussianSharpness::ApplyTolerances(const GPWTolerances&) const`, called by `GPW_Evaluator` right after the cross-cast; `NR_Evaluator` holds the values as instance state and the `static` env reads (incl. `CollocationEps()`, now a constexpr default 1e-10 for tests) are gone. **Deviation from "clone":** it MUTATES the caller-built basis through the const face (a clone of the virtual-diamond IBS stack was not worth it); safe because all GPW_Evaluators of a run get the same options. **Cache audit:** tolerance-dependent = `itsBoxTasks`, `itsBoxTaskOrder`, `itsIntegrateMemos`, `itsFieldHistory` (all cleared on a change); geometry-only (shells, reaches) stay. Unit tests: `GPWTolerances.*` (override layer; `ApplyTolerancesReachesThePairLoopsAndIsReversible`). ctest 1008/1008. REMAINING in step 5: none in code; step 6 (the input deck).

## 9. Step 6 -- the input deck, and the RETIREMENT of the env hooks (planned 2026-10-03, not started)

**Principle (user, 2026-10-03): once the deck exists, the environment must not be a second way to set a T1/T2 value.**  No
precedence rules, no "env overrides deck", no conflicts.  The deck is the single source; the env hooks are removed, not demoted.

**The deck.**  One serialized `RunSpec` (JSON): header (git hash + dirty flag, checksums of the basis / materials / seed data
files, deck schema version) + the T1 and T2 fields (`SolidCalcOptions`, `SCFParams`, `MeshParams`, `GPWTolerances`, the +U block,
the structure).  Defaults are filled in, so the RESOLVED deck is complete.  Every run writes its resolved deck beside its output
(`~/Code/<app>-runs/<Material>/<name>.json`, D-RUNDATA); a result re-runs from its embedded deck and WARNS on a version or
data-checksum mismatch.  ONE loader serves `gpwprobe --deck`, `scfrun`, the GUI and pybind (flag the binding owner; do not edit
`pybind/`).  Named presets (`parity-cp2k`) are decks too.

**Ad-hoc overrides without the environment, and WHAT IS THE RECORD (user, 2026-10-03).**  `--set key.path=value` on the CLI
(`--set tol.screenEps=1e-8`) is the only override path.  The file the user wrote is NEVER the record of what ran (it may lean on
defaults and omit the overrides).  Every run therefore writes TWO things: the **input deck** is left untouched; the
**`<name>.r<NNN>.json`** (`r001`, `r002`, ... -- a REVISION number, because one name is run many times) is the record: the input deck + every `--set` + every default filled in + the header (git hash and
dirty flag, data-file checksums, schema version) + a **provenance block** (the original command line, the input deck's path and
checksum, and per field whether it came from the deck, `--set` or a default).  *Reproduce* = `--deck X.r003.json` with no
`--set` (immune to later changes of the code's defaults, because they were filled in).  *Reproduce with one intentional change* =
`--deck X.r003.json --set tol.screenEps=1e-8`, which writes `X.r004.json`, so the lineage is a chain of self-contained
files.  The revision is the next free number in the run directory, claimed with an atomic exclusive create (parallel jobs cannot collide; a number is never reused or overwritten); the run's other outputs share the stem (`<name>.r003.log`), so deck and results pair by name; the provenance block records `parent: <name>.r002.json` when the run started from a resolved deck.  Stray environment variables that are still set are listed in the output header as ignored.

**Removal, in this order (each step green on ctest before the next):**

| step | what goes | replaced by |
|---|---|---|
| 6a | `ApplyEnvOverrides(GPWTolerances&)` and its call in `SolidCalculation.C:362`; the `QCHEM_BECKE_*` / `GPW_BECKE_*` reads in `BeckeXCParams`; `GPW_MGRID_ECUTS`, `GPW_SCREEN_EPS`, `GPW_FIELDSHARP`, `GPW_RELCUTOFF`, `GPW_DENSITY_EPS`, `GPW_VLOC_EPS`, `GPW_LOCALPP_RELCUTOFF`, `GPW_RELFIELDSHARP` | the typed fields (already there) written in the deck / `--set` |
| 6b | the `RunPolicy` route switches (`QCHEM_DM_LOWRANK`, `GPW_STREAM_FOLD`, `QCHEM_MIX_RHO_M`, `QCHEM_XC_DM_*`, `QCHEM_IMPOSE_SYMMETRY`, `QCHEM_BECKE_XC`, `GPW_DAWARE_SCREEN`, `QCHEM_U_EIGEN`) and `CP2K_COMPAT` | a deck `policy` block; `CP2K_COMPAT` becomes the `parity-cp2k` preset; the banner's Deviation table is unchanged |
| 6c | the A/B hatches that move numbers (`GPW_CONTRACT_CUBE`, `GPW_EXP_RECURRENCE`, `GPW_SPHERE_SCREEN`, `GPW_LONG_SWEEP`, `GPW_XC_DM_*`) | a deck `advanced` block, or DELETE the ones whose A/B is retired |
| 6d | `IntegrationTests/GPW/Harness.C` `EnvOverrides` (14 vars) and `CLIapps/gpwprobe.C`'s ~57 `<P>_*` vars | `--deck` + `--set`; the probes' geometry discriminators move into the deck's structure section.  Retire the `.cmd` files |
| 6e | the enforcement lint | extend `scripts/audit-settings-doc`: a `getenv("` under `src/`, `CLIapps/` or `IntegrationTests/` outside the resource reader (`QCHEM_OPENMP_THREADS`, `QCHEM_BLAS_THREADS`, `OMP_NUM_THREADS`, `KMP_BLOCKTIME`) and `qchem.Diagnostics` FAILS the ctest.  The settings reference page is then GENERATED from the registry |

**What deliberately stays in the environment:** tier 3 resources (threads, BLAS threads) and tier 4 diagnostics
(`QCHEM_DIAGNOSTICS=...`).  Neither changes a number, so neither can conflict with the deck; both are RECORDED in the output header
(bit-identical reproduction needs the same binary and thread counts).

**Transition (one release):** 6a-6d first make a still-set env var a loud WARNING (`<name> is retired; use --set <path>`) and
IGNORE it; deleting the reads happens at the following commit.  Silently honouring it would reintroduce exactly the conflict this
step exists to remove, so there is no "honour with notice" phase (the aliases of step 4 already had one).

**Order and coupling:** needs the deck schema (and D-STRUCTDATA step 2 for the PP set, V-CALCNET for where the option structs
live) before 6a can land; 6a-6c are small once the loader exists.  Anchors must not move: the defaults in the deck equal
today's values, so a deck-less run is unchanged.

**Step 6.1 DONE 2026-10-03 -- the deck's DATA LAYER** (`qchem.Deck`, `src/Calculation/Deck.C` + `Imp/Deck.C`; ctest +6, 1015/1015).
`ToJson`/`FromJson` for `SolidCalcOptions`, `SCFParams`, `MeshParams`, `GPWTolerances` (and the Hubbard manifolds): a writer emits EVERY field,
a reader keeps the default for a missing key and THROWS on an unknown key (naming the path and the legal keys -- a deck typo must not
run the default); enums by exact name; units are the library's (`U_Ha`, no hidden eV conversion).  `ApplySet("a.b.c=value")` (JSON value or bare
string, array indices), `ClaimRevision` (atomic `O_EXCL`, `<name>.rNNN.json`), `WriteRevision` (deck + header + provenance:
command line, overrides, input-deck checksum, code version, ignored env vars), `LoadDeck` (strips the header, warns on a code-version
mismatch, refuses a newer schema).  NOT serialized: `onIteration` (callback), the Hubbard `siteOps`/`greyOps` (derived).
**Step 6.2 DONE 2026-10-03 -- the structure section is a NAME (user ruling)**: `RunSpec{structure, solid, scf}`; `structure` is the key into
`materials.json` / `molecules.json` and nothing about it is restated (lattice, atoms, spin decoration, species all come from the name; a
stray `lattice` key is refused).  `Resolve(spec)` returns the `Material` and fills `solid.Nelec`/`species` when left to derive (stated values are
honoured: a charged cell, a different PP valence), so the resolved deck records what was used.  Unknown name -> throws listing the known names.
A MOLECULE name is refused for now (molecular decks: next).  Not offered: a lattice-constant override (`Materials::Get(name, a)` for EOS scans) --
add an `a` key only if a scan needs it from a deck; today a scan is a new entry or a `--set`-able key added then.  ctest +2.
**Step 6.3 DONE 2026-10-03 -- a deck RUNS** (`rundeck <deck.json> [--set path=value]... [--out DIR]`, `deck::Run`).  `RunSpec` gained `kmesh`, `basis{data,spherical}`
and `schedule[{accelerator,scf}]` (a deck gives `scf` OR `schedule`, never both).  `Run` resolves, CLAIMS+WRITES `<structure>.rNNN.json` BEFORE the SCF
(a crash still leaves its record), builds lattice + basis + `SolidCalculation`, converges, returns `RunOutcome` (energy only when converged).  Verified:
Si Γ deck energy == the hand-built run **bit-identically** (`EXPECT_DOUBLE_EQ`) and = -7.115063 (CP2K anchor); `--set` typo caught with the legal keys; a run
started from `r001` records `parent: ...r001.json`.  **INTERIM HONESTY**: until 6a, the library still honours env overrides, so `rundeck` lists each
set `GPW_*`/`QCHEM_*` (minus resources/diagnostics) in the record's `provenance.activeEnvironment` and prints a NOTE -- the record is never silently
different from the deck.  Known gaps: shell TRIM and the basis `Reader` choices are not in the deck yet (`MakeBasisLowQ`'s `GPW_BASIS_*`/`_TRIM` env
reads are harness-side, retire in 6d); no molecular deck; no restart (`saveStateTo` is serialized, `Restart` is not driven); `rundeck` output dir is `--out`
(use `scripts/rundir`).  ctest 1019/1019.
**NEXT**: 6a (retire the tolerance env hooks: `ApplyEnvOverrides`, `QCHEM_BECKE_*`), then 6b/6c, then 6d (harness/probe env decks -> decks), 6e (lint).  Earlier plan text: the loader in `gpwprobe --deck/--set` and `SolidCalculation` writing its revision, THEN the 6a-6e retirements.

**Step 6a DONE 2026-10-03 -- the tolerance env hooks are RETIRED** (ctest 1019/1019).  Removed: `ApplyEnvOverrides`, the `GPW_{VLOC_EPS,LOCALPP_RELCUTOFF,
RELFIELDSHARP,MGRID_ECUTS,SCREEN_EPS,FIELDSHARP,RELCUTOFF,DENSITY_EPS}` reads, and the `QCHEM_/GPW_BECKE_{NR,ALPHA,L,ROT,EPS}` reads in `BeckeXCParams`
(a negative argument = the default 40 / 2.0 / 29; the typed `MeshParams` fields are the setting).  `qchem::RetiredEnvironmentSet()` /
`WarnRetiredEnvironment()` (qchem.Environment, called from the `SolidCalculation` constructor and `rundeck`) report a still-set retired variable ONCE on stderr with the
deck key that replaces it and IGNORE it; `rundeck` records it in `provenance.ignoredEnvironment`.  Verified: `GPW_SCREEN_EPS=1e-3` leaves Si Etot at -7.115063428;
`--set solid.tolerances.screenEps=1e-3` moves it to -7.065399 and prints the non-default banner.  The Settings reference page lists them as retired with their deck keys.
Historical records and logs (`doc/Records/*`, `doc/logs`) still quote the old variables as what produced those numbers; that is correct as history and is NOT edited.
**RE-RUN OF `~/Code/qchem6-runs` (user, 2026-10-03: do it when the deck system is done):** surveyed 2026-10-03 -- LiMn2O4, MnO, MnO2, NiO, batch; the files are
`gpwprobe` CAMPAIGN logs (frozen/free chi, ACBN0 outer loop, anneal schedule with penalty, save/restart, shell-trim vet loops, structure-edit discriminators, k222 ladders).
So the re-run needs the deck sections those use that do not exist yet: `response` (chi0/chi/perturb/FD), `hubbardLoop` (ACBN0, tolU), `restart`, `trim`, and the
probes' structure edits (swap sublattice, shift) -- i.e. step 6d is the prerequisite; each log's `.cmd` / header is the recipe to translate.  Plan: after 6d, translate each
campaign to a deck, re-run, and compare to the banked number (a number that moves is an anchor to re-judge against an independent route, never to refresh -- Pins pin 10).

**Step 6b DONE 2026-10-04 -- the CP2K-deviation policy is a typed value, the env hooks are RETIRED** (ctest 1021/1021).  `RunPolicySpec` (`cp2kCompat` +
eight `optional<bool>` routes: `dmLowRank streamFold mixRhoM xcFromDM imposeSymmetry beckeXC dAwareScreen hubbardEigen`) is `SolidCalcOptions::policy` = the deck's
`solid.policy` block.  The facade INSTALLS it (`SetRunPolicy`) as the first act of `BuildBasis`, before any factory consults `theRunPolicy()`; `ReresolveRunPolicy` is
gone (`SetRunPolicy` / `ScopedRunPolicy` replace it in A/B arms).  Unset optional = "not stated" (umbrella / default decides); stated wins over the umbrella (rule kept, tested).
A record stores ONLY stated routes + the umbrella (so a later `--set solid.policy.cp2kCompat=true` still means what it says); the RESOLVED table is
`provenance.policyResolved`.  The banner now prints deck keys (`policy.streamFold=on*(stated)`).  Retired + reported (ignored): `CP2K_COMPAT`, `QCHEM_DM_LOWRANK`,
`GPW_STREAM_FOLD`, `QCHEM_MIX_RHO_M`, `QCHEM_/GPW_XC_DM_SOURCE`, `QCHEM_IMPOSE_SYMMETRY`, `QCHEM_BECKE_XC`, `GPW_DAWARE_SCREEN`, `QCHEM_U_EIGEN`.  Verified: `CP2K_COMPAT=1 rundeck`
runs the DEFAULT policy and says so; `--set solid.policy.cp2kCompat=true --set solid.policy.streamFold=true` gives `symmetry FREE [VETOED by policy]`, uniform XC, stream fold stated.
**The `parity-cp2k` PRESET is just `policy.cp2kCompat=true`** (a named preset file is a deck-include feature for later).  Harness / `gpwprobe` / `scripts/retake5a` use the interim
`GPW_PARITY=1` (harness `EnvOverrides`) until 6d.  `doc/Benchmark.md` §2 carries a dated note translating the old variable names (its historical commands are left as written).
**Known hazard (noted, not fixed):** the policy is process-global, consulted at BUILD time; two live runs with DIFFERENT policies in one process would see the last installed.
Sequential runs (every test and rundeck) are safe.  Remaining 6c: the A/B hatches (`GPW_CONTRACT_CUBE`, `GPW_EXP_RECURRENCE`, `GPW_SPHERE_SCREEN`, `GPW_LONG_SWEEP`, `QCHEM_XC_DM_MIX/BOOST`).

**Step 6c DONE 2026-10-04 -- the A/B hatches and the XC-feed controls** (ctest 1022/1022).  DECISIONS (mine, per "you decide"; each reversible from git):
* `GPW_CONTRACT_CUBE`, `GPW_EXP_RECURRENCE`, `GPW_SPHERE_SCREEN`, `GPW_LONG_SWEEP`: **REMOVED** (retired + reported as "removed", ignored).  Each switched between two routes that agree (~ulp / to the screen) and
  is a VERIFICATION instrument, not a run input, so it does not belong in a deck.  The A/B survives as unit-test hooks (`NR_Evaluator::ContractCubeOverride` (existing), new
  `ExpRecurrenceOverride`, `SphereScreenOverride`; test `GPWTolerances.TheRouteHooksAreUlpLevelAndReset`).  `GPW_LONG_SWEEP`'s kappa-sweep path stays as the non-Gaussian-PP fallback; only the env switch went.
* `QCHEM_XC_DM_MIX`, `QCHEM_XC_DM_BOOST` (and `GPW_` aliases): **PROMOTED** to `policy.xcDMMix` / `policy.xcDMBoost` (they are controls of the `xcFromDM` route and travel with it; they were process-lifetime
  `static` reads before, now per-run).  Deck test covers them.
* `GPW_COLLOC_MEMO`: stays an environment variable -- a RESOURCE (memory against time, never changes a number), documented with the thread knobs; 6e's lint must allow it.
Remaining env reads under `src/` after 6c: resources (`QCHEM_OPENMP_THREADS`+alias, `QCHEM_BLAS_THREADS`, `OMP_NUM_THREADS`, `KMP_BLOCKTIME`, `GPW_COLLOC_MEMO`) and the diagnostics registry.  NEXT: 6d (the harness / `gpwprobe` env-decks -> decks), then 6e (lint).

## 10. Step 6d DRAFT -- the deck sections the campaign runs need (drafted 2026-10-04; USER REVIEWS before any code)

**The unit.**  A deck is ONE run of ONE structure.  Everything in `gpwprobe`'s `RunTMO` (and the harness `EnvOverrides`) that is a *choice about the run* has a home below;
everything else is NOT a run input and leaves the environment by moving to where it belongs (§10.3).

### 10.1 Schema (existing keys unchanged; NEW = marked)

```jsonc
{
  "structure": "NiO_AFM2",                 // a name in materials.json / molecules.json; nothing about it is restated
  "kmesh": [2,2,2],
  "basis": {
    "data": "VALENCE_LOWQ_VA", "spherical": true,
    "trim": [ {"Z":28,"l":0,"alpha":0.06} ]            // NEW: a STATED trim (was <P>_TRIM=Z:l:alpha,...)
    // "vet": true                                      // NEW: the pin-22 vet-stage trim at solid.orthoTol (was <P>_VET=1); exclusive with "trim"
  },
  "solid": {                                // SolidCalcOptions -- unchanged: Nelec/species derived, multiplicity, seed, ortho, orthoTol, cutoffFactor,
    "multiplicity": 1,                      //   densityEcut, imposeSymmetry, greyImposition (IMPOSE=2), spinsShareFermi, momFromSeed, forceComplex, xcMesh{...},
    "imposeSymmetry": false,                //   tolerances, policy, hubbard[...] ...
    "hubbard": [ {"site":0,"l":2,"U_eV":4.0,"atomicRadial":true,"orthoAtomic":true, "Uirrep_eV":[3,4,5]} ]   // NEW SPELLING: U_eV, Uirrep_eV (see §10.4 Q1)
  },
  "schedule": [                             // the anneal: stages, each {accelerator, scf}; replaces <P>_ANNEAL / _ACC / _ANNEAL_PENALTY / the lone scf block
    {"accelerator":"Ladder","scf":{"smearingkT":0.01,"momSmearPenalty":0.1,"stopOnAccelExhausted":true}},
    {"accelerator":"GDM",   "scf":{"smearingkT":0.005}}
  ],
  "state": {                                // NEW (CK-1)
    "save": "auto",                         // "auto" = <outDir>/states/<stem>.h5 after every stage; or a path.  (was <P>_SAVE)
    "restartFrom": "NiO_AFM2.r003"          // a REVISION STEM (-> <outDir>/states/<stem>.h5, so lineage is by name) or a path.  Runs the schedule's FINAL stage only.  (was <P>_RESTART)
  },
  "then": [                                 // NEW: post-convergence actions on the converged calculation, run IN ORDER; each reports itself
    {"estimateHubbardU": {}},                                                      // was <P>_ACBN0=1
    {"hubbardLoop": {"maxOuter": 8, "tolU_eV": 1e-3}},                             // was <P>_ACBN0=n, _ACBN0_TOL; re-converges with the schedule's final stage
    {"independentResponse": {"nq": 2}},                                            // was <P>_CHI0=nq
    {"hubbardLinearResponse": {"perturb":[0,1], "maxIter":200, "restart":30, "tol":1e-8}},   // was <P>_CHI=1, _CHI_PERTURB, _CHI_MAXIT, _CHI_RESTART
    {"hubbardFiniteDifference": {"perturb":[0], "alphaHa": 0.005}}                 // was <P>_CHI_FD
  ]
}
```
Rules: a `then` entry whose prerequisite is unmet THROWS before the SCF starts (`hubbardLinearResponse` needs `solid.forceComplex:true`, a full k-mesh and manifolds; `independentResponse`
needs a k-mesh commensurate with `nq`) -- the old probes found out after the SCF.  Each `then` action's inputs and its result summary go in the revision record.  A run's `state.restartFrom`
revision becomes `provenance.parent` automatically.

### 10.2 Mapping: every variable `RunTMO` / the harness reads, and its deck home
| old variable(s) | deck key | note |
|---|---|---|
| `<P>_KMESH` | `kmesh:[n,n,n]` | |
| `<P>_ORTHO_TOL`, `_CUTOFF_FACTOR`, `_ECUT`, `_SHARED_MU`, `_MOM_SEED`, `_REAL`, `_IMPOSE` | `solid.orthoTol, cutoffFactor, densityEcut, spinsShareFermi, momFromSeed, forceComplex (REAL=0), imposeSymmetry+greyImposition` | `IMPOSE=0/1/2` is two booleans |
| `<P>_XC_UNIFORM`, `_XC_ECUT`, `_NR`, `_L` | `solid.xcMesh.cellKind, eCut, nRadial, angularDegree` | the 2026-09-27 "explicit NR/L pins Becke" rule is `cellKind:"Becke"` stated |
| `<P>_U`, `_U_RADIAL`, `_U_IRREP` | `solid.hubbard[]` entries (`U_eV`, `atomicRadial`, `orthoAtomic`, `Uirrep_eV`) | `every|atomic|ortho|orthofull` are four ways of writing the list; `orthofull` = the 8-manifold list with the U=0 spectators |
| `<P>_ALPHA, _KERKER_G0, _XC_CUSP, _PULAY, _PULAY_START, _MOM, _MOM_START, _MOM_PENALTY, _MOM_HOLD, _KT, _EPS, _MEASURE`, `GPW_<P>_NMAX`, `GPW_<P>_VERBOSE` | `schedule[].scf.{startingRelaxRo, kerkerG0, xcCuspDeficit, pulayDepth, pulayStart, useMOM, momStartIter, momSmearPenalty, momGuard.holePersistence, smearingkT, minDeltaRho, deltaRhoMeasure, NMaxIter, verbose}` | |
| `<P>_ANNEAL`, `_ACC`, `_ANNEAL_PENALTY` | `schedule[]` | one stage per kT; `stopOnAccelExhausted` is `scf.stopOnAccelExhausted` on the non-final stages |
| `<P>_VET`, `_TRIM`, `GPW_SPHERICAL`, `GPW_BASIS_SPAN` | `basis.{vet, trim, spherical, data}` | the "spherical d needs VA/SPH span" default+refusal becomes a deck validation error with the same message |
| `<P>_SAVE`, `_RESTART` | `state.{save, restartFrom}` | the FM arm's ".fm" suffix is gone: an arm is its own deck, its own stem |
| `<P>_ACBN0`, `_ACBN0_TOL`, `_CHI0`, `_CHI`, `_CHI_PERTURB`, `_CHI_MAXIT`, `_CHI_RESTART`, `_CHI_FD` | `then[]` | §10.1 |
| harness `GPW_MEASURE/EPS/NMAX/PULAY/PULAY_START/MOM/ACC/IMPOSE/SMEAR/VERBOSE/REAL/KERKER_G0/SEED/ORTHO` | the same `scf` / `solid` keys | `--set` replaces each |
| harness `GPW_PARITY` | `solid.policy.cp2kCompat` | |

### 10.3 What is NOT a run input, and where it goes (so the env can leave without a deck section)
* **Structure edits** (`<P>_SWAP_SUBLATTICE`, `_SWAP_ORDER`, `_SHIFT`): these are DISCRIMINATORS -- tests that the code is equivariant under relabelling/translation.  They become gtest
  cases (the `GPW_MnO.*` suite) calling the library, no environment.  The deck's structure stays a name.
* **Arms**: `_SKIP_AFM` / `_SKIP_FM` and the FM arm's `afm=false` decoration: an arm is a DECK.  Add `MnO_FM2` / `NiO_FM2` to `materials.json` (same cell, both Mn `spin:+1`) so the FM arm is
  `structure:"MnO_FM2", solid.multiplicity:11`.  The AFM-vs-FM comparison + PASS/FAIL checks (staggered moment, charge, ordering) are a GATE, not a run: they stay in `gpwprobe mno`/`nio` (which then RUN two decks
  and judge them) or move to ITMain.
* **Sweeps and ladders** (`SI_LADDER`, `GPW_KSHIFT`, `NAFGDM_*`, `GATE1_*`, `becke-ladder`): a sweep is a LOOP over decks, not a deck: `for n in ...; do rundeck base.json --set kmesh=[$n,$n,$n]; done`.  The probes that only
  loop retire; ones with checks stay as thin drivers over `deck::Run`.  A supercell is a materials entry (`Si_diamond_2x1x1`) if a campaign needs it -- not a deck key.
* **Instrumentation** (the `m(r)` point probe, `Instrumentation(arm, label)`): output only; stays in the drivers.

### 10.4 Questions for the user (my defaults in bold; each is a one-line change)
1. **Units of U in the deck**: **`U_eV`** (what the literature quotes and what `HubbardU(site,l,eV)` already takes) vs the current `U_Ha`.  Reader converts; the record writes `U_eV` back.  Changing the 6.1 key is a one-line rename + test.
2. **Restart by revision stem** (`"restartFrom":"NiO_AFM2.r003"`) with states under `<outDir>/states/` -- **yes** -- vs bare paths only.
3. **`then` as an ordered list of named actions** -- **yes** -- vs flags on the deck.
4. **FM/AFM as separate materials entries** -- **yes (`MnO_FM2`, `NiO_FM2`)** -- vs a deck key flipping the decoration (which would be a structure edit, against your "just a name" ruling).
5. **Gates stay as drivers** over `deck::Run` -- **yes**; the env-var knobs of the probes die with their sub-commands' migration, one probe at a time.

### 10.5 Build order (each step green, ctest count up)
6d.1 `U_eV` spelling + `basis.trim|vet` + `state` (save/restartFrom/lineage) -> 6d.2 `then[]` actions (+ pre-flight validation) -> 6d.3 `MnO_FM2`/`NiO_FM2` + translate `gpwprobe mno/nio` onto `deck::Run` (env knobs deleted, a
`.cmd`-style deck beside each banked log) -> 6d.4 harness `EnvOverrides` deleted (ITMain tests state their options in code; benchmark scripts use decks + `--set`) -> 6d.5 the remaining probes (ladder, ksweep, naf-smear,
becke-ladder, gate1) as deck loops or drivers -> **D-ENV-RERUN** (the `qchem6-runs` campaigns from decks) -> 6e lint.
