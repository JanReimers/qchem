# The input deck — how a run is described, run, recorded and reproduced

**Read this first if you are about to run, vary or reproduce a periodic (solid) calculation.**  Since 2026-10-04 (D-ENV, record `Records/EnvKnobInventory.md`) a run is ONE JSON file.  The environment no longer sets any value that
can move a number: no `GPW_*`/`QCHEM_*`/`MNO_*`/`NIO_*` recipe variables, no `gpwprobe` knob soup, no test-harness `EnvOverrides`.

## 1. Quick start

```bash
cd build/Release && ninja rundeck                                  # the one CLI
../../scripts/vet-decks                                            # every deck under decks/ parses, resolves and builds its lattice + basis (no SCF; seconds)
./CLIapps/rundeck ../../decks/MnO_AFM2_free.json --check           # vet one deck
./CLIapps/rundeck ../../decks/MnO_AFM2_free.json \
    --set solid.hubbard.0.U_eV=4 --set scf.NMaxIter=300 \
    --out $(../../scripts/rundir qchem6 MnO)                       # run it, with two intentional changes; the record goes where runs live (D-RUNDATA)
```

## 2. The rules (user rulings, 2026-10-03/04)

1. **The deck is the only way to set a calculation choice or numerical-method value.**  There is no override-by-environment, no precedence rule.  A retired variable that is still set is reported once on stderr, **ignored**, and recorded in the run's
   `provenance.ignoredEnvironment`.  What stays in the environment cannot change a number: threads (`QCHEM_OPENMP_THREADS`, `QCHEM_BLAS_THREADS`), `GPW_COLLOC_MEMO`, `QCHEM_DIAGNOSTICS=...`.  `scripts/audit-getenv` (ctest `GetenvAudit`) enforces it.
2. **The structure is just a NAME** — a key of `src/Structure/Data/materials.json` (`molecules.json` is not driven by decks yet).  Lattice, atoms, spin decoration and the pseudopotential vocabulary come from it and are never restated.  A supercell is a
   materials entry (`Si_diamond_2x2x2`), not a deck key.  Nothing FM/AFM-specific lives in the framework: an "arm" is built by the integration test that wants it.
3. **The file you wrote is never the record of what ran.**  Every run writes `<outDir>/<structure>.rNNN.json` BEFORE the SCF (so a crash still leaves it): the deck + every `--set` + every default filled in + a header + `provenance`
   (command line, input deck path + checksum, the `--set` strings, code version, the resolved CP2K-deviation table, ignored env vars, `restartedFrom`), and after the run `results` (converged, energy, the `postSCF` results).  `rNNN` is the next free number,
   claimed atomically (parallel jobs cannot collide; a number is never reused).  **To reproduce:** `rundeck X.r003.json` (no `--set`).  **To reproduce with one intentional change:** `rundeck X.r003.json --set key=value` → writes `X.r004.json`.
4. **Units:** energies a person quotes are **eV in the file** (`U_eV`, `Uirrep_eV`, `alpha_eV`); **atomic units in RAM**.  The writer picks the eV value that converts back to the same double, so a record re-runs bit-identically.
5. **A reader is strict**: an unknown key throws, naming the path and the legal keys (a typo must not silently run the default); a missing key keeps the default (so an old deck survives a new field).  A key beginning with `_` is a comment (`_doc`), ignored and not kept.
6. **Important gates are hard-coded integration tests, never a probe command line.**  A gate is a PASS/FAIL check on a result; the framework knows nothing about gates.  Tests call `qchem::deck::Run(spec, prov, dir)` (see `IntegrationTests/GPW/MnO.C`, `Si.C`).

## 3. The schema (see `src/Calculation/Deck.C` for the C++ types and `Imp/Deck.C` for the readers)

```jsonc
{
  "structure": "NiO_AFM2",                // REQUIRED: a name in materials.json
  "kmesh": [2,2,2],
  "basis": { "data": "VALENCE_LOWQ_VA",   // a BasisSetData name: SIPP_SR, VALENCE_LOWQ_{SR,SR2,SPH,VA,VB}, DZVP, ...
             "spherical": true,           // the spherical lattice view (contaminant-free d/f); refused with VALENCE_LOWQ_SR on a transition-metal cell
             "trim": [{"Z":28,"l":0,"alpha":0.06}],   // a STATED shell trim (exclusive with vet)
             "vet": true },                // the pin-22 vet-stage trim at solid.orthoTol
  "solid": {                              // SolidCalcOptions -- every field; Nelec/species/multiplicity default from the material
    "multiplicity": 1, "seed": "IonicSAD", "ortho": "CholeskyPivoted", "orthoTol": 1e-4, "accelerator": "Null",
    "imposeSymmetry": false, "forceComplex": false, "kShift": [0,0,0], "densityEcut": -1, "cutoffFactor": 2,
    "xcMesh": { "cellKind": "Auto|Uniform|Becke", "nRadial": 40, "angularDegree": 29, "beckeEps": 1e-6, ... },
    "tolerances": { "screenEps": 1e-10, "densityEps": 1e-10, "vlocEps": 1e-5, "fieldSharp": 0.667, "relCutoff": 0, "mgridEcuts": [], ... },
    "policy": { "cp2kCompat": false, "streamFold": true, ... },     // the declared CP2K deviations; ONLY stated routes are recorded
    "hubbard": [ { "site":0, "l":2, "U_eV":4.0, "atomicRadial":true, "orthoAtomic":true, "Uirrep_eV":[...], "alpha_eV":0 } ]
  },
  "scf": { "NMaxIter":200, "minDeltaRho":1e-6, "deltaRhoMeasure":"MixerResidual|MaxDeltaD", "pulayDepth":8, "pulayStart":5,
           "kerkerG0":1.0, "useMOM":false, "smearingkT":5e-3, ... },     // SCFParams, one stage ...
  "schedule": [ {"accelerator":"DIIS","scf":{...}}, {"accelerator":"GDM","scf":{...}} ],   // ... OR an annealed recipe (never both)
  "state":   { "save": "auto", "restartFrom": "NiO_AFM2.r003" },       // CK-1: auto = <outDir>/states/<stem>.h5; restartFrom = a revision stem (lineage) or a path; runs the FINAL stage only
  "postSCF": [                            // post-convergence actions on the converged calculation, IN ORDER; each reports itself; results go in the revision
    {"estimateHubbardU": {}},                                            // ACBN0 (U,J) from the last iterate
    {"hubbardLoop": {"maxOuter":8, "tolU_eV":1e-3}},                     // ACBN0 outer loop
    {"independentResponse": {"nq":2}},                                   // chi0 on an nq^3 q-mesh (nq must divide every kmesh division)
    {"hubbardLinearResponse": {"perturb":[0,1], "tol":1e-8, "maxIter":200, "restart":40}},   // the self-consistent q=0 chi (needs forceComplex)
    {"hubbardFiniteDifference": {"perturb":[0], "alpha_eV":0.027}}       // the FD cross-check of the above
  ]
}
```
Pre-flight (`Resolve`, before any revision is claimed) refuses a deck that cannot work: no Hubbard manifold for a `postSCF`, a linear response without `forceComplex`, an `nq` that does not divide the k-mesh, a `perturb` index past the manifolds, `vet`+`trim`, ... .
It deliberately does NOT refuse "imposeSymmetry + a response": whether an imposed group reduces a k-mesh depends on the group (NiO AFM-II 2×2×2: 8 → 8); the response throws at run time if the mesh was reduced (`--check` prints a note).

## 4. Where things are

| what | where |
|---|---|
| the deck types, readers/writers, `Run`, `Vet`, `--set`, revisions | `src/Calculation/Deck.C`, `Imp/Deck.C` (module `qchem.Deck`); tests `src/Calculation/tests/Deck.C` |
| the CLI | `CLIapps/rundeck.C` |
| shipped decks | `decks/` (recipes), `decks/bench/` (CP2K-parity benchmark rows; `scripts/retake5a` runs them), `decks/ladder/`, `decks/probes/`, `decks/campaigns/{NiO,MnO,MnO2,LiMn2O4}/` (see §5) |
| the retired-variable table (name → deck key) | `src/Common/Imp/Environment.C` `RetiredEnvironmentSet`; the settings reference page `src/Common/Settings.dox` (Doxygen) |
| what each old knob became | `Records/EnvKnobInventory.md` §6–§11 (esp. §10.2 mapping table, §11 the campaigns) |
| supercell structures | `materials.json` (`Si_diamond_2x1x1|2x2x1|2x2x2`) |
| the `scripts/` | `vet-decks`, `audit-getenv`, `audit-settings-doc`, `log2deck.py`, `retake5a`, `rundir` (run-data home) |

## 5. The campaign decks (`decks/campaigns/`)

Thirty decks, one per log in `~/Code/qchem6-runs`.  **Two are exact** translations of saved `.cmd` files (`mno_afm2_free_U0_{radEvery,ortho}_chi_20260930.json`).  **The other 28 are DRAFTS** rebuilt by `scripts/log2deck.py` from what each log printed
about itself; the deck's `_doc` lists what the banner did not carry and was ASSUMED (the convergence threshold, 1e-6).  All vet.  A draft must be judged by the person who ran the original before its numbers are compared with a banked one.

## 6. What is NOT in the deck system (yet)

Molecular (non-periodic) runs (`Calculation`/`AtomCalculation` facades, `molecules.json`) — a later increment.  `gpwprobe becke-ladder` (an instrument that scores quadrature rules on a frozen density).  The `CP2K`/`QE` oracle runs (their own inputs).
The structure-equivariance discriminators (swap sublattice, shift) are to become gtests.

## 7. Adding a setting (the rule that keeps this true)

A value that can move a number goes in `SolidCalcOptions` / `SCFParams` / `RunPolicySpec` (or the matching struct) **and** in `Imp/Deck.C`'s `ToJson`/`FromJson` (the reader's `Done()` will refuse the key until you do), with a default equal to today's behaviour so no anchor moves.
Never read it from the environment: `scripts/audit-getenv` fails the ctest run.  Add a deck test (`src/Calculation/tests/Deck.C`).
