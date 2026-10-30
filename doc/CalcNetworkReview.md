# Calculation network — design review (V-CALCNET, lane B, 2026-10-29)

Born at the top of `doc/` as the DESIGN step of `doc/QueueCampaignPlan.md` §1/§5 (lane B, before C4).  Retires to
`doc/Records/` when its migration (§5) is executed or re-planned.  Reading only: no code changed.
Scope: the facades `AtomCalculation`, `Calculation` (molecules), `SolidCalculation`; their options structs; `Materials`;
`Deck`; `SolidState`; and who builds the basis, the PP terms, the seed and the grids.  Questions for the user are in §6.

---

## 1. The network today (what client code sees, and who decides what)

### 1.1 The three facades

| | `AtomCalculation` | `Calculation` (molecules) | `SolidCalculation` |
|---|---|---|---|
| module | `qchem.AtomCalculation` | `qchem.Calculation` | `qchem.SolidCalculation` |
| ctor input | `Z, charge` + `AtomCalcOptions` + `SCFParams` + `AcceleratorOptions` | `Structure` (deep-copied) + `CalcOptions` [+ `SCFParams`] + `AcceleratorOptions` | `Lattice_3D` + **a caller-built molecular basis** (`shared_ptr<const Real_BS>`) + `SolidCalcOptions` + `SCFParams` \| `vector<SCFStage>` + `SolidAcceleratorOptions` |
| ctor converges? | yes, once | yes (default `SCFParams` unless the D-CALC-DOUBLE-SCF overload is used) | yes (single or annealed schedule) |
| result | `bool Converge`, accessors with no proof | `bool Converge`, accessors with no proof | `Outcome<Converged, SCFFailure>`: the answers live ONLY on the `Converged` proof |
| basis built by | the facade (radial exponent pools) | the facade (`Gaussian::Factory(json{basis-name,...})` [+ `SymmetryAdapt`]) | the CALLER builds the molecular basis; the facade builds the GPW / PW Bloch basis over it |
| basis trim (pin 22) | — | none | **the DECK** (`VetStageTrim` in `deck::Run`); the facade only REPORTS conditioning (`VetGpwConditioning`) |
| PP selected by | `pseudopotential` bool + `valence` (Zion) | `pseudopotential` bool + `ppValence` (one q for every atom) → `GetGTH` default or that q | `species` = `vector<pair<element, q>>`, from the deck/`Material` |
| PP data read by | `GetGTH` inside the Hamiltonian factory | `PPZion`/`PPSpecies` in the facade (`GetGTH`) + the Hamiltonian factory | the Hamiltonian factory (router) + the facade twice more (`GatherSharpness` for α_pp, `PseudoAtomicFunction` for +U projectors) |
| functional | `model` / `xc` | `model` | **hard-coded `"LDA"`** (12 lines with a `"LDA"` literal in `Calculation/Imp` + `ChargeDensity/Imp`) |
| grids | `AtomCalcOptions::mesh` | `CalcOptions::mesh` | the GPW basis ladder (`densityEcut`/`cutoffFactor`/`raster`) + the XC mesh, `Auto` RESOLVED by the facade (`ResolveXCMesh`) |
| seed | CoreGuess | SAD (DFT) / CoreGuess (HF) | `IonicSAD` default, built by the facade (`MakeSeedDensity`), library `atomic_valence_densities.json` keyed (Z, functional, Nelec) — **no q** |
| save / restart | — | — | `SaveState` / `saveStateTo` / static `Restart` → `Outcome<unique_ptr, RestartRefusal>`; HDF5 `qchem.SolidState` with REFUSE / WARM fingerprint |
| deck | — | — (molecular decks are "the next increment") | `qchem.Deck`: `RunSpec` → `Resolve` → `Run` |

The three are SIBLINGS with no common face.  That was deliberate when written (the GUI consumes `Calculation`; atoms
differ), and nothing needed one until now: the geometry optimizer is the first client that must drive molecules AND
solids through one door.

### 1.2 The data around them

```
   materials.json ──► Materials::Get(name) ──► Material{cell, species(element,q)}  ─┐
   molecules.json ──► StructureData (molecular decks: not yet)                      │
                                                                                    ▼
   deck JSON ──► deck::FromJson ──► RunSpec{structure name, kmesh, Basis{data|file, trim|vet, spherical, family},
                                            solid: SolidCalcOptions, scf|schedule, state, postSCF}
                    deck::Resolve: fills solid.species / Nelec from the Material; pre-flights postSCF
                    deck::Run:    Lattice_3D(cell,kmesh) → MakeGaussianBasis(.bsd) → [VetStageTrim] → SolidCalculation(…)
                                  or SolidCalculation::Restart(state)

   SolidCalculation::BuildBasis   (everything the saved-state fingerprint is judged against)
       SetRunPolicy · magnetic decoration · GPWFactory / PW Factory · conditioning REPORT · ResolveXCMesh(GatherSharpness:
       GetGTH per species) · Crystal_EC
   SolidCalculation::BuildRun     (everything that depends on the seed)
       +U site groups & PseudoAtomicFunction (GetGTH) · Hamiltonian::Factory(species, "LDA", xcMesh, vxcFit, hubbard)
       · accelerator · MakeSeedDensity (atomic_valence_densities.json) · SolidSCFIterator · banner
```

### 1.3 What the map shows (the placement defects)

1. **The PP set has no owner.**  It is a `vector<pair<string,int>>` passed by value through `Material` → `RunSpec` →
   `SolidCalcOptions` → `Hamiltonian::Factory`, and each consumer re-reads `gth_potentials.json` on its own
   (`GetGTH` is called from 6 production files outside qcPseudopotential).  The molecular facade expresses the same
   choice differently (a bool and one global q).  The seed library cannot be told q at all; the basis cannot be checked
   against it (D-BASIS-PP).  The functional — which SELECTS the PP entry — is a literal.
2. **The basis decision is split across three places**: the deck chooses the file and runs the vet trim; the solid
   facade builds the Bloch basis and only reports conditioning; the molecular facade has neither trim nor report.
   A `SolidCalculation` built without a deck (every ITMain GPW test) gets no pin-22 trim.
3. **"Setup" and "converge" are fused in the constructors.**  `SolidCalculation` already splits them privately
   (`BuildBasis` | `BuildRun` | `Converge`) because `Restart` needed the gap; the split is right and is the seam §3 uses.
4. **Positions are part of the object's identity.**  Every geometry-dependent object (Bloch basis, overlap, PP terms,
   Ewald, Becke mesh, the +U site groups) is built once from one `Structure`, and the fingerprint REFUSES a restart onto
   moved atoms.  There is no way to say "same run, atoms moved".
5. **Two failure conventions.**  The solid facade is Outcome-based (`Converged` proof, `SCFFailure`); the molecular and
   atom facades return `bool` and serve accessors unconditionally.

---

## 2. The rule this review applies

A run splits into three lifetimes, and each object belongs to exactly one:

| lifetime | what | mutability |
|---|---|---|
| **RUN SETUP** | the CHOICES and everything RESOLVED from them once: species + PP set, functional, basis RECIPE (file, per-species shells after the vet trim), grid parameters (resolved numbers, not `Auto`), XC-mesh kind, spin group + multiplicity, magnetic decoration, imposed group, +U manifolds (site, l, U), k-mesh, SCF recipe, run policy | immutable value; built by ONE explicit phase (pin 11, D-ENV "resolve once at the facade") |
| **GEOMETRY POINT** | positions → the instantiated basis (same recipe, moved centres), S, the Hamiltonian's geometry-keyed terms (V_loc, KB, Ewald, Becke weights, site rotations) | rebuilt per geometry; immutable while it lives |
| **SCF STATE** | the density (and C, ε, f), the accelerator / mixer history, the iterator | per SCF; the DENSITY outlives its geometry point and seeds the next one |

Today `SolidCalculation::Imp` holds all three in one bag; `BuildBasis` mixes setup decisions with geometry-point
construction.  Everything below follows from pulling those apart.

---

## 3. The face the geometry optimizer drives

### 3.1 The face

```cpp
//! What a geometry optimizer (or an MD / NEB / phonon driver) consumes: the energy, and later its gradient, at
//! positions it chooses, many times, over ONE fixed run setup.  Implemented for molecules and for solids.
class PotentialEnergySurface            // abstract, data-free (CLAUDE.md design rule)
{
public:
    virtual ~PotentialEnergySurface() = default;
    virtual const Structure& Reference() const = 0;   //!< species, atom order, (cell) -- the positions' frame
    //! Converge at \a R (Cartesian bohr, Reference()'s atom order; a solid's cell is fixed by the setup), seeded from the
    //! LAST converged density carried onto \a R.  Non-const: it advances the surface's state.  A non-convergent SCF is a
    //! legitimate outcome (Outcome); a request the setup cannot honour (forces from a term with no gradient, R breaking an
    //! imposed group, a basis that left its conditioning window) is an Outcome too -- the optimizer must stop, not crash.
    virtual Outcome<PESPoint, PESFailure> ConvergeAt(const std::vector<rvec3_t>& R, const PESRequest& want) = 0;
};

struct PESRequest { bool forces = false; bool stress = false; };   // stress: the variable-cell seam (later)
struct PESPoint                                                    // a VALUE: the optimizer keeps a history of these
{
    double                              energy;      //!< Ha
    std::optional<std::vector<rvec3_t>> forces;      //!< −dE/dR, Ha/bohr, present iff requested and supported
    size_t                              iterations;  //!< the SCF's, so a driver can see the warm start working
};
//! WHY a point failed.  SCFFailure (existing) + the setup refusals named above, each with its cause.
struct PESFailure { enum class Why { SCF, BasisConditioning, SymmetryBroken, NoGradient } why; std::string details; };
```

Two concrete implementers, each a SESSION over one setup:

```
   MoleculeSetup ─┐                          ┌─► MoleculeSession : PotentialEnergySurface
                  ├─ RunSetup (immutable) ───┤
   SolidSetup  ───┘                          └─► SolidSession    : PotentialEnergySurface
                                                     owns: current GeometryPoint, current SCF state, previous density
```

The single-point facades stay, as thin compositions:  `SolidCalculation(lat, …)` = `SolidSetup` + `SolidSession` +
one `ConvergeAt(reference)`, with today's `Converged` proof kept for the rich single-point answers (densities, +U
loop, responses).  `Calculation`'s public API does not change (the GUI and pybind consume it); it gains the
`Outcome` overload beside the `bool` one.  `AtomCalculation` has no geometry and stays outside this face; it shares
the PP set (§4) and nothing else.

### 3.2 What `ConvergeAt` does, and what it keeps

1. **Check R against the setup.**  Same atom count and species order.  If the setup IMPOSES a group, R must respect it
   to a stated tolerance — else `SymmetryBroken`: the optimizer symmetrizes its forces, so a broken R is the
   optimizer's defect, and lowering the symmetry (Jahn–Teller in Li₂Mn₂O₄) is a NEW setup, chosen by the user
   (question Q4).  The group itself is R1.0r's intersection (crystal ∩ k-mesh), fixed at setup.
2. **Build the geometry point** from the frozen recipe: the SAME shells (pin 22's trim was decided ONCE at the reference
   geometry and is NOT re-run — re-vetting per step makes the PES discontinuous), the same grids (raster N from C2's ε
   rule depends on the exponents, not on positions, so it is geometry-independent by construction — C2 should keep it
   that way), the same XC-mesh kind (a Becke mesh MOVES with the atoms; its weight derivatives are a force term),
   the same PP set.  Then measure λ_min(S(R)): if it left the window the trim was decided in (below `orthoTol` × a
   margin), refuse with `BasisConditioning` rather than silently change the basis (Q3).
3. **Carry the density.**  The AO density follows the atoms (D_ij is indexed by functions that moved with their
   centres), but it is not idempotent in the new metric S(R).  First increment: carry C, re-orthonormalize it in S(R)
   (Löwdin), rebuild D — the seed role is `WorkingDensity`, exactly a `BuildStage` carry today.  Later lever: history
   extrapolation over the last k geometries (CP2K's ASPC); the session owns that history, no other object does.  A PW
   session carries ρ(G) / C(G) directly (no metric change).
4. **Converge** with the setup's restart policy: kT = 0 restarts enter the Ladder at its GDM rung (CLAUDE.md, the NaF
   measurement) — the regime matrix (QueueCampaignPlan P4a) has the row "carried stage, changed Hamiltonian (geometry
   step)", UNMEASURED today; this is its customer.
5. **Forces, if asked** (the seam, §3.3).  Fixed U along the trajectory (no ACBN0 loop per step: Q5).

### 3.3 The forces seam (designed later, by the forces note — lane A, C3)

dE/dR at a converged density is a sum over TERMS plus the basis's Pulay term:

```
F_A = − Σ_terms ∂E_term/∂R_A |_D   − Tr(W ∂S/∂R_A)          W = energy-weighted density (Pulay; zero for PW)
```

The natural carrier is lane A's stateless Hamiltonian: `H.At(cd)` returns an immutable value with everything
density-derived (V_H, v_xc, W, the block matrices).  The seam is one CAPABILITY face on that value, reached by the
sanctioned abstract→abstract cross-cast — `GeometryGradient { std::vector<rvec3_t> dE_dR() const; }` — plus a
basis-side capability for ∂S/∂R, ∂T/∂R.  A term that does not implement it makes `ConvergeAt(forces)` return
`NoGradient` naming the term (a P3 precision contract: refuse with the cause, never a silent zero).  Until then the
optimizer can run on finite differences of `ConvergeAt` energies, which is also the forces' test oracle.

Forces-note inputs this review hands over: the Becke weight derivative (the mesh moves), the uniform-raster egg-box
(§2(b) of the campaign plan; P3 escalates to an alias-free raster for a force consumer), R1.0r (forces are dE/dR only
if ρ is symmetrized under the k-mesh's group), and the Löwdin carry (it is a seed, so it does not enter the
derivative).

---

## 4. Where PP data, basis provenance and the seed variant q are read

### 4.1 The options

| | read where | for | against |
|---|---|---|---|
| **A** (current plan of D-STRUCTDATA step 2) | the Calculation level parses `pseudopotentials.json` and hands `species` down | one place, simple | Calculation-level code then knows a file format; every consumer still re-reads `gth_potentials.json`; molecules stay separate |
| **B** | the Hamiltonian PP terms look up their own entry by element | data next to its main user | a hidden read at will (D-ENV forbids it); the seed, the basis check, the sharpness gather and the +U projector would each read it again and could disagree; the terms would have to know about the material's choice |
| **C** | `BasisSet::Orbital_PP_IBS` | — | it is an integral SERVICE (∫ χ V χ for a species field it is GIVEN); owning the choice of potential would make the basis decide physics.  Its job does not change |
| **D (recommended)** | qcPseudopotential reads; the RUN SETUP chooses; everything else is INJECTED | one reader, one choice, one value every consumer shares; testable; the same value for molecules and solids | one new value type and its plumbing (§5 steps 1–3) |

### 4.2 The recommendation (D) in detail

- **Data**: `src/Pseudopotential/Data/pseudopotentials.json`, beside `gth_potentials.json`: per element × PP family
  the q variants, their default, the basis provenance tag each `.bsd` block must declare, and the seed key.  Materials
  name a set (`"pp": "gth-lda-lowq"` or similar) and stop carrying `species` (D-STRUCTDATA step 2's ruling stands:
  PP data is a property of (element, PP), not of the material).
- **The one reader**: qcPseudopotential exports a value type, e.g. `PseudopotentialSet` — per Z: element, family,
  functional, q, Zion, the GTH parameters, the provenance tag, the seed key.  nlohmann stays inside the module's
  implementation unit.  It answers the questions its clients ask (`Zion(Z)`, `Potential(Z)`, `SeedKey(Z)`,
  `ExpectedBasisTag(Z)`) — named for the consumer (CLAUDE.md), no `GetTable()`.
- **The one choice**: the run SETUP resolves it once from (the material's/molecule's set name, a deck override per
  element, the functional) — the explicit phase.  A molecule without a named set gets each element's default, which
  is what `ppValence = 0` does today.
- **Injected into** (constructor arguments, never read at will): the Hamiltonian factory (replaces the
  `vector<pair<string,int>>` + `"LDA"` pair), the seed (`MakeSeedDensity` / `AtomicDensity` key on (Z, functional,
  q) — the folded "seed keyed by q" row), the basis factory (it CHECKS each block's declared tag and THROWS on a
  mismatch — D-BASIS-PP; a broken invariant no caller can act on, so `throw`, and `deck::Vet` catches it in seconds),
  the XC-mesh sharpness gather, the +U pseudo-atomic projector, `Material::Nelec()`, and the saved-state fingerprint.
- **Answer to the user's suspicion**: yes, the Calculation level is too high to READ the data — it reads nothing in
  D; it only NAMES the set.  But the terms are too low to CHOOSE it: the choice is per run (q7 vs q15 Mn is the same
  crystal), and the seed and basis must agree with the terms, which only a value handed to all of them guarantees.

### 4.3 Basis provenance

The `.bsd` block declares (family, q) in a machine-readable header line; the BasisSet library parses its own files
and exposes the tag only to the check above.  D-BASIS-REGISTRY's file-stem discovery + a content checksum go into the
run record and the fingerprint (the record stores the checksum, not only the path).  The setup keeps the basis RECIPE
— file checksum, per-species shell list after the trim — which is also exactly what the optimizer must freeze.

---

## 5. Migration — ordered, each with its gate

Every step lands on lane B as its own chunk-sized commit series; the gates use the controller's sweep +
`scripts/anchordiff` (CLAUDE.md).  Steps 1–5 are C4; 6–8 are new work that unblocks the forces path.

| # | step | gate | closes |
|---|---|---|---|
| 1 | `pseudopotentials.json` + `PseudopotentialSet` in qcPseudopotential; materials name a set; `Material::Nelec()` from it; `species` leaves `materials.json` | UT: every material resolves and its Nelec equals today's (an invariant over EVERY entry, no count pin); `Deck.EveryShippedDeckVets`; sweep: 0 anchors moved | CleanCode D-STRUCTDATA step 2 (main part) |
| 2 | Inject the set: `Hamiltonian::Factory`, `GatherSharpness`, `PseudoAtomicFunction`, the molecular `PPZion`/`PPSpecies` (the molecular `pseudopotential`/`ppValence` fields become a set name + per-element override, old fields kept as a deprecated alias for the GUI); the functional comes from the setup, not the literal | sweep bit-identical (`anchordiff --tol 1e-10`); a grep-lint: no `GetGTH` outside qcPseudopotential and the basis generator (`ValenceBasisGen`, `valgen`, `gth2upf` legitimately enumerate q variants); no `"LDA"` literal in `Calculation/Imp`, `ChargeDensity/Imp` | part of "functional recorded truthfully" (V); prerequisite for GGA |
| 3 | Seed keyed by q: `AtomicDensity` Get/Has take the key from the set (~9 sites, scoped in `Records/OpenWork_History5.md`) | UT: Li q1 vs q3 IonicSAD picks different entries; the sweep flags any test whose q disagrees with the library (fix the test's species, never the lookup) | the folded "seed keyed by q" row; with `valgen --seed-into` (writes q) |
| 4 | Basis provenance: header tag per `.bsd` block; factory checks it against the set and throws; registry by file stem + checksum in the record | UT: a q3 run on a q1 block throws, naming both; every shipped deck vets | CleanCode D-BASIS-PP, D-BASIS-REGISTRY |
| 5 | Fingerprint from the set + the recipe: PP identity (family, q per Z) and basis checksum REFUSE; functional recorded and WARM; positions WARM (with the max displacement printed) — the D-STATE-SCHEMA additive path first, so existing states still load | SolidState UTs: an old file restarts; an LDA→PBE restart is warm; a moved-atom restart is warm and says by how much | OpenWork "restart fingerprint incomplete" (in part); CK-3's geometry half |
| 6 | Extract `SolidSetup` (immutable): the decisions now in `BuildBasis` + the deck's vet trim move into ONE resolve phase; the D-ENV typed fields and the GPW tolerances (bound where the molecular basis is built — D-ENV-TOL option A) live there; `SolidCalculation` = setup + session + one converge | bit-identical sweep; an ITMain GPW test built WITHOUT a deck now gets the same trim as its deck twin (the defect in §1.3 item 2) | D-ENV-TOL; part of D-ENV step 5; part of C5's actuator (the trim gets one home) |
| 7 | The same for molecules: `MoleculeSetup`; `Calculation` keeps its API and gains the Outcome overload; molecular decks (`structure` from `molecules.json`) | M_* anchors bit-identical; a molecular deck vets and runs | D-CALC-DOUBLE-SCF migration (callers move to the setup); the "molecular decks" increment |
| 8 | `PotentialEnergySurface` + `MoleculeSession` / `SolidSession`, energy only, density carry by Löwdin | (a) `ConvergeAt(reference)` == the single-point energy, bit-identical; (b) `ConvergeAt(R+δ)` from the carry == a fresh run at R+δ to 1e-8 Ha, in fewer iterations; (c) E along a 20-point displacement path has no jump above the SCF tolerance (frozen trim); molecule (H₂O) AND solid (NaF Γ) | opens the optimizer; finite-difference forces available as the oracle |
| 9 | (forces note, lane A / C3) the `GeometryGradient` capability on `H.At(cd)` + the basis ∂S/∂R | analytic vs step-8 FD forces | — |

Not in this migration: AtomCalculation unification (no client needs it), variable cell (the `stress` seam only), the
PW save/restore (its own OpenWork row; a PW session carries in memory regardless).

---

## 6. Questions for the user — RULED 2026-10-29 ("Everything looks good.")

- **Q1 — PP placement: option D.**  qcPseudopotential READS, the run setup CHOOSES, everything else is HANDED the value.  RULED: agreed.
- **Q2 — the face.**  `PotentialEnergySurface::ConvergeAt(R, want)` over a fixed setup, single-point facades recomposed on top,
  `Calculation`'s API kept.  RULED: agreed, name accepted.
- **Q3 — basis leaving its conditioning window mid-relaxation.**  RULED: REFUSE that point (`BasisConditioning`).
- **Q4 — symmetry along a relaxation.**  RULED: keep the setup's imposed group (forces symmetrized, a broken R refused); symmetry
  lowering is a NEW setup the user chooses; a FREE relaxation stays available as a setup choice.
- **Q5 — +U along a relaxation.**  RULED: U held fixed at the setup's values.
- **Q6 — molecular PP fields** (`CalcOptions::pseudopotential` / `ppValence`: deprecated alias, or clean replace).  RULED: **clean
  replace** ("I vote for clean replace") -- a set name + a per-element override, NO deprecated aliases.  pybind is not edited by
  lane B; the exact `pybind/` lines referencing the old fields are listed in the step-2 report for the binding owner.

---

## 7. Steps 6–8 as designed (lane B, 2026-10-30; C4 steps 1–5 landed)

**One deviation from §3.1's drawing, and why.**  §3.1 drew `SolidCalculation` = setup + session + one `ConvergeAt`.
The composition runs the OTHER way here: a SESSION is a sequence of single-point calculations over ONE shared setup.
The three lifetimes of §2 are unchanged -- the setup is a shared immutable value; the single-point object
(`SolidCalculation` / `Calculation`) IS the geometry point plus its SCF state; the session owns only what crosses
geometries (the carried density, later the extrapolation history).  What this buys: the whole rich single-point API
(`Converged` proof, +U loop, responses, FD χ, `SaveState`) stays on the per-point object untouched, and a session
step is "build the next point from the setup, seeded with the carry" -- the same act as `Restart`, from memory
instead of disk (CK-3's geometry half falls out of it).

| type | module | holds | built by |
|---|---|---|---|
| `SolidSetup` | `qchem.SolidSetup` (new; `SolidCalcOptions` + its helpers MOVE here, re-exported by `qchem.SolidCalculation`) | the options as asked; the RESOLVED decisions: `PseudopotentialSet`, imposition (the policy veto applied), magnetic decoration, XC mesh (`Auto` resolved), (nUp, nDn), ortho (PW → Orthonormal), the accelerator options; the molecular basis at the reference geometry; (6b) its recipe | `SolidSetup::Resolve(lat, mol \| recipe, opts, acc)` -- the ONE explicit phase (pin 11); installs the run policy first |
| `GaussianBasisRecipe` (6b) | `qchem.SolidSetup` | source (registered data \| file), the trim -- STATED, or decided ONCE by the vet loop at the reference geometry and then frozen (§3.2 item 2) --, spherical view, the GPW tolerances (D-ENV-TOL option A) | `Make(structure)` checks PP provenance (closes the "Deck path only" gap) and builds; the deck's `MakeGaussianBasis` + trim code move here |
| `MoleculeSetup` (7) | `qchem.MoleculeSetup` (new; `CalcOptions` moves) | PP set or all-electron, the valence structure rule, (nUp, nDn) + the polarized promotion, the basis recipe (name, engine, angular, SALC on/off), model, resolved seed, resolved accelerator | `MoleculeSetup::Resolve(st, opts, acc)`; `Calculation` keeps its API, gains `Attempt(params) -> Outcome<Converged, SCFFailure>` beside `bool Converge` |
| `PotentialEnergySurface` (8) | `qchem.PotentialEnergySurface` (new) | the face of §3.1, verbatim; a wrong atom count is a broken invariant (throws), not a `PESFailure` | -- |
| `SolidSession`, `MoleculeSession` (8) | in their facades' modules | `shared_ptr<const Setup>`, the reference structure, the CURRENT point (`unique_ptr<SolidCalculation>` / `<Calculation>`), the carried density | `ConvergeAt(R)`: check R (count, species order; imposed group → `SymmetryBroken`) → structure at R (clone, set positions) → `Make` the basis from the frozen recipe → λ_min(S(R)) vs the setup's window (`BasisConditioning`, Q3) → a new point seeded `WorkingDensity` from the carry → `PESPoint` |

**The carry, first increment: D, not C.**  The solid side already has the vehicle: `RestoreDensity` rebuilds per-block
D over a basis of the same block shape.  A `SnapshotState(bs, wf)` (the in-memory half of `WriteSolidState`) makes the
step a restart from memory.  D carried verbatim onto S(R) is not idempotent there; the first fill makes it so (a seed,
it does not enter any derivative).  The Löwdin carry of C (§3.2 item 3) is the measured lever after gate (b) shows
the iteration count; the molecular side needs the same per-block D → `tDM_CD<double>` composition, built beside it.

**Increments, each its own commit series and gate** (serial sweeps + `anchordiff --tol 1e-10`):

- **6a** -- `SolidSetup` extracted; `BuildBasis` = `Resolve` + build the geometry point; no API change.  Gate: bit-identical.
- **6b** -- `GaussianBasisRecipe` + ctor overloads that take it; the deck builds through it (its trim code deleted);
  the tolerances bound in the recipe (D-ENV-TOL option A).  Gate: bit-identical; a recipe-built run with `vet` gets the
  deck twin's trim; a q mismatch throws without a deck.
- **7a** -- `MoleculeSetup` + `Attempt`.  Gate: M_* bit-identical.  **7b** -- molecular decks (`structure` from
  `molecules.json`): the `RunSpec` grows a molecular arm.  Gate: a molecular deck vets and runs.
- **8** -- the face + both sessions, energy only.  Gates (a)–(c) of §5 on H₂O and NaF Γ.

**As built (lane B, 2026-10-30)** -- `7d980306` 6a, `301f28fb` 6b, `34dc9205` 7a, `d5835f33` 8a, `2174cab5` 8b:
- 6a/6b as designed.  `SolidCalculation(setup, params | schedule)` and `Restart(path, setup, stage)` are the doors; the old
  constructors resolve and delegate.  `deck::Run` resolves one setup; its trim code is gone.
- 7a as designed, without `Attempt` (its failure type, `SCFFailure`, needs a home both facades can import -- OpenWork row).
- 8: the face; `SolidSession` refuses a broken imposed group (space-group / Shubnikov order), a basis the vet loop would
  trim further at `orthoTol` (Q3), and forces/stress; `MoleculeSession` refuses a change of the SALC abelian group.
  MEASURED (UTs): (a) the reference point is the single point bit for bit on both; (b) a carried step equals a fresh run --
  Si Γ, 0.05 bohr: 11 vs 13 iterations at Δρ/N 1e-8 (energies within the two stops' 1e-8 errors; at 1e-9 Si Γ sits on its
  Δρ quadrature floor, ρ_lost/N 1.8e-9, and neither run is meaningful); H₂O HF/DZVP: 25 vs 34, equal to 1e-12 Ha.
  Gate (c), NaF Γ, the Q3 UT and the Löwdin lever are the OpenWork row "PES / run-setup remainders".
