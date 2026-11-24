# Carrying the run policy on the run: who can hand each reader its policy

Design note for CleanCode row **D-MUTABLE-STATICS-RUNPOLICY** (lane B, 2026-11-23).  It is for the user to read BEFORE any code.
The user's ruling (2026-11-23, quoted in the row): *"For thePolicy I will follow your recommendation"*.  That means carrying
the resolved policy on the run, not in the process-wide `thePolicy()` (`src/Common/Imp/RunPolicy.C`).  This note takes the
inventory, groups the readers by which object can hand them the policy, and names the hard group.  The open work stays in
the row.

## 1. The inventory (main 341ed8e5a)

The row's "37 reads" counted every mention, including comments, imports and tests.  The **production call sites are 14, in 8
files**, and they read 8 routes plus the two XC-from-D controls:

| site | route read | what it decides |
|---|---|---|
| `Calculation/Imp/SolidSetup.C:432` | `SymmetryImposition` | `s->imposed` (the permission ANDed with the caller's flag) |
| `Calculation/Imp/SolidSetup.C:525` | `BeckeXC` | `ResolveXCMesh` (the Becke mesh or the uniform grid) |
| `Calculation/Imp/SolidCalculation.C:295` | `Banner()` | the run banner line |
| `Calculation/Imp/SolidCalculation.C:309` | `HubbardEigen` | the +U line of the report (print only) |
| `Calculation/Imp/SolidCalculation.C:361` | `XCFromDM` | the XC-source line of the report (print only) |
| `Calculation/Imp/SolidCalculation.C:3179` | `BeckeXC` | a mesh decision on the Auto cell kind |
| `BasisSet/Gaussian/Lattice/Imp/GPW_Evaluator.C:151` | `DAwareScreen` | which `LatticeScreener` each Bloch evaluator builds |
| `BasisSet/Lattice/Imp/BasisSet.C:168` | `StreamFold` | arming the orbit fold in the `GPW_BasisSet` constructor |
| `ChargeDensity/Imp/DensityMixer.C:97` | `MixRhoM` | (ρ,m) or (↑,↓) leaves in `ComposePeriodic` |
| `ChargeDensity/Imp/Factory.C:22` | `DMLowRank` | `DefaultRhoRoute()`: factored or direct ρ for every `IrrepCD_Factory` that does not state one |
| `ChargeDensity/Internal/Imp/DensitySampler_Singles.C:134` | `XCFromDM` | `UseDMSource()`, inside `ExactSourceOf`: may XC read ρ[D] from a mixed field |
| `…/DensitySampler_Singles.C:169` | `XCDMMixOverride` | the XC-from-D α override |
| `…/DensitySampler_Singles.C:170` | `XCDMBoost` | the XC-from-D α boost |
| `Hamiltonian/Internal/Imp/Hubbard.C:214` | `HubbardEigen` | `Hubbard_U`'s eigenvalue (Dudarev) or diagonal (CP2K) form |

Two readers ALREADY carry the policy on the run.  Each builds its own `RunPolicy(o.solid.policy)`:
`SolidCalculation.C:1390` (the saved state's fingerprint) and `Deck.C:1634` (the record's `policyResolved`).  They are the
pattern the rest would follow.

**The one writer** is `SolidSetup::Resolve` (`SolidSetup.C:252` and `:262`, `SetRunPolicy(opts.solid.policy)`).  The test side
has the RunPolicy unit test, one read in `IntegrationTests/GPW/Boxes.C:650`, and `ScopedRunPolicy` in `ChargeDensity/tests/KerkerMix.C`.

## 2. The readers, grouped by who can hand them the policy

### A. The setup: trivial (6 sites)
`SolidSetup` holds `opts.solid.policy` at both read sites, and `SolidCalculation` holds its setup.
- **Hand-off:** `SolidSetup` resolves `RunPolicy policy{opts.solid.policy}` once, keeps it as a member, and both files read
  `setup->policy`.  The banner then prints the policy the run HOLDS, which the row asks for.
- **Cost:** one commit, no signature beyond the setup.

### B. The basis: easy, one record (2 sites)
`DAwareScreen` (in the `GPW_Evaluator` constructor) and `StreamFold` (in the `GPW_BasisSet` constructor) are both read while
`GPWFactory` builds the Bloch basis, from a `GPWParams`.
- **Builders of that `GPWParams`:** the setup (`SolidSetup.C:607`) and `SolidCalculation.C:480, 1816, 2433`.  Six test sites
  call `GPWFactory` directly.
- **Hand-off:** two fields on `GPWParams`, defaulted to the qchem values so the direct test calls stay unchanged.
  `GPW_BasisSet` passes them on to each `GPW_IBS` / `GPW_Evaluator` it builds.
- **Precedent:** the `CollocationRoute`, bound per basis through `ToleranceBinding::BindRoute`, retired three process-wide
  overrides this way.

### C. The SCF: easy, one record (1 site)
`MixRhoM` is read in `ComposePeriodic`, which `KerkerMixerFactory` / `PulayMixerFactory` call from
`SCFIterator.C:1657` / `:1660`.
- **Hand-off:** a `rhoM` field on `KerkerParams` (`PulayParams` inherits it), filled from `SCFParams`, beside `XCCuspDeficit`.
  The XC cusp-deficit choice already travels this way (`DensitySampler_Singles.C:150`: *"decided when the MIXER is BUILT … so
  CP2K parity is a property of the constructed object"*).
- **Tests:** `KerkerMix.C`'s `ScopedChannelBasis` becomes a stated field.

### D. The Hamiltonian: moderate (4 sites, two routes)
- **`HubbardEigen`** is read once, in `Hubbard_U`'s constructor.  `Hamiltonians.C:626` builds it from the
  `HubbardManifold` list.
  - **Hand-off:** a form field beside the manifolds (or on each `HubbardManifold`), filled by the facade that builds the U
    list.  Easy.
- **`XCFromDM` and its two controls** are read in free functions in `DensitySampler_Singles.C`.  `ExactSourceOf(cd)` is the
  gate, and it is called PER ITERATION by all three sampler families: Singles, Pair through `HasExactSource` /
  `HasXCCorrection`, and Projected (`DensitySampler_Pair.C:63`, `DensitySampler_Projected.C:32`), plus the gradient sampler's
  refusal (`:793`).  The samplers are built by `MakeDensitySampler` (`Hamiltonians.C:222`, `:618`, and the recursion inside
  `ChargeDensity/Imp/DensitySampler.C`).
  - **Hand-off:** an `XCSourcePolicy {bool fromDM; std::optional<double> mix; double boost;}` given to `MakeDensitySampler`.
    Each sampler holds it, and the free functions take it as an argument.
  - `DensitySampler` is not a tier-1 hub, so this is plumbing, not a design discussion.  But it touches every sampler family
    and both Hamiltonian factory paths.
- **The alternative, rejected here:** let the mixer decide, as it does the cusp correction, by depositing the D source only
  when a route wants it.  That changes what `tDM_Sourced_CD::DMSource()` means.  That face is on the protected
  `ChargeDensity` hub (`scripts/hub-interfaces.snapshot:327`), and `FieldMixer.C:110` (pin 28) relies on the source being there.

### E. The density route: THE HARD ONE (1 site, every run kind)
`DMLowRank` is read by `DefaultRhoRoute()`.  Every `IrrepCD_Factory` call that does not state a route reaches it:

| caller | what it builds |
|---|---|
| `Orbitals/Internal/Imp/TOrbitals.C:219` (`MakeChargeDensity`) | every iteration's density, in atom, molecule and solid runs alike |
| `Calculation/Imp/SolidState.C:263-264` | a restart's density |
| `ChargeDensity/Imp/Seed.C:143, 151` | the seed |

(`TransitionDensity` already states `RhoRoute::Direct`.)

Why it is the hard one:
1. **It sits below every facade.**  `TOrbitals` holds a basis and a D, not a run.  To reach it, the route must be handed
   through the orbitals' construction chain, and that runs through the `WaveFunction` hub (`hub-interfaces.snapshot:528`), or
   the basis must answer it.  Either is a design discussion first (CLAUDE.md: never route around a hub).
2. **It leaks TODAY, not only in theory.**  Only `SolidSetup::Resolve` installs a policy.  `MoleculeCalculation` and
   `AtomCalculation` never do, so their densities take whatever route the last SOLID run in the process installed.  In a GUI
   session or a test binary, a `cp2kCompat` solid followed by a molecule runs the molecule on the Direct route.  The route is
   `movesState=true` (`Imp/RunPolicy.C:31`), so this is a state-moving leak, not a cosmetic one.  It is the one reader where
   the process-wide value is already giving a run a policy it did not state.

The candidate hand-offs, for the discussion:
- **(E1) The route on the basis**, by the `ToleranceBinding` / `CollocationRoute` precedent: a capability face that the Bloch
  and molecular blocks answer (`RhoRoute`), bound by each facade when it builds its basis.  `IrrepCD_Factory` asks the basis by
  an abstract→abstract cross-cast.  It is not a member on the `BasisSet` hub, but it is a new face that the factory consults,
  so it needs the user's agreement.  It also matches the factory's own comment, *"Policy belongs HERE, in the factory"*:
  the factory still decides, from what the basis was bound with.
- **(E2) The route on the orbitals**: each `TOrbitals` is constructed with a `RhoRoute` and passes it to `IrrepCD_Factory`.
  Clean at the leaf, but the route must thread through every `WaveFunction` factory, so it is a hub signature change.
- **(E3) The density factory stated at each facade**: `SolidState` and `Seed` state the route, and `TOrbitals` keeps the
  default.  This is cheap, but it leaves the per-iteration path, the one that matters, still reading the process.  **Not
  recommended.**

**Recommendation: E1.**  The basis is the object every run kind already builds per run and already binds per-run routes onto,
and binding it closes the molecule/atom leak as a side effect.

## 3. Order, if the user agrees

The groups are independent.  Each commit carries one group, with 0 movers and the banner reading the run's own policy:
1. **A** (setup): also makes the banner print the held policy.
2. **B** (basis), **C** (SCF), **D-Hubbard**: one commit each.
3. **D-XC source**: one commit.  It touches all three sampler families.
4. **E**: only after the design discussion picks E1 or E2.
5. **Retire** `thePolicy()`, `SetRunPolicy` and `ScopedRunPolicy` once nothing reads them.  The test A/Bs then state their
   policy on `RunOptions`.

Until step 5, `SetRunPolicy` stays.  Readers move off it one group at a time, and each group's commit removes its reads.

## 4. Addendum: group E as E2′, the density maker (user ruling 2026-11-23)

**The ruling.**
- E1 is WITHDRAWN.  The user: *"It sounds very odd to me things with names like RunPolicy and a RhoRoute enum should be stored
  and served up by the Orbitalbasis."*
- E2′ is chosen.  The user: *"I like E2'."*
- On where the factored/direct decision is made, the user: *"As kT rises the factoring route gets more and more inefficient
  (useless).  SO there is some high level decision making there that we don't want hard coded inside Orbitals and
  WaveFunction.  1) Density maker could make that decision, or 2) the decision could be pre determined and passed into the
  Density maker at construction time, or 3) Density maker could be an interface with multiple concrete implementations
  representing different policies.  The last is overkill if there just a single policy defined by one enum RhoRoute."*
- Taken: option **(2)**.  The setup decides, and the maker only carries the decision out.

### 4.1 The shape

**`ChargeDensity::DensityMaker`**: a CONCRETE class in `qchem.ChargeDensity.Factory`, with no virtuals.
- It holds the resolved `RhoRoute` and nothing else.
- Its one job is `template <class T> std::unique_ptr<tDM_CD<T>> Make(hmat_t<T> D, const tobs_t<T>* bs, const Irrep& qns) const`,
  which calls `IrrepCD_Factory(D, bs, qns, itsRoute)`.
- It decides nothing and never switches route.

**Who builds it: each SETUP, once.**  The setup holds `std::shared_ptr<const DensityMaker> densityMaker` on the base `Setup`, so
every run kind has one.
- `SolidSetup` resolves the route:
  1. a stated `policy.dmLowRank` wins;
  2. otherwise `cp2kCompat` gives `Direct`;
  3. otherwise the qchem default (`PivotedCholesky`), at every kT: there is no kT rule (§4.3, the user 2026-11-24).
- `MoleculeSetup` and `AtomSetup` take the qchem default; they have no policy block today.
- Each writes its decision and the reason under the run report's `policy` scope.

**The construction chain it travels.**  Each step DEMANDS the maker: no default and no null (CLAUDE.md, constructors pass
concrete objects).
1. the facades: `SolidCalculation.C:717, 1021, 1025`, `MoleculeCalculation.C:220`, `AtomCalculation.C:308, 317`;
2. → `tSCFIterator`'s three constructors (`SCFIterator.C:195, 213, 221`) and `SolidSCFIterator`;
3. → `WaveFunction::Factory` (`qchem.WaveFunction.Factory`);
4. → the polarized / unpolarized / composite wave functions;
5. → `tIrrepWF`'s constructor;
6. → `Orbitals::Factory(bs, ms, maker)` → `TOrbitalsImp(bs, ms, maker)`.

`TOrbitalsImp::MakeChargeDensity` then calls `itsMaker->Make(...)` in place of the free `IrrepCD_Factory`.

**As built (E-a):** the FACADE builds the wave function (`WaveFunction::Factory`, with its setup's maker) and hands it to the
iterator, which owns it.  The iterator never holds the maker: only a strategy-seeded iterator is given one, to make the seed's block
density (`MakeSeedDensity`).  Its constructors lost the ortho / ortho-tolerance / band-count arguments, which now go to
`WaveFunction::Factory`.  The maker is a small VALUE (it holds one enum), passed by const reference and copied into the orbitals.

**Other densities built from a D take the setup's maker too:**
- the restart (`SolidState.C:263-264`);
- the carried states (`MoleculeCalculation.C:464`, `AtomCalculation.C:517`);
- the seed (`Seed.C:143, 151`), whose callers pass it in.

**What goes away:** `DefaultRhoRoute()`, its `theRunPolicy()` read, and the `IrrepCD_Factory` overload without a route.  Every
remaining direct call states its route; `TransitionDensity` already does.

**The hub.**  The protected face of `qchem.WaveFunction` (`tWaveFunction`, `hub-interfaces.snapshot:528`) does not change, and
neither `qchem.WaveFunction.Factory` nor `qchem.Orbitals` is in the snapshot.  So `HubInterfaceAudit` stays green.  The commit
still cites this ruling, because the chain it changes runs through the hub's construction.

**Bare and test constructors.**
- Tests build SCF runs through the facades.  There is one direct `SCFIterator` construction in a test, and no test calls
  `WaveFunction::Factory` or `Orbitals::Factory`.
- That one test, and the 14 test calls to the route-less `IrrepCD_Factory`, STATE their route: a `DensityMaker` or an explicit
  `RhoRoute`.
- There is no default maker helper, which would only re-create the hidden default.

### 4.2 Order
- **E-a.**  Thread the maker, with the route resolved exactly as today: 0 movers.
  - The leak test (`D-RUNPOLICY-LEAK`, knownRed until then) turns green here: a molecule after a `cp2kCompat` solid takes its own
    setup's route.
  - One commit, which cites the ruling.
- **E-b.**  Take the collocation route's setter off the basis (§4.4): 0 movers.
- **E-c.**  The route at kT > 0 (§4.3), after its measurement and the controller's OK on the rule.  It is `movesState`, so its movers
  are reported and judged.

### 4.3 The route at kT > 0: measured first (the user's correction, 2026-11-23)
The first draft said "at kT > 0, Direct".  WITHDRAWN.  The user: *"“at kT > 0, Direct” it is not the cut and dry.  If kT is small
(relative to the gap?) then D still factors nicely with some small pivot cutoff.  The devil is in the details!"*

What decides the cost is the effective rank r/n of D at the pivot tolerance.  That depends on kT RELATIVE TO THE GAP (or to the
density of states near the chemical potential), not on kT alone.
- **The rigorous part:** the cost model, O(npts n r) factored against O(npts n^2) direct, written as its own function.
- **The empirical part:** any threshold on r/n or on kT/gap is a NAMED calibration, with its measurement record at its declaration.
- **Who decides when:** the setup chooses the STARTING route; it cannot know the gap before the SCF.  The maker reports the
  achieved r/n per block every iteration as health and never switches.  The policy layer may revise the route only at the next SCF
  BOUNDARY (a stage, a restart, a geometry step), reading that health, and reports the revision under `policy`.  (The controller's
  reading, 2026-11-23; to be confirmed with the measurement.)
- **The measurement (E-c, before any rule):** r/n and the per-iteration rho cost on both routes, against kT/gap and the pivot
  tolerance, on a GAPPED system at kT well below the gap (Si, or MnO Γ Shub smear) and on a metal (Li k444 Becke smear).  It goes to
  the controller, with the proposed rule, before E-c is coded.

**The result and the ruling (2026-11-24).**  Measured on Si k222 (Becke, gap 0.072 Ha), MnO AFM-II Γ (gap 0.053 Ha) and Li k444 (metal) at kT
0..2e-2 Ha, factored + census / factored / direct, 27 runs (full table: the lane B note `laneB-Ec-measurement.md`; logs
`~/Code/qchem6-runs/{Si,Li,MnO}/laneB_Ec/`).
- Both routes gave the same Etot to 1e-10 Ha in every pair: the factored route is exact (its pivots sit at LAPACK's roundoff floor).
- The thermal tail is SHORT: r counts the bands within ~32 kT of μ.  At kT 0.02 (kT/gap 0.3-0.4) r/n reached 0.16 on MnO and 0.42 on Si;
  no block in any case or iteration reached the existing per-block fallback (2r >= n).
- The pivot tolerance barely matters: r(1e-4) = r(1e-14) up to kT 5e-3, 15 vs 20 on MnO at 0.02 -- a looser cut is an accuracy trade
  worth <= 25 %, so there is NO pivot-tolerance knob.
- The factored route was never slower (direct ÷ factored 1.2-2.2 MnO, ~1.5 Si, ~1.0 Li at n = 5).
- On a uniform raster XC mesh the route never runs (`FactorPays` false): the choice matters only where rho is sampled through the
  orbital basis (Becke).

**THE RULING: NO kT RULE** (the user, 2026-11-24, asked "kT, item (c): do you accept no kT rule": *"yes I am convinced.  Thanks."*).  The
per-block test on the rank each density ACHIEVES is the whole decision: no setup switches the route on kT before the SCF, and no policy
revises it at an SCF boundary (the controller's boundary-revision reading above is WITHDRAWN).  Built as E-c (fe2826687): the test in the
CLAUDE.md shape, `FactoredRhoOpRatio` (rigorous r/n) × `kFactoredRhoMargin` (= 2, the named calibration, unexercised at its threshold,
this record and the ruling at its declaration) in `FactoredForwardPays`, `src/ChargeDensity/Internal/Imp/IrrepCD.C`; and each factored
block's achieved rank Accumulated as health into a `densityFactor` record (densities, factored, fallbacks by reason, Σr, Σn).

### 4.4 The collocation route: kept PRIVATE to the evaluator (option (2), controller 2026-11-23)
`CollocationRoute {contractCube, expRecurrence}` (78e621788) was BOUND ON THE MOLECULAR BASIS through a public setter,
`ToleranceBinding::BindRoute`, and read by its MnD evaluator's kernels (`UseContractCube()` / `UseExpRecurrence()`).  It is not a
run input; only unit tests ever set it.

The first draft proposed passing it PER CALL, beside the `LatticeScreener`.  Coding it showed why it was bound per basis: the
evaluator's box tasks and integrate-back memos are BUILT with the kernel, so a per-call route would have to key those caches on it,
across 7 `LatticeCollocation` virtuals in 4 implementing layers.  Taken instead (controller: *"two exact algorithms for the same
integral tables are an implementation detail inside face 2, like an ERI algorithm.  The objection was the public setter, not the
private member"*):
- `ToleranceBinding::BindRoute` and its `PG_Cart::BasisSet` / `Orbital_IBS` / `NR_Evaluator` forwarders are DELETED.
- `NR_Evaluator::itsRoute` stays PRIVATE, fixed at production; its declaration says why.
- The kernel A/B tests (`M_TransitionCollocation`, `M_CollocationPullback`, `GPWTolerances`) set it through a `forward.H` friend,
  `CollocationRouteTests`, on a fresh basis before any work.

The user's test for what a basis may serve (2026-11-23): *"A basis describes functions ... only indirectly through 1) op(r)/grad(r),
2) integral tables for over operators and 3) Symmetry info exchange.  Anything outside those 3 is subject to scrutiny."*

**Not proposed:** `ToleranceBinding::BindTolerances` (D-ENV-TOL).  Its `screenEps` / `fieldSharp` / ... define the integrals the
basis itself returns (S, T, V's lattice sums), so they are a property of the basis, not a client's route choice.  Flagged for
the user in case the same objection is meant to cover it.
