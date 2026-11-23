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
