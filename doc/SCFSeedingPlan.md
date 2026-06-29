# SCF Seeding — Design Spec

Improve how the SCF loop is seeded with an initial charge density. Two pain points, one
architecture. Written for a fresh session; everything below is grounded in the current code.

## 1. Current state (what actually happens today)

The seed enters as the **optional last argument** of the iterator:

```cpp
// src/SCFIterator/SCFIterator.C
tSCFIterator(const tbs_t<T>*, const ElectronConfiguration*, ham_t*, acc_t*, tDM_CD<T>* cd = 0);
```

`Initialize(cd)` (src/SCFIterator/Imp/SCFIterator.C:76) just forwards it:

```cpp
itsWaveFunction->DoSCFIteration(*itsHamiltonian, cd);   // cd == nullptr is legal
itsWaveFunction->FillOrbitals(0.0001);
itsCD = itsWaveFunction->GetChargeDensity();
```

So `cd == nullptr` ⇒ the Hamiltonian is built with **no density** ⇒ only the density-independent
terms (kinetic + external/nuclear) contribute ⇒ a **core-Hamiltonian guess**. There is **no SAD /
atomic-density / superposition machinery anywhere** (verified by grep).

Two callers, two behaviours:

- **Molecular** (`UnitTests/QchemTesterImp.C:72`): passes *no* `cd` ⇒ core guess. Clean, no
  boilerplate, works — but for ionic systems the core guess is far from the answer.
- **Plane-wave** (`UnitTests/PlaneWaveDFTUT.C:1060-1065`): hand-builds a uniform seed, every test:
  ```cpp
  hmat_t<dcmplx> D0 = blazem::zeroH<dcmplx>(n);
  for (size_t i=0;i<n;i++) D0(i,i)=double(Nelec)/double(n);   // D = (N/n) I  -> uniform rho
  auto* seed = new ChargeDensity::IrrepCD<dcmplx>(D0, pw, irr);
  cSCFIterator scf(bs.get(), &ec, ham, acc, seed);
  ```
  (The standalone `RunSCF` test-driver, line ~438, does the same as `R.rho[(0,0,0)] = Nelec/Omega`.)
  PW seeds uniform *on purpose* — the comment notes "Hartree+XC active from iteration 0, as real PW
  codes do." Core-guess (`cd=0`) would leave Hartree/XC off for iteration 0 (free-electron orbitals).

### The two problems
1. **UX / boilerplate.** The PW path makes every test/app hand-assemble a density matrix. The seed
   strategy should be a one-liner (or a sensible default), symmetric with the molecular path.
2. **Crude for charge-transfer.** Uniform (PW) and core-guess (molecular) are both far from an ionic
   ground state. NaF / CsI converge but **grind** — the SCF rebuilds the entire Na→F charge transfer
   from a structureless start.

## 2. Goals / non-goals

**Goals**
- A single `SeedStrategy` abstraction with a sensible **default per basis type**, so no caller
  hand-builds a density matrix.
- A **superposition-of-atomic-densities (SAD)** seed that captures atomic shell structure ⇒ far
  fewer SCF iterations, especially for ionic solids.
- Keep the converged answer **bit-identical**: the seed changes only the *path*, never the minimum.

**Non-goals (this work)**
- Not not new accelerators, not forces.
- Not perfect oxidation-state inference — IonicSAD is an opt-in Phase 3 with a simple heuristic.

## 3. Architecture

The seed is conceptually one function: *given {Structure, basis, electron config}, produce an initial
`tDM_CD<T>` (or null for core-guess)*. Introduce a strategy + factory:

```cpp
enum class SeedStrategy { CoreGuess, Uniform, SAD, IonicSAD };

template<class T>
tDM_CD<T>* MakeSeedDensity(SeedStrategy, const tbs_t<T>* bs,
                           const Structure* st, const ElectronConfiguration* ec);
```

- `CoreGuess` → `nullptr` (today's `cd=0`).
- `Uniform`   → the centralized version of the current PW boilerplate (`D=(N/n)I`).
- `SAD`       → superpose neutral atomic densities (§3.2).
- `IonicSAD`  → superpose **ion** atomic densities (§3.3).

The iterator picks a **default** when none is supplied. Add the strategy to the ctor (or a small
overload) with a basis-resolved default:

```cpp
tSCFIterator(..., SeedStrategy seed = SeedStrategy::Default);
// Default resolves: molecular -> CoreGuess (today), PW -> Uniform (today),
//                   then BOTH -> SAD once §3.2/§3.3 land.
```

This deletes the PW test boilerplate (tests pass a strategy enum, or rely on the default).

### 3.1 The AO/FT axis returns — do NOT invent a third pattern

SAD is "sum atomic densities at the nuclei." The **concept** is shared; the **container bifurcates by
basis type**, exactly like the existing splits ([[project_functionfitter_isp_split]],
[[project_ao_ft_projection_cleanup]], [[project_dft_fitting_boundary_pin]]):

- **Molecular (`double`, AO):** block-diagonal superposition of atomic **density matrices**. Solve
  each unique element once (atomic LDA SCF), drop its converged DM into that atom's diagonal block of
  the molecular DM; off-diagonal blocks zero. Textbook SAD ⇒ an `IrrepCD<double>`.
- **Plane-wave (`dcmplx`, FT):** sum atomic density **form factors** with structure factors,
  `rho(G) = Σ_atoms ρ_atom(|G|) e^{-iG·R_atom}` ⇒ a Fourier density (`IrrepCD<dcmplx>`). This is the
  **same structure-factor assembly the pseudopotential already uses** (`MakeNuclear` / the local-PP
  form factor in `Integrals_Pseudo`/`PW_Pseudo`) — feed the *atomic valence density* radial function
  through that path instead of the PP.

So realize `MakeSeedDensity` as a **polymorphic capability split across the two basis libraries**
(molecule_BS provides the AO/block-DM face, lattice_BS provides the FT/form-factor face), with the
ChargeDensity/iterator layer agnostic — mirroring the cross-cast capability pattern already in the
tree. **Do not add a variant or a third mechanism.**

### 3.2 SAD — atomic densities

A density matrix is **basis-specific** (`D_ab` only means anything relative to its `{χ_a}`), so the
molecular face needs an atomic density *in the molecular basis block*, not in some foreign atomic
basis. Two routes; **Route 1 is the default, Route 2 is the escape hatch.**

**Route 1 — solve-in-block (no conversion, the standard SAD trick).** The molecular basis restricted to
one centre *is* a perfectly good atomic basis. So for each unique element, run a **single-atom SCF using
the same basis-set spec the molecule uses** (a one-atom `Structure` + the molecular Gaussian basis →
reuses the molecular SCF machinery wholesale). The resulting `D` is **natively in the molecular block** —
scatter it onto every atom of that element; inter-atom blocks = 0. No DM conversion, because there was
never a basis mismatch. Solve once per unique element, cache.
- **Spherical averaging:** an open-shell single-atom solve gives an *aspherical* DM (the
  maximally-stretched-determinant orientation — cf. the boron p-pancake). Prefer **fractional /
  spherically-averaged occupations** (O's four 2p as 4/3 per component) so the seed carries no arbitrary
  orientation. *Pragmatic caveat (user):* the seed only needs to be "in the infield," so an aspherical
  DM is often **good enough** and the averaging can be skipped for easy cases. The one place it bites is
  **near-degenerate / symmetric / magnetic** systems (TM oxides, Jahn-Teller) — there an oriented seed
  can break symmetry or tip the SCF into a wrong basin (same bifurcation family as the M_DFT_Water
  convergence cliff). So: spherical for the robust default; aspherical acceptable as a cheap fast-path.

**Route 2 — project a foreign atomic density** (a tabulated/numerical radial density, or a richer atomic
basis than the molecule's). Least-squares **project** `ρ_atom` onto the molecular basis block — and the
machinery already exists: `FunctionFitter_Density` is exactly "fit a density onto a basis under the
Coulomb/overlap metric." Lossy (projection error) and unnecessary whenever Route 1 applies, but the way
to consume an external density.

**The SAD data file (eventually).** Route 2's natural backing store is a **basis-independent radial
atomic-density database** (per element × functional), shaped like `saito.json` / the GTH JSON — store
`ρ_atom(r)` on a **radial grid**, not a DM and not a Gaussian/spline expansion. *Why a grid, not an
expansion:* both consumers **quadrature it anyway** — the molecular Route-2 projection is
`Σ_g w_g ρ_atom(r_g) χ_a(r_g) χ_b(r_g)` on the mesh, the PW form factor is a radial transform
`4π ∫ ρ_atom(r) sinc(Gr) r² dr` — the *same* mesh-quadrature pattern the molecular PP (Path A) and XC
paths already run. So a grid feeds straight into existing machinery with **no fit step, no exponent
choices, no fit conditioning**, and it's the direct output of the atomic solve; the analytic-integral
advantage of an expansion doesn't pay off for a ballpark seed. **Store the grid-spec WITH the values**
(grid type + `r_min`/`r_max`/`N` for a log-radial grid, or explicit abscissae) — not bare numbers.
Consumers are on *different* grids (the Becke mesh, the FFT grid), so they **interpolate** the stored
radial density onto their own quadrature/FFT points; lossless-for-seed-purposes since a PP-smoothed
valence density is smooth. Keep the stored grid a standard log-radial, fine enough that interpolation
is a non-issue. One file then serves **both faces**: the **PW** face radial-FTs it → `ρ_atom(|G|)` form
factor; the **molecular** face projects it (Route 2) → block. That unifies AO + FT SAD onto a single
source and removes the per-run atomic solve. Order of adoption: Route 1 (solve-in-block) is the quick,
file-free start; the radial-grid file (`{grid-spec, ρ_atom(r)[, spin]}` per element × functional, JSON)
is the later optimization that also feeds the PW form factor.

**PW face:** the atomic **valence** radial density FT'd to `ρ_atom(|G|)`, summed with structure factors
(§3.1) — the *same* `ρ_atom(r)` the data file would hold (no basis-block issue: PW has no atom-centred
blocks to reconcile).

**Spin:** thread `Pol` — open-shell wants a spin-polarized atomic seed (atomic spin densities).

### 3.3 IonicSAD (Phase 3, opt-in)

Neutral SAD still leaves ~1 e⁻ of Na→F transfer for the SCF. IonicSAD seeds Na⁺ + F⁻ atomic
densities ⇒ near-converged start for NaF/CsI. Needs **oxidation states**: a simple
electronegativity/formal-charge heuristic, or user-supplied (could ride on the Structure alongside
the Zion/PP work — see [[project_hgh_params_table]]). Opt-in and user-overridable; the heuristic is
the only genuinely fuzzy piece in this whole plan, so isolate it.

## 4. Phased plan (each phase regression-safe)

- **Phase 0 — refactor, zero behaviour change.** Add `SeedStrategy` + `MakeSeedDensity`; move the PW
  uniform boilerplate into `MakeSeedDensity(Uniform,…)`; iterator default resolves PW→Uniform,
  molecular→CoreGuess. *Gate: every SCF test identical iteration count + energy.*
- **Phase 1 — SAD molecular** (block-diagonal atomic DM). *Gate: molecular DFT energies unchanged;
  iteration count ≤ today.*
- **Phase 2 — SAD plane-wave** (G-space form factors, reuse structure-factor assembly). Flip the PW
  default to SAD. *Gate: Si / NaF / CsI energies unchanged; iteration count strictly down.*
- **Phase 3 — IonicSAD** (opt-in, oxidation heuristic). *Gate: NaF/CsI iteration count down further;
  user-override path tested.*

## 5. Validation
- **Hard invariant:** converged total energy bit-identical across seed strategies (pin existing
  energy regressions; they must not move). The seed only changes the path.
- **Soft guards:** iteration-count assertions on NaF/CsI/Si with *generous* bounds (iteration count is
  sensitive to accelerator + tol; assert "≤ baseline", not an exact number).
- Keep `CoreGuess` reachable — some tests/edge cases may want the structureless start.

## 6. Risks / nuances
- **Total charge** must integrate to exactly `Nelec` (neutral SAD does by construction; ionic must
  balance the formal charges).
- **Representability:** PW seed populates only G-vectors present in the basis; molecular seed only the
  atom blocks present. Truncation must not lose charge (renormalize after projection).
- **Atomic-SCF cost** is negligible (≤ #species solves, cached) but must use the *same functional* as
  the target or the seed fights the SCF.
- **"Hartree active from iter 0":** SAD/Uniform give nonzero ρ at iter 0 (good); CoreGuess does not.
  SAD strictly dominates Uniform (it has shell structure), so SAD becomes the default once it lands.
- **Don't re-fork the AO/FT axis** — follow §3.1.

## 7. Open questions for the implementer
1. Atomic densities **live-solved** (always functional-consistent) vs **tabulated** (faster startup,
   like `saito.json`)? Recommend live-solved + cached; tabulate later only if startup cost bites.
2. Home for `MakeSeedDensity`: a tiny qcChargeDensity module for the dispatch, with the two faces in
   molecule_BS / lattice_BS (recommended, consistent with the tree) — confirm against the current
   library DAG before wiring.
3. IonicSAD oxidation source: heuristic vs Structure-carried formal charges vs explicit user arg.

## 8. Key file pointers
- Seed entry: `src/SCFIterator/Imp/SCFIterator.C:76` (`Initialize`), ctor `:61`.
- Molecular caller (clean default): `UnitTests/QchemTesterImp.C:72`.
- PW boilerplate to delete: `UnitTests/PlaneWaveDFTUT.C:1060-1065` (and `RunSCF` ~438).
- PW structure-factor assembly to reuse for SAD form factors: `PW_Pseudo` / `Integrals_Pseudo`
  (`src/Hamiltonian/Internal/Imp/PWTerms.C`, `src/BasisSet/Lattice_3D/Imp/PlaneWave_IBS.C`
  `MakeNuclear`/`PseudoG0Energy`).
- Density container: `src/ChargeDensity/Internal/Imp/IrrepCD.C` (`IrrepCD<T>`, the `tDM_CD<T>`).

---

## 9. RESOLVED DESIGN (decided in a design session; supersedes the Route-1-first assumption above)

> Phase 0 is **DONE & committed** (`732eace2`): `SeedStrategy{Default,CoreGuess,Uniform,SAD,IonicSAD}` +
> `MakeSeedDensity<T>(s,bs,st,ec)` in new module `qchem.ChargeDensity.Seed`; `Uniform`=`D=(N/n)I`;
> `Default` resolves `double`->`CoreGuess`, `dcmplx`->`Uniform`; iterator ctor's last arg is now
> `SeedStrategy seed=Default`; PW test boilerplate deleted. 134/134 green, iters/energies bit-identical.

### 9.1 Phase 1 = data-file SAD, **DFT-only** (route reversal)

We jump **straight to the §3.2 "SAD data file" route, NOT Route 1 (live atom solves).** Rationale:
relying on invisible nested SCF convergence inside every `SCFIterator` construction is fragile/optimistic;
a curated, inspectable, committed file is robust. The atom solves move **offline** (one-time generation),
not runtime.

**Layering win — this dissolves the orchestration problem.** A *live* atom solve needs the Hamiltonian
factory + accelerator + `SCFIterator` (top of the DAG) and the calc config (`Pol`/functional/mesh/basis-spec)
that lives only at the `Hamiltonian::Factory` caller (`QchemTester`/`M_DFT`), NOT in the iterator. The
data-file route needs none of that: the runtime path is **read-file -> interpolate -> fit -> assemble**, all
inside `qcChargeDensity`'s existing deps (`qcFitting`, `qcMesh`, `qcStructure`, `qcBasisSet`). So Phase 0's
`SeedStrategy` ctor and the `MakeSeedDensity(strategy,bs,st,ec)` signature **survive unchanged** — no
Hamiltonian, no nested SCF, no config threading.

**DFT-only.** Gate is "molecular DFT energies unchanged." A fitted-density seed CANNOT seed HF: `Ham_HF_*`
builds **both** J and K from the density **matrix** via `AccumulateDirect`/`AccumulateExchange` (the DM-only
methods a fit NA-asserts), so HF needs a real `D_ab` on iteration 0. DFT (`Ham_DFT_*`) consumes the density
only through `rho(r)` (XC, `FittedVxc`) and `<rho|c>` (fitted Coulomb, `FittedVee::DoFit` -> cross-cast to
`ProjectedDensity_AO`). So **HF keeps `CoreGuess`**; HF SAD (needs a DM seed = Route-1 territory, or a
fit->DM projection) is a later phase.

### 9.2 The seed object (the real new type-work)

The seed must wear **three faces**: `tDM_CD<double>` + `ProjectedDensity_AO` + `ScalarFunction`, and
**NA-assert the DM-only HF methods** (`AccumulateDirect`/`AccumulateExchange`), exactly like the dcmplx PW
path. Why: `DoSCFIteration(tHamiltonian<T>&, const tDM_CD<T>*)` requires a `tDM_CD`; `FittedVee::DoFit`
cross-casts it to `ProjectedDensity_AO` and asks `GetRepulsion3C=<rho|c>` + `FitGetConstraint=N`; `FittedVxc`
samples `operator()(r)`. A fitted density supplies all of these from its fit coefficients
(`<rho|c> = sum_d e_d <d|c>`, `rho(r)` directly) — **no DM**.

Today's `FittedCD` is only a `ScalarFunction` and fits **from** a `DM_CD` (`DoFit(const DM_CD&)`), so it is
NOT usable as-is. **New type needed:** a fitted-density that implements the three faces above (working name
e.g. `SeededCD`/`FitSeed_CD`), built by fitting a **scalar** radial function.

### 9.3 Assembly = per-atom split + CompositeCD (chosen over a single whole-structure fit)

1. Split the molecular CD-fit basis into **per-atom `FIT_CD_ABS` groups** (basis primitive; the PG basis
   owns it, structure-atom order — functions are already grouped by centre, see `IrrepBasisSet.C:65`).
2. For each atom, fit its **recentred** radial `rho_elem(r)` onto that atom's group via
   **`FunctionFitter_Scalar`** (the `ScalarFFClient`/`GetScalarFunction()` path — it takes a raw
   `ScalarFunction`). NOTE the Coulomb-metric `FunctionFitter_Density` is **unavailable** here: its RHS
   `GetRepulsion3C = sum_ab D_ab<ab|c>` needs a DM. So overlap-metric fit + **renormalize to the atom's
   electron count** (good enough for a seed).
3. Wrap each fit as the new fitted-`tDM_CD`; `Insert` into a `tComposite_CD<double>`. Renormalize the total
   to `Nelec` (§6 representability).

(`Seed()`/`MapSubBasis()` from earlier brainstorming are **dropped** — Route-1 artifacts.)

### 9.4 The data file + generator

- **File:** one neutral-atom radial density per element on a **log-radial grid**, JSON shaped like
  `saito.json`/the GTH JSON: store the **grid-spec WITH the values** (`{kind, r_min, r_max, N}` or explicit
  abscissae) `+` optional spin. Reader mirrors `GetGTH`. Consumers interpolate onto their own grid. The same
  file later feeds the **PW form factor** (Phase 2).
- **Functional-agnostic for now:** generate with **LDA**; other functionals read it and eat +/-1 iteration
  (the seed only changes the path, not the converged minimum — the hard invariant). Later: a few files for
  popular functionals; power users supply custom files.
- **Generator: extend `scfrun.C` itself** (one-atom DFT driver, dumps `rho(r)`). scfrun's `--model` enum is
  `HF|DHF|E1|DE1` with **no DFT** (DFT uses the other `Factory(Pol,st,alpha|ex,mesh,bs)` overload), so add a
  DFT path + `--xc <LDA|...>` / `--alpha <float>` + `--out <json>` (or `--dump-density`), reusing the existing
  flag vocabulary verbatim (`--Z --q --pol --basis --acc --accel --maxiter --minro/--minde/--virial/--minfd
  --relax`, and the `--flag value` / `--flag=value` parser).

### 9.5 Implementation order (each step regression-safe) -- ALL DONE

a. **DONE** (`7aca367a`) -- `scfrun --model LDA --out` generates `Data/atomic_densities.json` (H..Ne; each
   `rho` integrates to `Z`).
b. **DONE** (`a24fb90d`) -- `qchem.ChargeDensity.AtomicDensity`: `RadialDensity` reader + interpolation +
   `RecentredAtomicDensity` (3-D `rho(|r-R|)`); `AtomicDensityUT`.
c. **DONE** (`2084fcd3`) -- `CompositeFittedCD`: superposition of per-atom `rho(r)`; `op(r)` for Vxc;
   `GetRepulsion3C` derived on demand (overlap-fit onto the term's own basis, then Coulomb-project), so
   `op(r)` is the only state. Nominally `tDM_CD` with `assert(false)` HF/DM stubs (the deferred
   `tChargeDensity`/`tDM_CD` ISP, "option 1", is the clean fix).
d. **OBVIATED** -- the on-demand fit in (c) removes the need for a per-atom `FIT_CD_ABS` split.
e/f. **DONE** (`65b06892`) -- `MakeSeedDensity(SAD)` builds the `CompositeFittedCD` from the structure;
   `Structure` threaded via a new `SCFIterator` ctor arg + `QchemTester::SetSeedStrategy`; M_DFT (N2,
   water cart+sph) opt into SAD. **Gate met:** converged energies match the CoreGuess baseline (water
   bit-identical; N2/Sph to ~1e-10, the 20-iter-cap floor). SAD is DFT-only; molecular `Default` stays
   `CoreGuess`, DFT opts in. 138/138 green.

> **The ISP (option 1) is DONE** (`007b9c93`). `tDM_CD<T>` was split into a DFT-only `tChargeDensity<T>`
> base (rho(r)/GetTotalCharge/Version/ReScale -- the Fock-build face the framework now takes) + a `tDM_CD<T>`
> derived face adding the matrix-only ops (DM_Contract, MixIn/GetChangeFrom, HF AccumulateDirect/Exchange).
> Energy/loop keep `tDM_CD`; the 3 HF sites (Vee/Vxc/VxcPol) `dynamic_cast` to `tDM_CD`. `CompositeFittedCD`
> is now a clean `tChargeDensity` with **zero `assert(false)` stubs**. Release 138/138 + Debug 117/117 green.

### 9.6 Deferred

Route-1 live solve; per-element caching (moot — file replaces solves); per-functional file entries; HF SAD;
spherical averaging (use the file's reference density as-is); the `tDM_CD::Seed()` virtual.

### 9.7 Phase 2 (plane-wave SAD) progress

- **(a) DONE** (`3682c14e`) — pseudo-valence density generator. **A cubic plane-wave box is the wrong tool
  for a spherical atom** (breaks symmetry, fights the open-shell 3p degeneracy, O(n³), didn't converge) —
  so reuse the robust **spherical all-electron** atom solver and **extract the valence shells**:
  `scfrun --model LDA --valence <n>` sorts the converged atom's occupied orbitals by energy, accumulates
  whole orbitals until `<n>` electrons, and dumps `rho_val(r)=Σ occ·|φ(r)|²` (kind=valence) to
  `Data/atomic_valence_densities.json`. Si q4: 17 iters, charge=4.0000, 4πr²ρ peaks at r=2.08 bohr. The
  core is all-electron-peaked vs a true pseudo density — accepted (a seed only needs the bonding region in
  the ballpark). The **true pseudo-valence radial PP solver is deferred to the Molecular-PPs project**
  (a real-space Atom-PP solver belongs there).
- **(b) TODO** — `FourierSeedCD`: a `tChargeDensity<dcmplx>` providing `GetFourierDensity()=ρ(G)` (the
  structure-factor sum `Σ_atoms ρ_atom(|G|) e^{-iG·R}`, FT analog of `CompositeFittedCD`), consumed by
  `PW_Hartree`/`PW_XC` via the existing `FourierDensity` cross-cast.
- **(c) TODO** — flip the PW default to SAD; gate Si energy unchanged + **iters strictly down** (PW tests
  converge-and-stop, so the win is visible here). NaF/CsI get a partial (neutral-SAD) win now; the
  charge-transfer payoff is Phase 3 IonicSAD.
