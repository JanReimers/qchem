# Molecular Pseudopotentials — Interface Divergence & Harmonization Findings

Companion to `doc/MolecularPseudopotentialPlan.md` (user-owned; not edited). This is **goal 2** of the
Molecule_PP project: *having got it working (goal 1), record the interface divergences it exposed between
the molecular (`src/BasisSet/Molecule`) and plane-wave (`src/BasisSet/Lattice_3D`) sides, and the concrete
path to harmonizing them* — in service of the long-term goal of hoisting structure-neutral code up into
`src/BasisSet`.

Status at time of writing: Si₂ molecular PP converges through the `qchem::Calculation` facade
(`A_PP.Si2_PP_U.LargeSeparation`), validated by (a) `Enn == Zion²/R` exactly and (b) `E(Si₂) ≈ 2·E(Si)`
additivity at large R. Full UTMain green.

## 1. The headline divergence: *where* PP assembly lives

| | assembly site | reached via | term(s) |
|---|---|---|---|
| **PW** (`Lattice_3D`) | the **basis** — `Integrals_Pseudo<dcmplx>` on `PlaneWave_IBS` | `dynamic_cast` basis→capability, `MakeLocalPotential`/`MakeSeparablePotential` | one: `PW_Pseudo` |
| **molecular** (`Molecule`) | the **term** — `MakeMolecularMesh`+`WeightedOverlap` inline | nothing (term owns the mesh) | two: `PP_Local`, `PP_NonLocal` |

`Integrals_Pseudo<T>` is realized **only** by the PW basis (`<dcmplx>`); there is **no `<double>` impl**. The
molecular path never touches it. So the plan's §1/§3 vision ("molecular realizes `Integrals_Pseudo<double>`")
was quietly superseded by assembly-in-the-term. Net: two assembly sites, nothing hoistable.

## 2. The key realization: the mesh already lives behind the `Fit_ABS` network

The molecular "assembly-in-the-term" is **not** the molecular convention — it is `PP_Local` being the odd
one out. Molecular XC assembly is *already* structure-neutral and mesh-hidden:

- `Fit_IBS` owns `qcMesh::Mesh itsMesh` + `SetMesh(Structure, MeshParams)` (`src/BasisSet/Fit_IBS.C`).
- `FIT_SF_ABS::Overlap(const Sf& f)` projects an **arbitrary scalar field** onto the fit basis over that
  mesh — the fit basis does the quadrature.
- `FunctionFitter_Scalar<T>` / `FunctionFitter_Density<T>` are structure-neutral templates with a molecular
  (Becke-mesh) impl **and** a `FourierFunctionFitter` (G-space) impl. `FittedVxc` holds a fitter, hands it a
  `ScalarFunction`, gets a matrix back, and **never sees a mesh**.

So the tension "the molecular *orbital* basis would have to become mesh-aware to realize
`Integrals_Pseudo<double>`" is a **false framing**: mesh-awareness belongs in the `Fit_ABS`/fitter network,
which both PW and molecular already populate.

## 3. Harmonization target

Route PP assembly through the `Fit_ABS` network the way XC already does. Two neutral primitives are needed
beside the existing fitter (each with a molecular Becke-mesh impl and a PW G-space impl):

- **(a) scalar-field → operator matrix** `⟨χᵢ|V(r)|χⱼ⟩` (local PP). ✅ **DONE.** The right abstraction turned
  out to be the **mesh on the geometry**, not a field-operator capability on the basis. The mesh is a virtual
  on `Structure`: **`Structure::CreateIntegrationMesh(mp)` (pure virtual)** — each geometry owns its most
  efficient mesh: `Atom` → single-centre radial×angular, `Molecule` → multi-centre Becke, `UnitCell` →
  uniform/unit-cell-Becke (throws until implemented; PW never asks — it owns its G-grid). No central `if/case`
  dispatch (SOLID). `PP_Local::CalculateMatrix` calls that virtual directly and quadratures with the generic
  `qcMesh::WeightedOverlap(mesh, basis, V_loc)` — so the term is geometry-neutral with **no `dynamic_cast`**
  and no basis-side indirection. V_loc is **static + smooth**, so this is raw quadrature (no fit), correctly
  distinct from the density/potential *fitter*. Bit-identical (same Becke mesh, same `WeightedOverlap`).
  *(An earlier iteration routed this through a `Mesh_Integrated_IBS`/`MeshIntegratorSource` basis capability;
  it was removed as over-abstraction — PW has its own PP term, so nothing needed the basis to manufacture the
  integrator once the mesh lives on the geometry. `qcStructure` already links `qcMesh`, and `qchem.Mesh`
  imports only `qchem.Types`, so coupling `Structure` to it costs ~nothing.)*
- **(b) scalar-field → projection vector** `⟨χᵢ|β_p Yₗₘ⟩` (nonlocal KB). ✅ **DONE.** Same pattern as (a): the
  projection primitive `qcMesh::Overlap(mesh, basis, β_p·Y_lm) → b` already existed and is generic, so the only
  coupling was the mesh source. `PP_NonLocal::CalculateMatrix` now calls `theStructure->CreateIntegrationMesh(mp)`
  directly (was inline `MakeMolecularMesh`) and builds each projector's `b` via `qcMesh::Overlap`, accumulating
  the rank-1 `D|b⟩⟨b|`. Geometry-neutral, no `dynamic_cast`. Bit-identical (`Si_PP_U` = the KB test: −3.3369 vs
  local-only −11.85). Both PP terms now assemble on the geometry's mesh with generic `qcMesh` quadrature — no
  inline `MakeMolecularMesh`, no basis-side indirection.

With both PP integral types routed on `Structure::CreateIntegrationMesh` + generic `qcMesh` quadrature, the
molecular local + KB pseudopotential assembly is fully geometry-neutral. `Integrals_Pseudo<T>` stays the
PW-only (G-space) capability; a future PW-vs-molecular PP unification is a separate, optional step. The broader
"fold Fitting behind `Fit_ABS`" refactor (for the density/potential *fit* path, cf. §3.1) remains the recommended
framework increment before Path B (semilocal ECPs).

## 4. Multi-species molecular PP (DONE) + the convergence follow-up it exposed

Multi-species molecular PP now works end-to-end: `Ham_PP_U(st, {{elem,q}…})` + `Hamiltonian::Factory(st, species, …)`
build one `MultiSpecies_Local` + one `MultiSpecies_Separable` per-Z router (mirroring the PW `BuildFromGTH`), and
the facade (`MakeValenceStructure` per-atom `Z−Zion`, `PPSpecies` distinct-species list) routes any molecule —
single species is the 1-element case (Si tests still bit-identical through this path). The terms already indexed
on the atoms' `itsZ`, so each atom gets its own local + KB potential and its own `Zion` for the ion-ion.
Validated by `A_PP.OSi_PP_U.MultiSpeciesRouting`: **`Enn = Zion_O·Zion_Si/R = 6·4/R` exactly** (a mis-route
would give 16/R or 36/R) — a convergence-independent proof of per-species routing.

**Follow-up this exposed (real, separate work):** absolute hetero-molecule *energies* are not yet trustworthy.
The molecular-PP SCF via the `Calculation` facade does **not fully converge** for harder closed-shell atoms
(e.g. O 2p⁴): the O–Si additivity `E(OSi) ≈ E(O)+E(Si)` failed by ~8 Ha, and even single Si only oscillates to
~1e-6 (Si2's pins are a deterministic snapshot). Two root causes, both orthogonal to the routing feature:
1. **Valence basis quality** — the ad-hoc `sipp.bsd` exponents are not GTH-optimized; a tight function under a
   smooth GTH PP can over-bind. Production needs proper GTH-optimized molecular bases (SZV/DZVP-GTH).
2. **SCF convergence** — the molecular facade drives PP with a core-guess seed + DIIS-from-start (`EMax=100`);
   harder atoms need a better seed (SAD-for-PP) and/or accelerator tuning, and the closed-shell `Molecule_EC`
   treatment of an open-shell atom (O) is itself a poor state. On the battery critical path (transition metals
   are harder still), so this convergence + basis hardening is the next substantive Molecule_PP increment.

### 3.1 HARD CONSTRAINT: the fit basis is GENERATED BY the orbital basis (the *process* matters)

The fit basis must be produced by the orbital basis through its own factory methods:

```cpp
virtual FIT_CD_ABS* CreateCDFitBasisSet (const Structure*, const qcMesh::MeshParams&) const;  // density-fit basis
virtual FIT_SF_ABS* CreateVxcFitBasisSet(const Structure*, const qcMesh::MeshParams&) const;  // potential-fit basis
```

The **molecular** path honours this: `FittedVee`/`FittedVxc`/`Ham_PP_U` call `bs->CreateXxxFitBasisSet(...)`,
get back a *distinct* `FIT_*_ABS`, and build the fitter on it — the orbital basis is the FACTORY of its fit
basis (a real, separately-tuned auxiliary basis).

The **PW** path VIOLATES this and must NOT be copied:
- `tBasisSet<dcmplx>::CreateCDFitBasisSet` / `CreateVxcFitBasisSet` are `assert(false)` stubs
  (`src/BasisSet/Imp/BasisSet.C`).
- `PW_Hartree`/`PW_XC` (`src/Hamiltonian/Internal/Imp/PWTerms.C`) instead **default-construct** a
  `FourierFunctionFitter` and hand it the *orbital* basis directly — hardcoding "fit basis ≡ orbital basis"
  and bypassing the generating process.

**The clean target:** PW *implements* `CreateCDFitBasisSet`/`CreateVxcFitBasisSet` (the DFT/`Band_FT_IBS`
side, with the composite delegating as the `<double>` path already does). PW is free to *return* itself, or a
more efficient tuned-{G} fit basis — but the caller always goes THROUGH the factory method. `PW_Hartree`/
`PW_XC` then obtain their fitter exactly as `FittedVee`/`FittedVxc` do: `MakeXxxFitter(bs->CreateXxxFitBasisSet(...))`.
The `assert(false)` stubs die; the process is uniform even when PW's answer is trivial.

This constraint governs the whole harmonization: whenever a fit/auxiliary structure is needed, obtain it from
the orbital basis via its factory method — never assume the orbital basis IS the fit/aux basis. (For the raw
field→operator primitive (a) there is no fit basis at all — direct `⟨χ|V|χ⟩` quadrature; the concern applies
to any density/potential *fit* the assembler performs.)

## 4. Smaller divergences found

- **`PseudoG0Energy` — ELIMINATED this session.** The PP-specific G=0 alignment `(N/Ω)·Σₐ α` moved off the
  `Integrals_Pseudo` basis interface into the `PW_Pseudo` term (which reads Ω from `UnitCell::GetCellVolume()`
  and α from the model it already owns). `Integrals_Pseudo<T>` is now two clean universal matrix methods;
  PW energies bit-identical.
- **`PW_IonIon` vs molecular `Vnn` — a T-template candidate.** Both are energy-only ion-ion terms delegating
  to `NuclearRepulsion(st, zionOf)`; they differ only in scalar type (`chmat_t` vs `rsmat_t`) and
  finite/periodic (already handled inside `NuclearRepulsion`). This session unified the **Zion callback**
  onto `Vnn` (one term serves all-electron via an identity default and PP via a Z→Zion map, mirroring
  `PW_IonIon`). A future `IonIon<T>` term could collapse the two entirely.
- **Electron-count coupling (wart).** `Ham_PP_U`/`FittedVee` read the valence count from
  `st->GetNumElectrons()`, forcing the "atom charge = Z−Zion" encoding (real atoms carry charge 0). The
  facade absorbs it via `MakeValenceStructure`, but a cleaner design passes the count/EC explicitly rather
  than deriving it from structure net-charges.
- **Basis provenance.** The molecular Factory is Gaussian-file-only, and all-electron `.bsd` files (e.g.
  `dzvp`) carry core shells that are inconsistent with a PP. A valence-only `sipp.bsd` was authored for
  validation; a proper GTH-optimized molecular basis family (SZV/DZVP-GTH) is the production follow-up.
- **Single-species only.** The facade PP path requires all atoms to be the same element. The
  `MultiSpecies_{Local,Separable}Potential` routers + per-Z `zionOf` (already built and used by PW) are the
  next increment for hetero molecules.

## 5. What landed this session (goal 1 + the cheap goal-2 wins)

1. `PseudoG0Energy` eliminated (A0 closeout; PW bit-identical).
2. `Vnn` gained an optional Z→Zion callback (all-electron default = identity) — one ion-ion term for both.
3. `Ham_PP_U` explicit ctor now takes the **combined** `LocalPotential` (real-space view for `PP_Local` +
   `Zion` for `Vnn`) and adds `Vnn(Zion)` — harmless (Enn=0) for a lone atom, correct for a molecule.
4. `qchem::Calculation` facade: `{.pseudopotential=true}` (single-species), valence-ion structure, PP front
   door, DIIS-from-start.
5. `sipp.bsd` valence Si Gaussian basis + Factory registration; `A_PP.Si2_PP_U.LargeSeparation` anchor.
