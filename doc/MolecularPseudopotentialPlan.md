# Molecular Pseudopotentials — Design Spec

Bring pseudopotential support to the **molecular (Gaussian, real-space) basis**, reusing the
pseudo-wall built for plane waves. Beyond the feature, this is the forcing function that proves the
PP interfaces are molecule/solid agnostic — and it drove several interface refinements below.
Written for a fresh session; grounded in the current code.

## 0. Two paths (do A, then B)

- **Path A — GTH/HGH in a Gaussian basis** (the CP2K-flavoured choice). Reuses the *existing*
  KB-separable + local models verbatim. Adds **zero new integral types** to the molecular basis. Its
  job is to *prove* `Integrals_Pseudo<double>` is real and bury G-space as a PW-only detail.
- **Path B — standard semilocal ECPs** (def2-ECP / LANL2 / Stuttgart via libcint). These use the
  `V_l(r) P_l` form (radial × angular-momentum projector) — **not** KB-separable, a genuinely new
  integral type (Type-2 angular). Path B *is* "Topic 2: a new PP type that forces model neutrality"
  done right, and it unlocks the published heavy-atom ECP libraries (Co/Ni/Mn for batteries).

Physics note that sets the order: plane-wave codes universally use **KB-fully-separable** (GTH are
*natively* separable). The older **semilocal** form (common radial `V_l(r)`, separate angular `P_l`)
is *not* obsolete — molecular Gaussian codes never left it; the published ECP libraries *are*
semilocal. So Path A is the cheap wall-validation; Path B makes us a real molecular-ECP consumer.

## 1. `Integrals_Pseudo<T>` — shrink to two universal methods

Current interface (`src/Pseudopotential/Integrals_Pseudo.C`):
```cpp
virtual hmat_t<T> MakeLocalPotential    (const Structure*, const LocalPotential&)    const=0;  // KEEP
virtual hmat_t<T> MakeSeparablePotential (const Structure*, const SeparablePotential&) const=0; // KEEP
virtual double    PseudoG0Energy         (const Structure*, const LocalPotential&, double N) const=0; // ELIMINATE
```
The two matrix methods are already **"model in → `hmat_t<T>` out"** — no `FourierMap`, no G-vectors,
no Ω in the signatures. So the FourierMap-is-internal proof already holds at the contract level; the
molecular impl is a drop-in returning `hmat_t<double>` (= `rsmat_t`; molecular PP blocks are
real-symmetric).

### Eliminate `PseudoG0Energy` (don't relocate it)
The PW impl (`PlaneWave_IBS.C:250`) is `(N/Ω)·Σ_a loc.FormFactorG0(Zₐ)`. **α comes from the model**
(`LocalPotential::FormFactorG0`); **Ω is just a number.** So it isn't assembly at all — the term can
compute it from `{Structure, model, N}`, which it already holds. Therefore:
- Delete `PseudoG0Energy` from `Integrals_Pseudo`. The capability becomes the two pure matrix methods —
  fully universal, zero leak, zero orphan.
- `PW_Pseudo` (term) computes the alignment itself. Molecular: no alignment term exists; `FormFactorG0`
  is simply never called (nothing returns 0.0 because nothing asks).
- **Verify:** Ω reachable term-side via `Structure → UnitCell → det`. If yes, term needs no basis. If
  Ω is locked in the basis, expose a *neutral* `CellVolume()` getter (plain geometry, not a PP method) —
  the PP-specific formula still moves term-side.

## 2. The model-level reciprocal-space bias (the real find)

Both `LocalPotential::FormFactor(Z,G2)` and `SeparablePotential::Projector(Z,p,q)` expose their radial
functions **only through the reciprocal view**. That's the *same* reciprocal-space bias `FourierMap`
had on the density side — the **AO/FT axis, one level down in the PP models** (cf.
[[project_ao_ft_projection_cleanup]]). A radial function has two spectral views, `V(r)` and `ṽ(G²)`;
the models expose only `ṽ`. Molecular assembly needs `V(r)`.

### Fix: each radial function exposes BOTH views; the basis consumes the one matching its space.

**`LocalPotential` — ISP into three facets** (it has no neutral structure, so it splits cleanly):

| facet | methods | consumer |
|---|---|---|
| `LocalPotential_Q` (reciprocal) | `FormFactor(Z,G²)`, `FormFactorG0(Z)` | PW assembly + PW G0-alignment |
| `LocalPotential_R` (real) | `Vloc(Z,r)` | molecular mesh assembly |
| charge | `Zion(Z)` | Ewald (`NuclearRepulsion`) — orthogonal to both views |

Select the view by the **same `T`** that selects the matrix type (the existing `conditional_t`
pattern, `FourierDensityBase<T>`):
```cpp
template<class T> struct LocalView;
template<> struct LocalView<dcmplx>{ using type = LocalPotential_Q; }; // PW        -> reciprocal
template<> struct LocalView<double>{ using type = LocalPotential_R; }; // molecular -> real
// Integrals_Pseudo<T>::MakeLocalPotential(const Structure*, const typename LocalView<T>::type&)
```
The compiler now enforces the wall: a reciprocal-only model **is-a** `LocalPotential_Q` and *cannot be
passed* to molecular assembly — capability expressed as a type, no NA-assert. Combined model is the
(harmless) diamond: `class LocalPotential : public virtual LocalPotential_Q, public virtual
LocalPotential_R {}`; GTH/HGH derive from the combined (closed forms in both spaces); a single-view
model derives from just its face.

**`SeparablePotential` — different ISP shape.** It *has* a neutral structural core
(`NumProjectors`/`Coefficient`/`AngularMomentum`, shared by both views); only the radial **leaf**
splits (`Projector(Z,p,q)` vs `BetaR(Z,p,r)`). So it's "shared structural base + per-view radial
leaf," not a three-way. Keep the two models asymmetric; don't force a uniform template.

### Module distribution (open; decide at impl)
One module per face so consumers import only what they touch:
- `qchem.Pseudopotential.LocalPotential.Reciprocal` → `LocalPotential_Q` (PW basis imports only this)
- `qchem.Pseudopotential.LocalPotential.Real` → `LocalPotential_R` (molecular basis imports only this)
- `qchem.Pseudopotential.LocalPotential` → combined + charge + concrete models (Ham term imports this)

Maximal BMI hygiene (PW never pulls in the real-space view, vice versa) vs module proliferation — the
tree already splits this fine, so it's in keeping.

## 3. Path A assembly = Becke-mesh quadrature (no libcint, no new integral type)

With the real-space view, the molecular PP block is **the exact quadrature the XC path already runs**:
```
⟨χᵢ|V_loc|χⱼ⟩ = Σ_g w_g χᵢ(r_g) χⱼ(r_g) Vloc(r_g)   // same shape as ⟨χᵢ|V_xc|χⱼ⟩
⟨χᵢ|β_p⟩      = Σ_g w_g χᵢ(r_g) β_p(r_g)             // overlap-with-a-function
```
Two things make it clean: (1) a good PP is **smooth** (erf local form is finite at r=0, no nuclear
cusp) — the mesh nails it, *better* than all-electron nuclear attraction; (2) the PP matrix is
**static** (density-independent) — quadratured **once**, not per SCF iteration. So Path A reuses the
Becke mesh + collocation already present for molecular DFT. libcint is **not** needed for Path A.

## 4. Path B = semilocal model + libcint

Add `SemilocalPotential` (`V_l(Z,r)` per channel + angular `P_l`) as a **molecular-only** model and a
`MakeSemilocalPotential` method (PW doesn't implement it — the mirror of `PseudoG0Energy` being
PW-only). The angular projector `P_l` is nonlocal in angle and awkward on a mesh — this is exactly
where libcint's analytic Type-2 ECP integral earns its keep. Note the satisfying symmetry of the final
capability picture: a **universal core** (local + KB-separable) plus **two honestly-asymmetric
extensions** — G0-alignment leans PW, semilocal leans molecular; neither is an orphan.

## 5. Term side — per-atom Coulomb-or-PP

On the molecular side the PP **replaces** bare nuclear attraction for pseudized atoms while
all-electron atoms keep full `Ven`. That's a *per-atom* choice — the same multi-species/`Zion` routing
already in `MultiSpecies_LocalPotential`/`MultiSpecies_SeparablePotential`, lifted into a unified
"external potential" term. First molecular appearance of the molecule/solid-agnostic external term.

## 6. Phasing & validation
- **A0:** eliminate `PseudoG0Energy` (verify Ω term-side); ISP-split `LocalPotential`/`SeparablePotential`
  with `LocalView<T>`; add the real-space views to GTH/HGH. *Gate: all PW tests bit-identical.*
- **A1:** molecular `Integrals_Pseudo<double>` via Becke-mesh quadrature; per-atom Coulomb-or-PP term.
  *Gate: a molecular GTH calc (e.g. a heavy-atom dimer) converges to a sane energy; all-electron atoms
  unaffected.*
- **B:** `SemilocalPotential` + libcint Type-2; consume a published ECP (def2-ECP). *Gate: vs a
  reference molecular-ECP energy.*

## 7. Key file pointers
- Capability: `src/Pseudopotential/Integrals_Pseudo.C` (shrink to 2 methods).
- Models: `src/Pseudopotential/LocalPotential.C` (3-facet ISP), `SeparablePotential.C` (core+leaf ISP).
- PW G0 impl to delete/relocate: `src/BasisSet/Lattice_3D/Imp/PlaneWave_IBS.C:250`.
- Molecular mesh/quadrature to reuse: the Becke mesh (`src/Mesh`) + the XC quadrature path
  (`FittedVxc`/`FittedVee` in `src/Hamiltonian/Internal/Imp/`).
- Term-side routing model: `MultiSpecies_*Potential` in the two model files.
