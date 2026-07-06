# qcMesh1 — clean-room design sketch

A from-scratch replacement for qcMesh.  Rationale (see qcMeshUpgradePlan.md): the cruft is
*structural* (storage model, ISP-violating function interfaces, integrator-as-a-class, dead code),
while the *numerics* (radial/angular grids, Becke, Gauss-Legendre) are sound.  When the damage is
scaffolding and the value is transplantable kernels, a clean room beats an in-place refactor.

**THE ONE DISCIPLINE:** transplant the numerical kernels VERBATIM (radial log-grid + stiffness handling,
angular schemes, GL nodes).  If you find yourself re-deriving a quadrature formula, stop and copy it.
The DFT oracle + pinned molecular-energy anchors are the regression backstop.

## Scope

- **IN (geometry-free primitives):** Mesh (points+weights), 1D radial, angular, single-center product,
  free-function quadrature, the pointwise function interfaces, shared 1D Gauss-Legendre.
- **OUT (geometry-aware → stays a Structure-side consumer):** Becke molecular mesh, the future
  FLAPW cell mesh.  These need atom positions; they CONSUME qcMesh1, they aren't part of it.
  (A companion from-scratch Becke rewrite is sketched at the end.)

## Do NOT carry over (the old qcMesh's "what NOT to do")

- R/W stored as `vector<tuple<r,w>>` (AoS) — invert to separate arrays (SoA).
- `ScalarFunction`/`VectorFunction` carrying `operator()(const Mesh&)` — ISP violation; drags Mesh into
  every pointwise field.  Split it out.
- `MeshIntegrator` as a class holding a `Mesh*` — make these free functions.
- `Mesh::Clone()` — not needed once Mesh is a value type (verify no polymorphic-copy caller first).
- The `Repulsion`/`Repulsion3C` (1/r12) integrals — they "can't work" (Coulomb goes through Hartree/
  Poisson, not quadrature).  Drop them (verify no live caller first).
- The complex `MeshParams` requiring unused fields — typed struct with defaults instead.

## Core interfaces

```cpp
// ---- Mesh: a concrete VALUE type = quadrature points + weights, stored SEPARATELY (SoA). -----------
// Points are streamed by evaluators (phi(r)); weights by integrators.  No polymorphism, no Clone.
class Mesh
{
public:
    Mesh() = default;
    const vec_t<rvec3_t>& Points () const { return itsR; }   // r_i      (pick the real vec-of-vec3 typedef)
    const rvec_t&         Weights() const { return itsW; }   // w_i
    size_t                size   () const { return itsW.size(); }

    void Append(const rvec3_t& r, double w);                 // builders push points
    void ShiftOrigin(const rvec3_t& o);                      // r_i += o  (place a center)
private:
    vec_t<rvec3_t> itsR;
    rvec_t         itsW;
};

// ---- Pointwise fields: NO Mesh dependency (that was the ISP sin). ----------------------------------
template <class T> class ScalarField                         // e.g. rho(r), vxc(r), -Z/|r-R|
{
public:
    virtual ~ScalarField() = default;
    virtual T         operator()(const rvec3_t&) const = 0;
    virtual vec3_t<T> Gradient  (const rvec3_t&) const = 0;
};

template <class T> class BasisField                          // the basis as a vector field [phi_i(r)]
{
public:
    virtual ~BasisField() = default;
    virtual vec_t<T>     operator()(const rvec3_t&) const = 0;   // [ phi_i(r) ]
    virtual vec3vec_t<T> Gradient  (const rvec3_t&) const = 0;   // [ grad phi_i(r) ]
};
```

## Quadrature — free functions, and the integral-type vocabulary

`MeshIntegrator`-the-class becomes plain functions taking `const Mesh&` first.  Crucially, the physics
stays with the CALLER (it supplies the `ScalarField` V); qcMesh1 only knows `Sum_i w_i (...)`.  Note how
`Inv_r1`/`Inv_r2`/`IntegralPotential` all COLLAPSE into one `WeightedOverlap(mesh, basis, V)` — the
caller passes V = 1/r, 1/r^2, or vxc.  That small, physics-free set IS the integral-type taxonomy the
pseudo-wall leans on.

```cpp
// integral f d3r = Sum_i w_i f(r_i)
template <class T> T Integrate(const Mesh&, const ScalarField<T>&);

// <a_i | a_j>  and  <a_i | b_j>
template <class T> smat_t<T> Overlap(const Mesh&, const BasisField<T>&);
template <class T> mat_t <T> Overlap(const Mesh&, const BasisField<T>& a, const BasisField<T>& b);

// <a_i | V | a_j>  -- subsumes Inv_r1 (V=1/r), Inv_r2 (V=1/r^2), and the DFT IntegralPotential (V=vxc).
template <class T> hmat_t<T> WeightedOverlap(const Mesh&, const BasisField<T>&, const ScalarField<double>& V);

// <grad a_i | grad a_j>  (the kinetic <p^2> block) and the gradient cross terms if needed.
template <class T> hmat_t<T> KineticGrad2(const Mesh&, const BasisField<T>&);
template <class T> mat_t <T> Grada_b     (const Mesh&, const BasisField<T>& a, const BasisField<T>& b);
```

## Radial / Angular / single-center product

```cpp
class RadialMesh    { public: virtual const rvec_t& R() const=0; virtual const rvec_t& W() const=0; };
class AngularMesh   { public: virtual const vec_t<rvec3_t>& Dirs() const=0; virtual const rvec_t& W() const=0; }; // sum W = 4pi
// concretes (TRANSPLANTED numerics): MHL/Log/Linear radial; Gauss/GaussLegendre/EulerMaclaren angular (Lebedev later)

// Single-center tensor product -> Mesh.  GEOMETRY-FREE (no atom positions): r_i*Omega_j, w = w_i^rad * w_j^ang.
Mesh ProductMesh(const RadialMesh&, const AngularMesh&);
```

## Shared 1D Gauss-Legendre

One implementation, used by the radial meshes, the PW real-space grid, AND the Atom BSpline integrator
(today each rolls its own; the BSpline one is `gauleg.f` from Numerical Recipes).  Port gauleg -> C++/blaze;
acceptance test = nodes/weights vs a known table + exact integration of polynomials up to degree 2n-1.

```cpp
struct GaussLegendre { rvec_t x, w; GaussLegendre(int n, double a, double b); };
```

## Typed MeshParams (NOT json) + the XC-grid override

Typed + defaulted so the user sets only what they need (C++20 designated initializers, compiler-checked):

```cpp
enum class RadialKind  { MHL, Log, Linear };
enum class AngularKind { Gauss, GaussLegendre, EulerMaclaren };
struct MeshParams
{
    RadialKind  radial   = RadialKind::MHL;   int nRadial  = 30;
    AngularKind angular  = AngularKind::Gauss; int nAngular = 12;
    // ...only the LIVE knobs; nothing the user must fill in but never uses.
};
// PW real-space XC grid (qcMeshUpgradePlan item 3): the basis builds its grid from a MeshParams whose
// DEFAULT is the current AutoGrid (4*m_max+1) -> override is opt-in, no breaking change.
```

## Companion: from-scratch Becke molecular mesh (Structure side)

A clean Becke rewrite is likely faster than debugging the current MoleculeMesh bug (coincident-atom
dimer != atom).  Lives in qcStructure (needs atom positions); CONSUMES qcMesh1.

```cpp
// qcStructure: Becke-partitioned multi-center integration mesh.
Mesh MakeMolecularMesh(const Molecule&, const MeshParams&);
//   for each atom a:  Mesh m = ProductMesh(radial, angular);  m.ShiftOrigin(R_a);
//                     reweight each point by the Becke cell function  p_a(r) / Sum_b p_b(r)
//                     append to the combined mesh
```

- **Likely bug culprit:** the confocal-elliptical coordinate mu_ab = (|r-R_a| - |r-R_b|)/R_ab is 0/0 when
  R_ab -> 0 (coincident atoms).  Handle the degenerate/near-coincident pair explicitly.
- **Acceptance test (the bug's mirror):** a homonuclear dimer at separation -> 0 must integrate to the
  single-atom result; and a well-separated dimer's total must equal 2x the atom (additivity).  Plus the
  existing molecular-DFT energy anchors (the oracle) gate it.

## Migration order (incremental, green-gated, oracle-backstopped)

1. Stand up qcMesh1: `Mesh` (SoA) + one RadialxAngular + `ProductMesh` + `Integrate` + GL.  Smoke test.
2. Migrate the PW IBSs first (cleanest consumer): hand-rolled `Integral`/`UniformGrid` -> qcMesh1.
3. From-scratch Becke `MakeMolecularMesh` + the dimer test; migrate molecular DFT to it.
4. Migrate the fitted potentials (the gnarly PW#2 "big job") -- do LAST, deliberately, oracle-gated.
5. Delete old qcMesh once nothing imports it.

Keep milestone 1 TINY and prove it green before widening -- the rewrite-stall guard.
```
