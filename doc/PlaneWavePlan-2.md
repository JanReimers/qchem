# Code Review 
Following the pipline Si SCF DFT test, line 1016 in UnitTests/PlaneWaveDFTUT.C

## UnitTests/PlaneWaveDFTUT.C

- Line 1020: Amat construction not user friendly.  Should read: UnitCell cell=Cubic(h); //Off other types Orthorhombic(h1,h2,h3), Triclinc(...) etc.
- Line 1021: Lattice_3D        lat(cell, ivec3_t(1,1,1)); The ivec3_t(1,1,1) tells me there is only one k point in the BZ integration mesh!!  This supposed to be somethng like ivec3_t(10,10,10) for 10x10x10 inegration mesh.
- Line 1023: The irrep basis set should come from a factory, and be held as an abstract Orbtial_1E_IBS*  If that breaks the code below then we need to fix the interfaces, not cast the pointer.
- Line 1027: Evenutally we ned support fractional coordinates and the type system should be able to distinguish dimensionless, a.u. (and later Angtsrom.)
- Line 1028 to 1030: THis looks a job for the factory or PlaneWave_IBS constructor.  The structure (lat) should have all the nessecary information.
- Line 1036: PW_BasisSet should be defined in src/BasisSet/Lattice_3D/BasisSet.C and should come from a factory, and held as a BasisSet* in the unit test code.  The PW_BasisSet consructor should have enough information available to build a list of PlaneWave_IBS instances.
- Lines 1040-1045: This all belongs in src/Hamiltonian/Internal/Imp/Hamiltonians.C, and a new Ham type needs to be declared in src/Hamiltonian/Internal/Hamiltonians.C
- Lines 1051-1053: THis does not belong in a high level SCF integration test, but right now I don't know where to move it.

## src/BasisSet/Lattice_3D/PlaneWave_IBS.C
- Line 40: ivec3_t& N (feeds itsN) is not used other than for console output.
- Line 40: kIndex, on the atom side you see the all IBSs are constructed from a sym_t& irrep.  And then inside the constructor we call a helper (Getl(irrep)) in qcSymmtery to pry specific info out of the abstract reference. So finally the ideal construct sig is:  PlaneWave_IBS(const ReciprocalLattice& recip, const sym_t& irrep, double Ecut);
The BasisSet constructor Takes a Lattice_3D structure, make a reciprocal lattice and then figures the required list of BZ k vectors. (In future when we get space groups these can be trimmed and weighted). Ideally, this should be the only k loop in the whole framework (I think, we shall see:)
- Line 47,69,75,85: Who consumes these?  It looks like private or protected operations.  Remember local (non integration) unit tests should have friend access.
SOLID ISP: these four smell like a separate interface consumed by something other than the Hamiltonian terms.
- Line 57: This looks like a construction time thing.  At a minimum it should be protected and only called by derived classes.

## Ideal psuedo code for the Si SCF DFT test

This is just defining a goal ... some of this may not work in the end ... but I looks really nice and clean to me.  Even a non programmer can adjust this and have fun.  I left the rest unspecified, this is a good chunk to think about for now.

```c++
{
    using namespace qchem::BasisSet::Lattice_3D;
    const double a=10.26; // a.u.
    CubicUnitCell cell(0.5*a); // 0.5* is for ???
    cell.AddAtom(14,{0,0,0});
    cell.AddAtom(14,{.25,.25,.25}); //Fractional coordinates.
    Lattice_3D lat(cell,{10,10,10}); //define BZ k-grid
    BasisSet* bs=Factory(Type::PW,lat); //Polymorphic so we are forced to use a pointer.
    const int  Nelec=8;
    Crystal_EC  ec(bs, Nelec); //EC constructor can loop over bs->GetIrreps()
    using enum qchem::Hamiltonian::Type;
    cHamiltonianImp* ham=qchem::Hamiltonian::Factory(MyFancyPseudoXYZ,lat);
    .
    .
    .

    // Be carefull what we delete in case the SCFiterator owns bs or ham ... I can't remember!!

}
```
# Higher level goal(s)

We need to copy src/BasisSet/Molecule/Evaluators/Evaluator.C over to src/BasisSet/lattice_3D/Evaluators/Evaluator.C and adjust the namespace and module name accordingly. Then we are free to (grudginly) tune these interfaces to meet the requirements for all the Lattie_3D IBS classes.

In src/BasisSet/Molecule/IrrepBasisSet.C we can see that the Orbital_1E_IBS class supports both single element evaluators (is1E_Evaluator)  and full matrix Evaluators (isM_1E_Evaluator) for Make{Overlap,MakeKinetic, MakeNuclear}.  But the Molecule example is either all single element or all matrix.  For PW it appears we have a mixture.  Can we function foward inline single element calls for all three Make{Overlap,MakeKinetic, MakeNuclear}?  We can still capture the diagonal flavour for Overap and Kinetic in Orbital_1E_IBS<is1E_Evaluator>. SOmething like:

```c++
template <is1E_Evaluator E> class Orbital_1E_IBS<E>  : public virtual ::BasisSet::Orbital_1E_IBS<dcmplx>
{
    ...
    chmat_t MakeOverlap() const
    {
        const E& e=Cast();
        chmat_t S=blazem::zeroH<dcmplx>(e.size());   // hmat_t(n) does NOT zero the off-diagonals
        for (size_t i:e.indices()) S(i,i)=e.MakeOverlap(i,i);
        return S;
    }
... //etc.
};
```
Start migrating the PlaneWave_IBS class so that all code for the is1E_Evalualtor concept is in there. And PlaneWave_IBS inherits PlaneWave::Evaluator  so that the const E& e=Cast(); line above works.

See src/BasisSet/Molecule/Evaluators for examples.

At the end of each evaluator class we should have a statement like:
```c++
static_assert(is1E_Evalualtor<Evaluator>);
```
If the is1E_Evalualtor is impossible and we must use isM_1E_Evalualtor that should also be fine.  The goal is put the least possible burden on Evaluator developer, which means is1E_Evalualtor is prefered.  Then our collection of IBS basis sets is maximally extensible.

Currently all the pseudo potential code lives in the PlaneWave_IBS class. Naievly I would expect there to be one or more Vpsuedo classs in src/Hamiltonian/Internal/Terms.C where all the pseudo potential construction code should live.  But the IBS knows how the generate the matrices <a|XXX|b> required for Vpsuedo to build itself.  This is the tought part, putting a wall where IBS knows as little as possible about what type of Vpsuedo we are using and Vpsuedo is not aware that the PW IBS has a bunch of G vectors inside.  So what is XXX?  For atoms XXX is just 1/r very general, and molecules we forced to be more specific, XXX is Sum(-Zi/ri).  For pseudos it is presumably more complicated.

# Dangling (put a pin in it) items.

## Plane-wave / T-templating pins (accumulated 2026-06-24, after Stage-3 commit 975aa61e)

Grouped by project-module (internal library). These are the pinned/hanging items from the
WaveFunction->SCFIterator complexification milestone. Memory: [[project_plane_wave_plan]].

### qcHamiltonian
- **bare->r\* rename** (remove transitional aliases): bare names (`Static_HT`, `Dynamic_HT`,
  `Hamiltonian`, `Static_CC`, `Dynamic_CC`, `DM_CD`, `Composite_CD`, `obs_t`, ...) are still
  `using X = rX` shims so ~200 existing sites didn't churn. Do the codebase-wide rename. NOTE the
  `rHamiltonian` rename must also dissolve the `qchem::Hamiltonian` namespace-vs-type shadow
  (the `::Hamiltonian` / `qchem::Hamiltonian::` qualifications scattered around are working around it).
- **newCD double-duty redesign**: `Dynamic_HT_Imp` caches by Irrep and `newCD(cd)` does two jobs
  (cache-clear + refit-trigger). The cache-freshness bug was fixed (commit 582f1d05: GetMatrix clears
  on cd-change), but the deeper smell remains — split the two responsibilities. See
  [[project_hamiltonian_dynamic_cache_bug]].

### qcChargeDensity
- **DM_CD interface-slimming**: the dcmplx plane-wave density NA-asserts the HF/fit-only methods
  (`AccumulateDirect`/`AccumulateExchange`, `GetRepulsion3C`, `FitGetConstraint`). These are real
  Gaussian-basis facilities bolted onto the shared `tDM_CD<T>` base (`AccumulateDirect/Exchange`
  pure-virtual) + the inherited `DensityFFClient`. Carve the HF/fit contract off the core density
  interface so the complex path doesn't have to stub them.

### qcBasisSet
- **Complex integral cache**: `theGlobalCache` is `IntegralsCache<double>*` only. `Integrals_Overlap<dcmplx>::Overlap()`
  currently uses a per-instance `static map<this, hmat_t<dcmplx>>` stand-in buffer (keeps the
  no-data-in-interface diamond rule). Build a real complex `IntegralsCache` (a `theGlobalCacheC`, or
  make the cache T-aware) and route the complex Overlap/Kinetic/etc. through it. Also revisit the
  `BasisSet<dcmplx>::Create{CD,Vxc}FitBasisSet` NA-asserts when fitting reaches the complex path (it won't soon).
- **DFTPotential_IBS<T>** provisional name — confirm/rename the matrix-delivery G-space-potential interface.

### qcMesh / qcLattice_BS
- **qcMesh owns quadrature**: PlaneWave_IBS's hand-rolled `Integral(f)=Sum Omega/Npts*f` + UniformGrid
  belongs in qcMesh (Mesh carries points+weights; integral = Sum w_i f(r_i)). Then migrate all PW IBSs.
- **Reconcile mesh ownership**: NEW model = the basis owns the mesh (correct); OLD code has fitted
  potentials owning the Mesh (likely wrong) — revisit after the qcMesh refactor.

### qcSCFAccelerator
- **Complex DIIS/GDM accelerator**: only `tSCFAcceleratorNull<dcmplx>` exists for the complex path. With
  no DIIS damping, the density matrix drifts marginally within a degenerate occupied manifold (Si Gamma_25';
  energy flat to 1e-9 but |delta rho| won't reach 1e-7). A complex DIIS would fix the drift and tighten convergence.

### qcElConfig
- **Crystal_EC as a module**: currently a test-local struct in PlaneWaveDFTUT.C (GetN(Bloch-irrep)=Nval,
  UsesAufbau=false). Promote to a real `src/ElectronConfigurations/Crystal_EC.C` module.

### qcWaveFunction / qcLattice_BS (deferred milestone — NOT a pin, the next big step)
- **Stage 4: k = Bloch irrep**. Make the WaveFunction per-irrep loop range over k-points so the BZ
  Sum_k w_k IS the framework's irrep loop, reproducing the 2x2x2 Si result (gap ~0.71 eV) through the
  real SCFIterator. The single-k path is done; this generalizes it to multi-k.

### Deferred PW physics (out of scope of the integration)
- FFT swap (prototype direct G-space sums are fine), Ewald ion-ion / absolute total energy,
  symmetry-reduced k-mesh, polarized LSDA (spin-polarized plane-wave path).

### qcStructure
- We need to change the typedef cl_t->st_t  (cl_t a left over from when Structure was called Cluster). Same for valriable names like cl or cluster
  
- Structure is going to need something like 
  bool isFinite() const=0;
  For Lattice_3D this returns false.  
  The Vnn Hamiltonian term will use this to decide if and Ewald sum is rquired.

### qcFitting
- Look at the typedefs bs_t is too vague. If it is a Fit_IBS then just say so, call it fitbs_t.