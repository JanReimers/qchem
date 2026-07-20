// File: BasisSet/Lattice_3D/Imp/GPW_IBS.C  GPW_IBS implementation (ctors + identity).
module;
#include <cassert>
#include <iostream>
#include <memory>
#include <string>

module qchem.BasisSet.Lattice_3D.GPW_IBS;
import qchem.Symmetry.Factory;              // BlochFactory (the convenience ctor + the k=0 fit-basis irrep)
import qchem.Symmetry.Lattice_3D.BlochQN;   // Symmetry::Lattice_3D::Getk (pry k out of the abstract Bloch irrep)
import qchem.BasisSet.Internal.DB_Cache;    // theCache<dcmplx>() -- process-wide cache for the static PP matrices
                                            // (qcLattice_BS is BasisSet-family, so it may peek at qcBasisSet Internal)
import qchem.BasisSet.Lattice_3D.Evaluators.PW;  // PW_Grid_Evaluator (the fit basis IS-A one; cross-cast target)

namespace qchem::BasisSet::Lattice_3D
{

GPW_IBS::GPW_IBS(const UnitCell& cell, const sym_t& irrep,
                 std::shared_ptr<const BasisSet::Real_BS> mol, double densityEcut, CellImages images,
                 double cutoffFactor)
    : BasisSet::IrrepBasisSetImp<dcmplx>(irrep)
    , GPW_Evaluator(std::move(mol), cell, densityEcut, Symmetry::Lattice_3D::Getk(irrep),
                    images==CellImages::HomeCellOnly, cutoffFactor) // irrep IS k
{}

// Convenience: build the Bloch irrep from BZ-grid indices and delegate to the primary constructor.
GPW_IBS::GPW_IBS(const UnitCell& cell, const ivec3_t& N, const ivec3_t& kIndex,
                 std::shared_ptr<const BasisSet::Real_BS> mol, double densityEcut, CellImages images,
                 double cutoffFactor)
    : GPW_IBS(cell, Symmetry::BlochFactory(N,kIndex), std::move(mol), densityEcut, images, cutoffFactor)
{}

// The DFT fit-basis factory: a plane-wave fit basis over GPW's OWN density grid (k=0 Gamma; the density is
// cell-periodic).  Both CD and Vxc share the one density grid this increment (relCutoff refinement deferred).
//
// GPW uses an ABSOLUTE densityEcut grid (a Gaussian basis has no single plane-wave orbital bandwidth to scale),
// so mp.relCutoff -- the CP2K REL_CUTOFF the Hamiltonian derives from the functional's GridCutoffFactor(), the
// GGA fit-grid densifier -- is NOT yet wired here (unlike PlaneWave_IBS, which builds its Vxc grid at
// Ecut*relCutoff).  For LDA relCutoff==1 so this is exact; GUARD it loudly so a future GGA-on-GPW attempt fails
// at this seam instead of silently quadraturing v_xc/grad(rho) on the LDA-grade grid.  Wiring it = build a
// separate, denser Vxc grid at densityEcut*relCutoff (mirroring the PW Vxc line).  See doc/GPWPlan.md.
BasisSet::cFIT_CD_ABS* GPW_IBS::CreateCDFitBasisSet(const Structure*, const qcMesh::MeshParams& mp) const
{
    assert(mp.relCutoff<=1.0 && "GPW: relCutoff>1 (GGA fit-grid refinement) not wired; densityEcut is absolute");
    return new PlaneWaveFit_IBS(GPW_Evaluator::DensityGrid(), Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)));
}
BasisSet::cFIT_SF_ABS* GPW_IBS::CreateVxcFitBasisSet(const Structure*, const qcMesh::MeshParams& mp) const
{
    assert(mp.relCutoff<=1.0 && "GPW: relCutoff>1 (GGA fit-grid refinement) not wired; densityEcut is absolute");
    return new PlaneWaveFit_IBS(GPW_Evaluator::DensityGrid(), Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)));
}

// The external-PP capability.  Local: G-space form-factor assembly (the model's FormFactor is used directly,
// no cross-cast -- mirrors the PW path).  Separable: the KB projectors need the real-space face, cross-cast.
// Both PP matrices are STATIC across an SCF but rebuilt PER k-BLOCK, so they go through the process-wide cache
// (theCache, keyed by BasisSetID + Structure::ID -- exactly the Nuclear() pattern): a multi-k / IBZ-vs-full-mesh
// run then reuses a k-block's PP across GPW_IBS instances instead of re-quadraturing it.  The build is the
// cache-miss `make` lambda; the outer Make* name is the Integrals_Pseudo override the term calls.
hmat_t<dcmplx> GPW_IBS::MakeLocalPotential(const Structure* cl, const Pseudopotential::LocalPotential& loc) const
{
    return theCache<dcmplx>().Get(IntegralsCache_Base::I2n::LocalPP, this, cl->ID(),
        [this,cl,&loc]{ return GPW_Evaluator::MakeLocalPP(cl, loc); });
}

// The CP2K local-PP split (doc/GPWPlan.md 0e-PP): the LONG (softened-Coulomb) matrix rides the smooth
// density-grid integrate-back (MakeLocalPPLong -- no sharp-field sweep); the SHORT (compact poly-Gaussian)
// matrix rides the sharp-field local-PP sweep (MakeLocalPP restricted to FormFactorShort).  Distinct cache
// keys keep them from colliding with each other or the full LocalPP.
hmat_t<dcmplx> GPW_IBS::MakeLocalPotentialLong(const Structure* cl, const Pseudopotential::LocalPotential& loc) const
{
    return theCache<dcmplx>().Get(IntegralsCache_Base::I2n::LocalPPLong, this, cl->ID(),
        [this,cl,&loc]{ return GPW_Evaluator::MakeLocalPPLong(cl, loc); });
}

hmat_t<dcmplx> GPW_IBS::MakeLocalPotentialShort(const Structure* cl, const Pseudopotential::LocalPotential& loc) const
{
    // NOTE (increment 2, doc/GPWPlan.md 0e-PP): the ANALYTIC short assembly (GPW_Evaluator::MakeLocalPPShort)
    // is built and finite-validated, but NOT yet wired into production -- it is exact, whereas the grid short
    // carries a ~0.5 Ha band-limiting error that CANCELS the grid long's (V_loc is smooth, so V_short(G) ~
    // -V_long(G) beyond the grid cutoff).  So analytic-short + grid-long misses the grid-calibrated gate;
    // both pieces must go analytic together.  Keep the grid short until the analytic long lands.
    return theCache<dcmplx>().Get(IntegralsCache_Base::I2n::LocalPPShort, this, cl->ID(),
        [this,cl,&loc]{ return GPW_Evaluator::MakeLocalPP(cl, loc, GPW_Evaluator::LocalPart::Short); });
}

hmat_t<dcmplx> GPW_IBS::MakeSeparablePotential(const Structure* cl, const Pseudopotential::SeparablePotential& nl) const
{
    auto* sepR=dynamic_cast<const Pseudopotential::SeparablePotential_R*>(&nl);
    assert(sepR && "GPW MakeSeparablePotential: the KB model must provide the real-space projector face (SeparablePotential_R)");
    return theCache<dcmplx>().Get(IntegralsCache_Base::I2n::SeparablePP, this, cl->ID(),
        [this,cl,sepR]{ return GPW_Evaluator::MakeSeparablePP(cl, *sepR); });
}

// The DFT 3-centre tables over the REQUESTED fit basis's grid (doc/GPWPlan §0e).  The fit basis \a c that the
// Hartree/XC term hands us is the one CreateCD/VxcFitBasisSet produced -- a PlaneWaveFit_IBS, which IS-A
// PW_Grid_Evaluator carrying the density-fit {G}/grid policy.  Cross-cast to that grid and build the tensor on
// it, so we RETURN THE REQUESTED TABLE rather than overriding the caller's fit-grid choice with the block's own
// (the shared EPW_Orbital_DFT_IBS mixin dropped \a c).  Bit-identical while the factory wraps DensityGrid();
// the seam is what lets the fit grid diverge (the deferred GGA Vxc densification) without touching these.
G_ERI3 GPW_IBS::MakeRepulsion3C(const cFIT_CD_ABS& c) const
{
    const auto& grid = dynamic_cast<const PW_Grid_Evaluator&>(c);   // throws bad_cast on a non-grid fit basis (loud)
    return GPW_Evaluator::Repulsion3CTensor(std::make_shared<const PW_Grid_Evaluator>(grid));
}
G_ERI3 GPW_IBS::MakeOverlap3C(const cFIT_SF_ABS& c) const
{
    const auto& grid = dynamic_cast<const PW_Grid_Evaluator&>(c);
    return GPW_Evaluator::Overlap3CTensor(std::make_shared<const PW_Grid_Evaluator>(grid));
}

std::string GPW_IBS::BasisSetID() const
{
    return Name()+GPW_Evaluator::IDFragment();   // Name + "|mol=..|k=..|cell=..|dEcut=.."
}

std::ostream& GPW_IBS::Write(std::ostream& os) const
{
    return os << Name() << " IBS: " << GetNumFunctions() << " periodic Gaussians, " << GetSymmetry();
}

} //namespace
