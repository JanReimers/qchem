// File: BasisSet/Lattice_3D/Imp/GPW_IBS.C  GPW_IBS implementation (ctors + identity).
module;
#include <cassert>
#include <cstdlib>   // std::getenv/std::atof (GPW_CDFIT_SCALE density-fit densification)
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

// The DFT fit-basis factory: a plane-wave fit basis over GPW's density grid (k=0 Gamma; the density is
// cell-periodic).  Both CD and Vxc share one grid this increment.
//
// DENSITY-FIT DENSIFICATION (doc/GPWPlan §0e step 2).  `DensityGrid()` (densityEcut = cutoffFactor*alpha_max)
// is scaled to resolve the ORBITAL exponent alpha_max -- but the density is the orbital PRODUCT (exponent
// ~2*alpha_max), which needs a FINER grid.  For a hard atom (F, alpha_max=40) the product is badly aliased on
// DensityGrid() (measured neg-frac 0.5, negCharge -9 e), collapsing the XC integral.  The fix is a denser fit
// grid HERE -- NOT a bigger densityEcut (which would coarsely inflate the orbital-scale grid, the local-PP
// sweep and the Hartree assembly too).  This is effective only because the fit-grid thread-through (758b92a8)
// makes MakeRepulsion3C/MakeOverlap3C build the tensor over the grid the fit basis carries.  GPW_CDFIT_SCALE
// (default 1 = DensityGrid(), bit-identical) sweeps the factor while we calibrate what a given alpha_max needs.
// (This subsumes the deferred GGA relCutoff densification: same seam, a different reason -- density resolution
// rather than the gradient -- so the mp.relCutoff<=1 GGA guard stays until v_xc's grad path also lands.)
//
// TODO (OOD cleanup, deferred): the {G}-grid POLICY (Gmax = cutoffFactor*alpha_max, and the density-product
// densification) should be STATED here, in CreateCDFitBasisSet, not read from GPW_Evaluator::DensityGrid() --
// the evaluator's job is low-level (integrals/IDs/op(r)), not deciding the charge-density fit grid.  The
// evaluator should expose alpha_max + the reciprocal lattice as data, and this function should build the grid
// from an explicit, readable policy ("4*alpha_max", "2x for the product", ...).
namespace
{
PW_Grid_Evaluator FitGrid(const PW_Grid_Evaluator& g0)
{
    double scale=1.0; if (const char* s=std::getenv("GPW_CDFIT_SCALE")) scale=std::atof(s);
    return scale>1.0 ? PW_Grid_Evaluator(g0.Recip(), rvec3_t(0,0,0), g0.Ecut()*scale) : g0;
}
} //anon
BasisSet::cFIT_CD_ABS* GPW_IBS::CreateCDFitBasisSet(const Structure*, const qcMesh::MeshParams& mp) const
{
    assert(mp.relCutoff<=1.0 && "GPW: relCutoff>1 (GGA fit-grid refinement) not wired; densityEcut is absolute");
    return new PlaneWaveFit_IBS(FitGrid(GPW_Evaluator::DensityGrid()), Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)));
}
BasisSet::cFIT_SF_ABS* GPW_IBS::CreateVxcFitBasisSet(const Structure*, const qcMesh::MeshParams& mp) const
{
    assert(mp.relCutoff<=1.0 && "GPW: relCutoff>1 (GGA fit-grid refinement) not wired; densityEcut is absolute");
    return new PlaneWaveFit_IBS(FitGrid(GPW_Evaluator::DensityGrid()), Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)));
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
