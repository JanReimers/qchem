// File: BasisSet/Gaussian/Lattice/Imp/GPW_IBS.C  GPW_IBS implementation (ctors + identity).
module;
#include <cassert>
#include <iostream>
#include <map>       // MakeProjectorMatrixByL (the per-l KB diagnostic)
#include <memory>
#include <stdexcept>   // std::invalid_argument (the transition pair's bra check)
#include <string>
#include <vector>
#include <cmath>       // std::round (the little-group test, S3)

module qchem.BasisSet.Gaussian.Lattice.GPW_IBS;
import qchem.BasisSet.Lattice_IBS;          // Lattice::ToScalar (the Step-3 exact TRIM narrow; identity for dcmplx)
import qchem.Symmetry.Factory;              // BlochFactory (the convenience ctor + the k=0 fit-basis irrep)
import qchem.Symmetry.Lattice_3D.BlochQN;   // Symmetry::Lattice_3D::Getk (pry k out of the abstract Bloch irrep)
import qchem.Symmetry.Lattice_3D.BlochOperationRep; // BuildBlochOperationRep (OrbitalImage, KAxisPlan P4)
import qchem.Symmetry.Lattice_3D.SpaceGroup; // DirectOp {W|τ} (the ctor's IBZ raster ops param type)
import qchem.BasisSet.Internal.DB_Cache;    // theCache<dcmplx>() -- process-wide cache for the static PP matrices
                                            // (qcMolecule_BS is BasisSet-family, so it may peek at qcBasisSet Internal)
import qchem.BasisSet.PlaneWave.Evaluators;  // PlaneWave::PW_Grid_Evaluator (the fit basis IS-A one; cross-cast target)
import qchem.BasisSet.PlaneWave.PlaneWaveFit_IBS;

namespace qchem::BasisSet::Gaussian
{

template <class T> tGPW_IBS<T>::tGPW_IBS(const UnitCell& cell, const sym_t& irrep,
                 std::shared_ptr<const BasisSet::Real_BS> mol, std::shared_ptr<const GpwRasterLadder> ladder, CellImages images,
                 std::vector<Symmetry::Lattice_3D::DirectOp> directOps, const GPWTolerances& tol)
    : BasisSet::IrrepBasisSetImp<T>(irrep)
    , GPW_Evaluator(std::move(mol), cell, std::move(ladder), Symmetry::Lattice_3D::Getk(irrep), irrep->IsReal(),
                    images==CellImages::HomeCellOnly, tol) // irrep IS k; IsReal() = the exact TRIM fact (doc/RealComplexPlan.md Step 1)
{
    SetSymmetryOps(std::move(directOps));   // the evaluator's T1 {G} fold
}

// A bare block: its own ladder from the same builder the factory uses, then the primary constructor.
template <class T> tGPW_IBS<T>::tGPW_IBS(const UnitCell& cell, const sym_t& irrep,
                 std::shared_ptr<const BasisSet::Real_BS> mol, double densityEcut, CellImages images,
                 double cutoffFactor, double ladderFactor, RasterFields rasterFields, const GPWTolerances& tol)
    : tGPW_IBS(cell, irrep, mol,
              MakeRasterLadder(cell, *mol, {.densityEcut=densityEcut, .cutoffFactor=cutoffFactor, .ladderFactor=ladderFactor,
                                            .rasterFields=rasterFields, .tol=tol}),
              images, {}, tol)
{}

// Convenience: build the Bloch irrep from BZ-grid indices and delegate to the bare-block constructor.
template <class T> tGPW_IBS<T>::tGPW_IBS(const UnitCell& cell, const ivec3_t& N, const ivec3_t& kIndex,
                 std::shared_ptr<const BasisSet::Real_BS> mol, double densityEcut, CellImages images,
                 double cutoffFactor, double ladderFactor, RasterFields rasterFields, const GPWTolerances& tol)
    : tGPW_IBS(cell, Symmetry::BlochFactory(N,kIndex), std::move(mol), densityEcut, images, cutoffFactor, ladderFactor, rasterFields, tol)
{}

template <class T> mat_t<T> tGPW_IBS<T>::MakeRepMatrix(const Symmetry::SymOp& g) const
{
    const rvec3_t k =Symmetry::Lattice_3D::Getk(this->GetSymmetry());
    const rvec3_t d =Symmetry::Lattice_3D::ImageOfK(g, k)-k;
    if (std::fabs(d.x-std::round(d.x))+std::fabs(d.y-std::round(d.y))+std::fabs(d.z-std::round(d.z))>1e-10)
        throw std::logic_error("tGPW_IBS::MakeRepMatrix: the op maps k onto another k -- it is not in this block's little group");
    // M in the block's OWN functions: its shells and its Bloch phases in the form they were built with (exact ±1 at a TRIM k,
    // the irrep's exact fact) -- so a real TRIM block's M is real, and the project's ONE exact narrow takes it to T.
    try
    {
        return Lattice::ToScalar<T>(Symmetry::Lattice_3D::BuildBlochOperationRep(GetAoShells(), GPW_Evaluator::Cell().GetCellMatrix(),
                                                                                  g, k, 1e-6, this->IsReal()));
    }
    catch (const std::runtime_error& e) {throw std::logic_error(std::string("tGPW_IBS::MakeRepMatrix: ")+e.what());}
}

template <class T> mat_t<dcmplx> tGPW_IBS<T>::OrbitalImage(const Symmetry::SymOp& g) const
{
    return Symmetry::Lattice_3D::BuildBlochOperationRep(GetAoShells(), GPW_Evaluator::Cell().GetCellMatrix(), g,
                                                        Symmetry::Lattice_3D::Getk(this->GetSymmetry()));
}

// THE DFT FIT BASES -- the {G} sets the density (CD) and the XC potential (Vxc) are represented / quadratured on.
// There is exactly ONE density FFT grid in GPW: GPW_Evaluator::DensityGrid() (== itsFFT_R_G_Grids), the {r}<->{G}
// FFT pair whose cutoff is densityEcut = cutoffFactor*alpha_max (the ctor's DENSITY-scale policy: cutoffFactor is
// sized to resolve the orbital PRODUCT exp(-2*alpha_max r^2), NOT a single orbital -- GPW orbitals are analytic
// Gaussians with no grid at all).  Both fit bases are built over that one grid:
//
//   CD  fit basis {G}_rho  (density representation): the density rho = Sum D_ij chi_i chi_j is expanded here, so
//        {G}_rho must resolve the product 2*alpha_max -- which IS what cutoffFactor calibrates.  {G}_rho =
//        DensityGrid().  (Under-resolving it aliases rho into large spurious negative lobes -> the XC collapse,
//        doc/GPWPlan §0e step 2: F needed ~8*alpha_max for negCharge -9 e -> -0.03 e.)
//   Vxc fit basis {G}_vxc  (XC quadrature): v_xc(rho) is sampled here.  {G}_vxc = relCutoff * {G}_rho, since the
//        nonlinearity rho^{4/3} carries more bandwidth than rho.  LDA: relCutoff == 1 (mild) -> {G}_vxc =
//        {G}_rho = DensityGrid().  A GGA (grad rho) sets relCutoff > 1 -> a denser Vxc grid; GUARDED (assert)
//        until that path is wired.
//
// So a reader sees BOTH bases here, and for LDA both are the one density grid.  (History: a temporary
// GPW_CDFIT_SCALE knob that forked a SECOND, denser grid is RETIRED -- resolving the product is cutoffFactor's
// job, one grid.  doc/GPWPlan §0e step 2.)
template <class T> BasisSet::cFIT_CD_ABS* tGPW_IBS<T>::CreateCDFitBasisSet(const Structure*, const qcMesh::MeshParams&) const
{
    // {G}_rho = DensityGrid() (cutoffFactor*alpha_max, resolving the density product); no relCutoff on the CD grid.
    return new PlaneWave::PlaneWaveFit_IBS(GPW_Evaluator::DensityGrid(), Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)),
                                "densityFit");
}
template <class T> BasisSet::cFIT_SF_ABS* tGPW_IBS<T>::CreateVxcFitBasisSet(const Structure*, const qcMesh::MeshParams& mp) const
{
    // A BLOCK builds only its OWN fitted representation -- here the plane-wave {G} one.  Choosing BETWEEN
    // representations (delta vs fitted) is the WHOLE-SET factory's job, one level up in tBasisSet, since
    // the delta basis is a property of the run's quadrature and not of any one Bloch block.
    // {G}_vxc = relCutoff * {G}_rho.  LDA relCutoff==1 => == DensityGrid(); a GGA's denser grid is not wired yet.
    assert(mp.relCutoff<=1.0 && "GPW: relCutoff>1 (GGA denser Vxc grid) not wired -- the LDA Vxc grid = the CD grid");
    // The SAME basis as CreateCDFitBasisSet's (same {G} ball, same BasisSetID): since 2026-11-09 it stores no crystal ops -- the
    // run's group star-averages the raster through G_RasterStarField -- so the two instances are identical.
    // ROLE = "vxcFitGrid", NOT "xcQuadrature" (2026-08-20).  This object is the PLANE-WAVE Vxc fit grid;
    // whether it is the run's XC QUADRATURE is the Hamiltonian's route decision, which a basis factory
    // cannot see.  It used to announce itself as "xcQuadrature" unconditionally, so a Becke run printed
    // `grids.xcQuadrature kind Becke` from the route and then `grids.xcQuadrature kind Uniform` from here
    // -- the same report key answering the reader's question two different ways on one run.  The route now
    // owns that key exclusively (Hamiltonians.C, both branches).
    return new PlaneWave::PlaneWaveFit_IBS(GPW_Evaluator::DensityGrid(), Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)),
                                "vxcFitGrid");
}

// The species-field integral service.  Both matrices are STATIC across an SCF but rebuilt PER k-BLOCK, so
// they go through the process-wide cache (theCache, keyed by BasisSetID + Structure::ID -- exactly the
// Nuclear() pattern): a multi-k / IBZ-vs-full-mesh run then reuses a k-block's matrices across GPW_IBS
// instances instead of re-quadraturing them.  The build is the cache-miss `make` lambda; the outer Make* is
// the Orbital_PP_IBS override the term calls.
//
// The RANGE picks the route AND the cache key (doc/GPWPlan.md 0e-PP): the LONG (softened-Coulomb) matrix rides
// the smooth density-grid integrate-back (MakeLocalPPLong -- no sharp-field sweep); the SHORT matrix is the
// ANALYTIC 3-centre Gaussian lattice sum (LatticeSum1E::MakeLocalGaussian, no grid; step (b), 2026-07-22 --
// safe to wire ONLY because step (a) first made the grid LONG standalone-exact, the absolute kappa rule with
// e^{-kappa/2} pair tails: the old grid short's ~0.5 Ha band-limit error used to CANCEL the grid long's, so
// exact-short + sloppy-long missed the gate.  Cross-validated by GPW.LocalPPKappaSelfConverged; the grid sweep
// remains the fallback inside MakeLocalPPShort for a field without the closed-Gaussian face).  Full is the
// plain sweep over the whole field.  Distinct keys keep the three from colliding.
template <class T> hmat_t<T> tGPW_IBS<T>::MakeSpeciesFieldMatrix(const Structure* cl, const SpeciesRadialField& f, FieldRange rng) const
{
    switch (rng)
    {
    case FieldRange::Long:
        return theCache<T>().Get(IntegralsCache_Base::I2n::LocalPPLong, this, cl->ID(),
            [this,cl,&f]{ return Lattice::ToScalar<T>(GPW_Evaluator::MakeLocalPPLong(cl, f)); });
    case FieldRange::Short:
    {
        // Its own key, without the grid facts (GPW_Evaluator::LocalPPShortIDFragment): a rebuilt basis reuses the block.
        struct ShortKey : DBCacheClient
        {
            std::string id; size_t dim;
            ShortKey(std::string i, size_t d) : id(std::move(i)), dim(d) {}
            std::string BasisSetID() const override {return id;}
            size_t      CacheDim  () const override {return dim;}
        } key(Name()+GPW_Evaluator::LocalPPShortIDFragment(f), this->CacheDim());
        return theCache<T>().Get(IntegralsCache_Base::I2n::LocalPPShort, &key, cl->ID(),
            [this,cl,&f]{ return Lattice::ToScalar<T>(GPW_Evaluator::MakeLocalPPShort(cl, f)); });
    }
    default:
        return theCache<T>().Get(IntegralsCache_Base::I2n::LocalPP, this, cl->ID(),
            [this,cl,&f]{ return Lattice::ToScalar<T>(GPW_Evaluator::MakeLocalPP(cl, f)); });
    }
}

template <class T> hmat_t<T> tGPW_IBS<T>::MakeProjectorMatrix(const Structure* cl, const SpeciesProjectorSet& nl) const
{
    auto* sepR=dynamic_cast<const SpeciesProjectorSet_R*>(&nl);
    assert(sepR && "GPW MakeProjectorMatrix: the projector set must provide the real-space radial face (SpeciesProjectorSet_R)");
    return theCache<T>().Get(IntegralsCache_Base::I2n::SeparablePP, this, cl->ID(),
        [this,cl,sepR]{ return Lattice::ToScalar<T>(GPW_Evaluator::MakeSeparablePP(cl, *sepR)); });
}

template <class T> std::map<int,hmat_t<T>> tGPW_IBS<T>::MakeProjectorMatrixByL(const Structure* cl, const SpeciesProjectorSet& nl) const
{
    // Diagnostic face (doc/SphericalLatticePlan.md I0): built on demand, NOT DB-cached -- its consumer
    // (the Ven_PP_NonLocal per-l energy print) caches the result term-side for the run's lifetime.
    auto* sepR=dynamic_cast<const SpeciesProjectorSet_R*>(&nl);
    assert(sepR && "GPW MakeProjectorMatrixByL: the projector set must provide the real-space radial face (SpeciesProjectorSet_R)");
    std::map<int,hmat_t<T>> out;
    for (auto& [l,m] : GPW_Evaluator::MakeSeparablePPByL(cl, *sepR)) out.emplace(l, Lattice::ToScalar<T>(m));
    return out;
}

template <class T> mat_t<T> tGPW_IBS<T>::MakeAtomCentredOverlap(const BasisSet::AtomCentredSet& set) const
{
    const auto* g=dynamic_cast<const BasisSet::AtomCentredSet_Gaussian*>(&set);
    if (!g) throw std::invalid_argument("GPW MakeAtomCentredOverlap: the set must carry the closed-Gaussian face (AtomCentredSet_Gaussian)");
    return Lattice::ToScalar<T>(GPW_Evaluator::MakeAtomCentredOverlap(*g));
}
template <class T> hmat_t<T> tGPW_IBS<T>::MakeSetOverlap(const BasisSet::AtomCentredSet& set) const
{
    const auto* g=dynamic_cast<const BasisSet::AtomCentredSet_Gaussian*>(&set);
    if (!g) throw std::invalid_argument("GPW MakeSetOverlap: the set must carry the closed-Gaussian face (AtomCentredSet_Gaussian)");
    return Lattice::ToScalar<T>(GPW_Evaluator::MakeSetOverlap(*g));
}

// The DFT 3-centre tables over the REQUESTED fit basis's grid (doc/GPWPlan §0e).  The fit basis \a c that the
// Hartree/XC term hands us is the one CreateCD/VxcFitBasisSet produced -- a PlaneWave::PlaneWaveFit_IBS, which IS-A
// PlaneWave::PW_Grid_Evaluator carrying the density-fit {G}/grid policy.  Cross-cast to that grid and build the tensor on
// it, so we RETURN THE REQUESTED TABLE rather than overriding the caller's fit-grid choice with the block's own
// (the shared Lattice::Orbital_DFT_IBS mixin dropped \a c).  Bit-identical while the factory wraps DensityGrid();
// the seam is what lets the fit grid diverge (the deferred GGA Vxc densification) without touching these.
template <class T> PairProjection<dcmplx> tGPW_IBS<T>::MakeRepulsion3C(const cFIT_CD_ABS& c) const
{
    const auto& grid = dynamic_cast<const PlaneWave::PW_Grid_Evaluator&>(c);   // throws bad_cast on a non-grid fit basis (loud)
    return GPW_Evaluator::Repulsion3CTensor(std::make_shared<const PlaneWave::PW_Grid_Evaluator>(grid));
}
template <class T> PairProjection<dcmplx> tGPW_IBS<T>::MakeOverlap3C(const cFIT_SF_ABS& c) const
{
    const auto& grid = dynamic_cast<const PlaneWave::PW_Grid_Evaluator&>(c);
    return GPW_Evaluator::Overlap3CTensor(std::make_shared<const PlaneWave::PW_Grid_Evaluator>(grid));
}

// The INSTANCE-scoped 3C accessors (declaration doc in GPW_IBS.C): GPW's matrix-free tensor closures
// capture this evaluator's mutable per-SCF state (CollocMemo D-screen, stream caches), so they must die
// with the run -- the base class's process-wide DBCache route leaked the D-screen into every later
// identical run (the 2026-08-18 first-run anomaly).  Same lookup-or-make contract, instance lifetime.
template <class T> const PairProjection<dcmplx>& tGPW_IBS<T>::Overlap3C(const cFIT_SF_ABS& c) const
{
    auto [it,fresh]=its3Cs.try_emplace("O|"+c.BasisSetID());
    if (fresh) it->second=MakeOverlap3C(c);
    return it->second;
}
template <class T> const PairProjection<dcmplx>& tGPW_IBS<T>::Repulsion3C(const cFIT_CD_ABS& c) const
{
    auto [it,fresh]=its3Cs.try_emplace("R|"+c.BasisSetID());
    if (fresh) it->second=MakeRepulsion3C(c);
    return it->second;
}

// THE TRANSITION PAIR (B2): the fit basis IS-A grid evaluator (the MakeRepulsion3C cast), and its content ID keys
// the evaluator's cached ladder.  Every method checks the bra first -- a pair over two AO sets would collocate
// silently wrong (finding 6).
template <class T> void tGPW_IBS<T>::CheckTransitionPair(const BasisSet::Transition_DFT_IBS& bra, size_t rows, size_t cols) const
{
    if (bra.TransitionBasisID()!=TransitionBasisID())
        throw std::invalid_argument("GPW_IBS: a transition pair over two DIFFERENT AO sets (bra "+bra.TransitionBasisID()
                                    +", ket "+TransitionBasisID()+") -- both must be Bloch sums of ONE molecular block "
                                    "(doc/LinearResponsePlan.md §3d finding 6; a per-k ortho drop does this, pin 22)");
    const size_t n=this->GetNumFunctions();
    if (rows!=n || cols!=n) throw std::invalid_argument("GPW_IBS: the transition δD is not bra x ket");
}
template <class T> ΔGq_Map tGPW_IBS<T>::TransitionRepulsion(const cFIT_CD_ABS& c, const BasisSet::Transition_DFT_IBS& bra,
                                                            const mat_t<dcmplx>& dD, const rvec3_t& q) const
{
    CheckTransitionPair(bra, dD.rows(), dD.columns());
    const auto& grid=dynamic_cast<const PlaneWave::PW_Grid_Evaluator&>(c);   // as MakeRepulsion3C (throws: loud)
    return GPW_Evaluator::TransitionRepulsionField("R|"+c.BasisSetID(), grid, dD, q);
}
template <class T> mat_t<dcmplx> tGPW_IBS<T>::TransitionPotential(const cFIT_CD_ABS& c, const BasisSet::Transition_DFT_IBS& bra,
                                                                  const ΔGq_Map& V) const
{
    CheckTransitionPair(bra, this->GetNumFunctions(), this->GetNumFunctions());
    const auto& grid=dynamic_cast<const PlaneWave::PW_Grid_Evaluator&>(c);
    return GPW_Evaluator::TransitionPotentialField("R|"+c.BasisSetID(), grid, V);
}
template <class T> cvec_t tGPW_IBS<T>::TransitionOnGrid(const cFIT_SF_ABS& c, const BasisSet::Transition_DFT_IBS& bra,
                                                        const mat_t<dcmplx>& dD, const rvec3_t& q) const
{
    CheckTransitionPair(bra, dD.rows(), dD.columns());
    const auto& grid=dynamic_cast<const PlaneWave::PW_Grid_Evaluator&>(c);
    return GPW_Evaluator::TransitionRawField("O|"+c.BasisSetID(), grid, dD, q);
}
template <class T> mat_t<dcmplx> tGPW_IBS<T>::TransitionGridAdjoint(const cFIT_SF_ABS& c, const BasisSet::Transition_DFT_IBS& bra,
                                                                    const cvec_t& v, const rvec3_t& q) const
{
    CheckTransitionPair(bra, this->GetNumFunctions(), this->GetNumFunctions());
    const auto& grid=dynamic_cast<const PlaneWave::PW_Grid_Evaluator&>(c);
    return GPW_Evaluator::TransitionRawAdjoint("O|"+c.BasisSetID(), grid, v, q);
}

template <class T> std::string tGPW_IBS<T>::BasisSetID() const
{
    return Name()+GPW_Evaluator::IDFragment();   // Name + "|mol=..|k=..|cell=..|dEcut=.."
}

template <class T> std::ostream& tGPW_IBS<T>::Write(std::ostream& os) const
{
    return os << Name() << " IBS: " << this->GetNumFunctions() << " periodic Gaussians, " << this->GetSymmetry();
}

template class tGPW_IBS<double>;   // the REAL TRIM block (Step 3; narrowed 1E/PP, complex fit side)
template class tGPW_IBS<dcmplx>;   // the general-k block (== the GPW_IBS alias, unchanged behaviour)

} //namespace
