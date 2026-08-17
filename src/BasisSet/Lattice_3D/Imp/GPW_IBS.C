// File: BasisSet/Lattice_3D/Imp/GPW_IBS.C  GPW_IBS implementation (ctors + identity).
module;
#include <cassert>
#include <iostream>
#include <map>       // MakeSeparablePotentialByL (the per-l KB diagnostic)
#include <memory>
#include <string>
#include <vector>

module qchem.BasisSet.Lattice_3D.GPW_IBS;
import qchem.Symmetry.Factory;              // BlochFactory (the convenience ctor + the k=0 fit-basis irrep)
import qchem.Symmetry.Lattice_3D.BlochQN;   // Symmetry::Lattice_3D::Getk (pry k out of the abstract Bloch irrep)
import qchem.Symmetry.Lattice_3D.SpaceGroup; // DirectOp {W|τ} (the ctor's IBZ raster ops param type)
import qchem.BasisSet.Internal.DB_Cache;    // theCache<dcmplx>() -- process-wide cache for the static PP matrices
                                            // (qcLattice_BS is BasisSet-family, so it may peek at qcBasisSet Internal)
import qchem.BasisSet.Lattice_3D.Evaluators.PW;  // PW_Grid_Evaluator (the fit basis IS-A one; cross-cast target)
import qchem.SymmetrizeMesh;                     // MakeInvariant/FoldMesh (the §6a W1 invariant XC quadrature)

namespace qchem::BasisSet::Lattice_3D
{

GPW_IBS::GPW_IBS(const UnitCell& cell, const sym_t& irrep,
                 std::shared_ptr<const BasisSet::Real_BS> mol, double densityEcut, CellImages images,
                 double cutoffFactor, RasterPolicy raster, double ladderFactor,
                 std::vector<Symmetry::Lattice_3D::DirectOp> directOps, RasterFields rasterFields,
                 std::vector<Symmetry::Lattice_3D::SymOp> magneticOps)
    : BasisSet::IrrepBasisSetImp<dcmplx>(irrep)
    , GPW_Evaluator(std::move(mol), cell, densityEcut, Symmetry::Lattice_3D::Getk(irrep), irrep->IsReal(),
                    images==CellImages::HomeCellOnly, cutoffFactor, raster, ladderFactor,
                    rasterFields) // irrep IS k; IsReal() = the exact TRIM fact (doc/RealComplexPlan.md Step 1)
    , itsMagneticOps(std::move(magneticOps))
{
    SetSymmetryOps(std::move(directOps));   // ONE storage (the evaluator): T1 {G} fold + Vxc-raster star ops
}

// Convenience: build the Bloch irrep from BZ-grid indices and delegate to the primary constructor.
GPW_IBS::GPW_IBS(const UnitCell& cell, const ivec3_t& N, const ivec3_t& kIndex,
                 std::shared_ptr<const BasisSet::Real_BS> mol, double densityEcut, CellImages images,
                 double cutoffFactor, RasterPolicy raster, double ladderFactor, RasterFields rasterFields)
    : GPW_IBS(cell, Symmetry::BlochFactory(N,kIndex), std::move(mol), densityEcut, images, cutoffFactor,
              raster, ladderFactor, {}, rasterFields)
{}

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
BasisSet::cFIT_CD_ABS* GPW_IBS::CreateCDFitBasisSet(const Structure*, const qcMesh::MeshParams&) const
{
    // {G}_rho = DensityGrid() (cutoffFactor*alpha_max, resolving the density product); no relCutoff on the CD grid.
    return new PlaneWaveFit_IBS(GPW_Evaluator::DensityGrid(), Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)),
                                "densityFit");
}
BasisSet::cFIT_SF_ABS* GPW_IBS::CreateVxcFitBasisSet(const Structure*, const qcMesh::MeshParams& mp) const
{
    // The DELTA-fit route builds NO fit basis at all (user 2026-08-01 -- a zero-function pseudo-basis
    // was a null-object smell): its quadrature comes from CreateXCQuadrature below.  This factory
    // serves the PLANE-WAVE fit route.
    // {G}_vxc = relCutoff * {G}_rho.  LDA relCutoff==1 => == DensityGrid(); a GGA's denser grid is not wired yet.
    assert(mp.relCutoff<=1.0 && "GPW: relCutoff>1 (GGA denser Vxc grid) not wired -- the LDA Vxc grid = the CD grid");
    // The Vxc fit basis carries the τ=0 direct ops so its raster STAR-AVERAGES itself under IBZ (empty = no-op).
    return new PlaneWaveFit_IBS(GPW_Evaluator::DensityGrid(), Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)),
                                "xcQuadrature", SymmetryOps());
}

BasisSet::XCQuadrature GPW_IBS::CreateXCQuadrature(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    // The DELTA-fit quadrature (doc/SymmetryUpgradePlan.md §6a).  On a §3-imposed run this basis
    // carries the crystal ops (ctor-injected, like the raster ops PlaneWaveFit_IBS gets):
    //  - Becke grid: the SITE-ADAPTED mesh (W2b) -- rep atoms carry site-group-invariant angular sets,
    //    partners the op-rotated copies -- invariant BY CONSTRUCTION, no group-average growth pass;
    //  - uniform grid (the cross cell): group-average it invariant (W1's MakeInvariant -- a no-op
    //    dedup when the grid is already commensurate).
    // Either way the geometry-fixed orbit fold feeds the engine's per-iteration rho star-average.
    assert(cl && "GPW: the XC quadrature needs the Structure to build its mesh");
    const auto& dops = SymmetryOps();
    if (dops.empty())
        return {std::make_shared<const qcMesh::Mesh>(cl->CreateIntegrationMesh(mp)), {}};

    // On a MAGNETICALLY imposed run (S3) the σ-carrying Shubnikov set drives the quadrature; its
    // SPATIAL parts equal the evaluator's DirectOps (the factory filled both from one group), so the
    // mesh/fold machinery is unchanged -- σ only ADDS the per-op spin actions and the odd-field flags
    // the engine's (ρ,m) channel-pair star-average consumes.
    std::vector<Symmetry::Lattice_3D::SymOp> ops;
    if (!itsMagneticOps.empty()) ops = itsMagneticOps;
    else { ops.reserve(dops.size()); for (const auto& op : dops) ops.push_back({op.W, op.tau}); }
    const double tol = qchem::kMeshMatchTol;   // the ONE mesh-point coincidence tolerance (R2.12)
    const Matrix3D<double>& A = Cell().GetCellMatrix();

    // Ask the CELL for a mesh invariant under the ops.  WHICH strategy realises that (site-adapted by
    // construction for Becke, group-average for uniform) is the cell's business -- this basis only needs
    // the invariance, which is the T2 precondition for the fold below.
    qcMesh::Mesh mesh = Cell().CreateIntegrationMesh(mp, ops);
    auto fold = FoldMesh(mesh, A, ops, tol);
    std::cout << "[Becke grid] imposed symmetry (" << ops.size() << " ops): invariant mesh "
              << mesh.size() << " points in " << fold.Reps()
              << " orbits; rho star-averaged each iteration (doc/SymmetryUpgradePlan.md 6a)" << std::endl;
    if (itsMagneticOps.empty())
        return {std::make_shared<const qcMesh::Mesh>(std::move(mesh)), std::move(fold)};

    // The σ tags (parallel to ops -- the fold's edge opIndex indexes them) and the odd-field zero
    // flags: FRACTIONAL mesh coordinates for the torus-metric fixed-point test.
    std::vector<Symmetry::SpinAction> sigmas;
    sigmas.reserve(ops.size());
    for (const auto& op : ops) sigmas.push_back(op.sigma);
    std::vector<rvec3_t> frac;
    frac.reserve(mesh.size());
    const Matrix3D<double> Ainv = Invert(A);
    for (size_t i=0;i<mesh.size();++i) frac.push_back(Ainv*mesh.Points()[i]);
    auto flags = Symmetry::Lattice_3D::FlipFixedPointsPeriodic(frac, ops, tol);
    size_t nFixed=0; for (char c : flags) if (c) nFixed++;
    std::cout << "[Becke grid] Shubnikov: " << nFixed << " flip-fixed mesh points (m==0 there exactly)"
              << std::endl;
    return {std::make_shared<const qcMesh::Mesh>(std::move(mesh)), std::move(fold),
            std::move(sigmas), std::move(flags)};
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
    // ANALYTIC short assembly (doc/GPWPlan.md 0e-PP step (b), 2026-07-22): exact 3-centre Gaussian lattice
    // sums (LatticeSum1E::MakeLocalGaussian), no grid.  Safe to wire ONLY because step (a) first made the
    // grid LONG standalone-exact (the absolute kappa rule, e^{-kappa/2} pair tails) -- the old grid short's
    // ~0.5 Ha band-limit error used to CANCEL the grid long's, so exact-short + sloppy-long missed the gate.
    // Cross-validated against the kappa-ruled grid short by GPW.LocalPPKappaSelfConverged; falls back to the
    // grid sweep inside MakeLocalPPShort for a model without the closed-Gaussian face.
    return theCache<dcmplx>().Get(IntegralsCache_Base::I2n::LocalPPShort, this, cl->ID(),
        [this,cl,&loc]{ return GPW_Evaluator::MakeLocalPPShort(cl, loc); });
}

hmat_t<dcmplx> GPW_IBS::MakeSeparablePotential(const Structure* cl, const Pseudopotential::SeparablePotential& nl) const
{
    auto* sepR=dynamic_cast<const Pseudopotential::SeparablePotential_R*>(&nl);
    assert(sepR && "GPW MakeSeparablePotential: the KB model must provide the real-space projector face (SeparablePotential_R)");
    return theCache<dcmplx>().Get(IntegralsCache_Base::I2n::SeparablePP, this, cl->ID(),
        [this,cl,sepR]{ return GPW_Evaluator::MakeSeparablePP(cl, *sepR); });
}

std::map<int,hmat_t<dcmplx>> GPW_IBS::MakeSeparablePotentialByL(const Structure* cl, const Pseudopotential::SeparablePotential& nl) const
{
    // Diagnostic face (doc/SphericalLatticePlan.md I0): built on demand, NOT DB-cached -- its consumer
    // (the Ven_PP_NonLocal per-l energy print) caches the result term-side for the run's lifetime.
    auto* sepR=dynamic_cast<const Pseudopotential::SeparablePotential_R*>(&nl);
    assert(sepR && "GPW MakeSeparablePotentialByL: the KB model must provide the real-space projector face (SeparablePotential_R)");
    return GPW_Evaluator::MakeSeparablePPByL(cl, *sepR);
}

// The DFT 3-centre tables over the REQUESTED fit basis's grid (doc/GPWPlan §0e).  The fit basis \a c that the
// Hartree/XC term hands us is the one CreateCD/VxcFitBasisSet produced -- a PlaneWaveFit_IBS, which IS-A
// PW_Grid_Evaluator carrying the density-fit {G}/grid policy.  Cross-cast to that grid and build the tensor on
// it, so we RETURN THE REQUESTED TABLE rather than overriding the caller's fit-grid choice with the block's own
// (the shared EPW_Orbital_DFT_IBS mixin dropped \a c).  Bit-identical while the factory wraps DensityGrid();
// the seam is what lets the fit grid diverge (the deferred GGA Vxc densification) without touching these.
Projector3<dcmplx> GPW_IBS::MakeRepulsion3C(const cFIT_CD_ABS& c) const
{
    const auto& grid = dynamic_cast<const PW_Grid_Evaluator&>(c);   // throws bad_cast on a non-grid fit basis (loud)
    return GPW_Evaluator::Repulsion3CTensor(std::make_shared<const PW_Grid_Evaluator>(grid));
}
Projector3<dcmplx> GPW_IBS::MakeOverlap3C(const cFIT_SF_ABS& c) const
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
