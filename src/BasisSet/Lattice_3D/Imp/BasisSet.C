// File: BasisSet/Lattice_3D/Imp/BasisSet.C  Plane-wave basis-set container + factory implementation.
module;
#include <cassert>
#include <cmath>     // lround (fractional k-point -> integer BZ-grid index); std::fabs (conditioning)
#include <iostream>  // std::cout (the run-start GPW grid diagnostic)
#include <memory>    // std::shared_ptr / std::move (the GPW basis owns the molecular Gaussian basis)
#include <sstream>   // irrep label for the pre-flight basis.perIrrep rows
#include <algorithm> // std::min (min singular value)
#include <vector>    // std::vector (the k-block list + the space-group atom basis)
module qchem.BasisSet.Lattice_3D.BasisSet;
import qchem.BasisSet.Internal.BasisSetImp;   // BasisSetImp<dcmplx> (the generic list-of-IBS container)
import qchem.BasisSet.Lattice_3D.GPW_IBS;     // GPW_IBS (the periodic-Gaussian block GPW_BasisSet owns)
import qchem.Reporting;                        // route the grid diagnostic into the run report when one is open
import qchem.Symmetry.Factory;                // BlochFactory (the Bloch irrep per k)
import qchem.Symmetry;                         // Spin, Irrep (the Bloch block identity for the pre-flight)
import qchem.Symmetry.Lattice_3D.SpaceGroup;   // SpaceGroup::Detect + AtomSite (crystal ops from the cell)
import qchem.Symmetry.Lattice_3D.BZReduction;  // ReduceToIBZ / IBZMesh (fold the MP mesh to the irreducible wedge)
import qchem.UnitCell;                          // UnitCell::GetCellMatrix / ToFractional (the space-group input)
import qchem.Structure;                         // Atom (itsZ / itsR -- the fractional basis for SpaceGroup)
import qchem.LASolver;                         // PivotedCholeskyDrops (the conditioning detector)
import qchem.Blaze;                            // blazem::eigen (grid-free overlap conditioning)
import qchem.Matrix3D;                          // Matrix3D (space-group cell matrix)
import qchem.Types;

namespace qchem::BasisSet::Lattice_3D
{

// ONE plane-wave block per Brillouin-zone k-point: the basis ctor is the single place that enumerates k, so
// the framework's per-irrep loop (each k IS a Bloch irrep) becomes the BZ sum Sum_k w_k.  The KMesh carries
// the points + weights (uniform 1/Nk for an unreduced grid; symmetry-reduced points/weights will plug in
// here later).  N=(1,1,1) -> a single Gamma block.
PW_BasisSet::PW_BasisSet(const ::qchem::Lattice_3D& lat, double Ecut)
{
    ReciprocalLattice recip=lat.Reciprocal();
    const ivec3_t N=lat.GetLimits();
    for (const auto& kp : lat.MakeKMesh())
    {
        ivec3_t ik(std::lround(kp.k.x*N.x), std::lround(kp.k.y*N.y), std::lround(kp.k.z*N.z));
        auto* pw=new PlaneWave_IBS(recip, Symmetry::BlochFactory(N, ik, kp.weight), Ecut);
        Insert(pw);                                 // BasisSetImp takes ownership (no PP: the Vpseudo
                                                    // Hamiltonian term owns the pseudopotential model)
    }
}

Complex_BS* Factory(Type type, const ::qchem::Lattice_3D& lat, double Ecut)
{
    assert(type==Type::PW);
    return new PW_BasisSet(lat, Ecut);
}

// GPW: one Bloch block of periodic Gaussians per Brillouin-zone k-point, built from the molecular basis over
// the cell's atoms -- the exact mirror of PW_BasisSet above (the ctor is the sole place that enumerates k, so
// the framework's per-irrep loop becomes the BZ sum Sum_k w_k).  N=(1,1,1) -> a single Gamma block.  The
// molecular basis is shared (shared_ptr) across the k-blocks.  Lattice images are eps-converged series
// summed inside the molecular seam (no radius exists); homeCellOnly is the finite "molecule in a box" MODE
// (image-free by definition -- every k-block identical), used by the finite==lattice gates.
// The k-mesh as (integer grid index, BZ weight) pairs -- either the FULL Monkhorst-Pack grid or, when
// p.reduceBZ, its irreducible wedge folded under the crystal point group (τ=0 symmorphic ops only -- the safe
// W-only guard; doc/GPWPlan1.md item 3).  IBZ needs the density symmetrized over each star to be exact; the
// caller is responsible for that (this only reduces which blocks are BUILT + carries the star weights).
struct KBlock { ivec3_t ik; double weight; };
// The τ=0 reciprocal point ops of the crystal (the W-only guard), or {} when we are NOT folding (reduceBZ off)
// -- the fold and the density star-average are one package, so a non-folded run carries no ops (trivial {E}).
static std::vector<Matrix3D<double>> ReciprocalOps(const ::qchem::Lattice_3D& lat, const GPWParams& p)
{
    namespace SL = qchem::Symmetry::Lattice_3D;
    if (!p.reduceBZ) return {};
    const UnitCell&  cell = lat.GetUnitCell();
    const Structure& st   = cell;   // iterate atoms via the PUBLIC Structure interface (GetAtom is protected)
    std::vector<SL::AtomSite> basis;
    for (Atom* a : st)
        basis.push_back({a->itsZ, cell.ToFractional(a->itsR)});
    return SL::SpaceGroup::Detect(cell.GetCellMatrix(), basis).ReciprocalPointOps(/*timeReversal*/true,
                                                                                  /*symmorphicOnly*/true);
}

static std::vector<KBlock> BuildKBlocks(const ::qchem::Lattice_3D& lat, const GPWParams& p,
                                        const std::vector<Matrix3D<double>>& ops)
{
    namespace SL = qchem::Symmetry::Lattice_3D;
    const rvec3_t kShift=p.kShift;
    const ivec3_t N=lat.GetLimits();
    std::vector<KBlock> blocks;
    if (!p.reduceBZ)
    {
        for (const auto& kp : lat.MakeKMesh(kShift))
            blocks.push_back({ivec3_t(std::lround(kp.k.x*N.x-kShift.x), std::lround(kp.k.y*N.y-kShift.y),
                                      std::lround(kp.k.z*N.z-kShift.z)), kp.weight});
        return blocks;
    }
    // IBZ: fold the MP grid under the (pre-computed) τ=0 reciprocal ops.  Every k-resolved BAND quantity obeys
    // f(Uk)=f(k) so the star weight is exact for the eigenvalue/occupation sum; the DENSITY needs the star
    // symmetrization (the composite ctor-injects the SAME ops, so a folded run is symmetrized automatically).
    SL::IBZMesh ibz = SL::ReduceToIBZ(N, kShift, ops);
    std::cout << "[IBZ] " << ibz.FullSize() << " k-points -> " << ibz.points.size()
              << " irreducible (point group |ops|=" << ops.size() << "); Σw=" << ibz.WeightSum() << std::endl;
    for (const auto& q : ibz.points) blocks.push_back({q.index, q.weight});
    return blocks;
}

GPW_BasisSet::GPW_BasisSet(const ::qchem::Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                           const GPWParams& p)
    : itsReciprocalOps(ReciprocalOps(lat, p))   // computed ONCE; exposed via GetReciprocalPointOps for the density
{
    const rvec3_t kShift=p.kShift;
    const ivec3_t N=lat.GetLimits();
    const GPW_IBS* first=nullptr;
    for (const auto& kb : BuildKBlocks(lat, p, itsReciprocalOps))
    {
        // Build the Bloch irrep WITH its BZ weight (star weight under IBZ) and the primary sym_t ctor -- the
        // weight carries the Sum_k w_k so the BZ-summed charge/energy are per-cell, not xNk.
        auto* b=new GPW_IBS(lat.GetUnitCell(), Symmetry::BlochFactory(N, kb.ik, kb.weight, kShift),
                            mol, p.densityEcut, p.images, p.cutoffFactor, p.raster, p.ladderFactor);   // mol shared across k-blocks
        if (!first) first=b;
        Insert(b);
    }
    // GRID DIAGNOSTIC (doc/GPWPlan §0e): every stored grid, once per run.  When a run report is open the
    // orchestrator now emits grids EXPLICITLY (EmitGpwGrids) AFTER its conditioning pre-flight, so a singular
    // basis aborts before the grid ladder is ever built (fail-fast -- the report surfaced that building grids
    // here forced EnsureLevels before the eigen-analysis could veto the basis).  So only the legacy unbracketed
    // console path stays here; ReportGrids also forces the grid build, which is exactly what we must NOT do early.
    if (first && report::Depth() == 0) first->ReportGrids(std::cout);
}

// Compact irrep label for the report (symmetry symbol + spin arrow), via the Streamable Write().
static std::string GpwIrrepLabel(const Irrep& q) { std::ostringstream os; q.Write(os); return os.str(); }

void EmitGpwGrids(const Complex_BS& bs)
{
    // The grids are k-independent, so the first GPW block speaks for all (and forces EnsureLevels HERE --
    // after the pre-flight, so a vetoed basis never reaches this).
    for (auto* b : const_cast<Complex_BS&>(bs).Iterate<GPW_IBS>()) { b->EmitGridsReport(); return; }
}

size_t VetGpwConditioning(const Complex_BS& bs)
{
    namespace rpt = qchem::report;
    rpt::Set("nFunctions", (long)bs.GetNumFunctions());   // basis-level scalar (cursor at the open "basis" section)
    size_t total = 0;
    for (auto* b : const_cast<Complex_BS&>(bs).Iterate<GPW_IBS>())
    {
        const hmat_t<dcmplx>& S = b->Overlap();           // ANALYTIC Bloch overlap -- no grids needed
        const std::string label = GpwIrrepLabel(b->GetIrrep(Spin::None));
        {
            rpt::Row row("perIrrep");
            rpt::Set("irrep",      label);
            rpt::Set("nFunctions", (long)b->GetNumFunctions());
            rvec_t d; mat_t<dcmplx> U; blazem::eigen(S, d, U);       // ascending eigenvalues of the Hermitian S
            const double mn = d[0], mx = d[d.size()-1];
            double msv = std::fabs(d[0]); for (double v : d) msv = std::min(msv, std::fabs(v));
            rpt::Set("lambdaMin", mn);
            rpt::Set("lambdaMax", mx);
            rpt::Set("cond", msv > 0 ? mx/msv : 0.0);
        }
        for (size_t idx : qchem::PivotedCholeskyDrops<dcmplx>(S))    // redundant AOs (empty on a healthy basis)
        {
            rpt::Row r("removed");
            rpt::Set("irrep", label);
            rpt::Set("index", (long)idx);
            ++total;
        }
    }
    return total;
}

Complex_BS* GPWFactory(const ::qchem::Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                       const GPWParams& p)
{
    return new GPW_BasisSet(lat, std::move(mol), p);
}
Complex_BS* GPWFactory(const ::qchem::Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                       double densityEcut, rvec3_t kShift, CellImages images, double cutoffFactor)
{
    return new GPW_BasisSet(lat, std::move(mol),
                            GPWParams{densityEcut, cutoffFactor, RasterPolicy::BallOnly, images, kShift});
}

} //namespace
