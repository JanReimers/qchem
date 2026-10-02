// File: BasisSet/Lattice/Imp/BasisSet.C  Plane-wave basis-set container + factory implementation.
module;
#include <cassert>
#include <functional> // VetStageTrim's basis maker
#include <stdexcept>  // VetStageTrim's refusals
#include <cmath>     // lround (fractional k-point -> integer BZ-grid index); std::fabs (conditioning)
#include <iostream>  // std::cout (the run-start GPW grid diagnostic)
#include <memory>    // std::shared_ptr / std::move (the GPW basis owns the molecular Gaussian basis)
#include <sstream>   // irrep label for the pre-flight basis.perIrrep rows
#include <variant>   // the mixed child-slot visits (Step 3c-2)
#include <algorithm> // std::min (min singular value)
#include <cstdlib>   // std::getenv/std::atoi
#include <vector>    // std::vector (the k-block list + the space-group atom basis)
module qchem.BasisSet.Lattice.BasisSet;
import qchem.Symmetry.Lattice_3D.Fold;   // RequireClosedGroup -- the imposition gate
import qchem.RunPolicy;   // theRunPolicy().StreamFold() -- the T3.2 fold, declared with the deviations (N5)
import qchem.BasisSet.Internal.BasisSetImp;   // BasisSetImp<dcmplx> (the generic list-of-IBS container)
import qchem.BasisSet.Gaussian.Lattice.GPW_IBS;     // Gaussian::GPW_IBS (the periodic-Gaussian block GPW_BasisSet owns)
import qchem.BasisSet.Gaussian.Lattice.GPW_Evaluator; // Gaussian::GPW_Evaluator (the shared grid face the mixed EmitGpwGrids visit casts to)
import qchem.BasisSet.Orbital_1E_IBS;         // Real_OIBS (the molecular orbital block -- the stream-fold cross-cast)
import qchem.BasisSet.Gaussian.Lattice.LatticeSum1E;  // SetStreamSymmetryOps (the T3 route (b) stream fold, §6b)
import qchem.Reporting;                        // route the grid diagnostic into the run report when one is open
import qchem.Symmetry.Factory;                // BlochFactory (the Bloch irrep per k)
import qchem.Symmetry;                         // Spin, Irrep (the Bloch block identity for the pre-flight)
import qchem.Symmetry.Lattice_3D.SpaceGroup;   // SpaceGroup::Detect + AtomSite (crystal ops from the cell)
import qchem.Symmetry.Lattice_3D.BZReduction;  // ReduceToIBZ / IBZMesh (fold the MP mesh to the irreducible wedge)
import qchem.UnitCell;                          // UnitCell::GetCellMatrix / ToFractional (the space-group input)
import qchem.Structure;                         // Atom (itsZ / itsR -- the fractional basis for SpaceGroup)
import qchem.LASolver;                         // PivotedCholeskyDrops (the conditioning detector)
import qchem.BasisSet.AoShellSource;           // the shells a Bloch block is built from (VetStageTrim)
import qchem.Blaze;                            // blazem::eigen (grid-free overlap conditioning)
import qchem.Matrix3D;                          // Matrix3D (space-group cell matrix)
import qchem.Types;

namespace qchem::BasisSet::Lattice
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
        auto* pw=new PlaneWave::PlaneWave_IBS(recip, Symmetry::BlochFactory(N, ik, kp.weight), Ecut);
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
// p.imposeSymmetry, its irreducible wedge folded under the crystal point group (τ=0 symmorphic ops only -- the safe
// W-only guard; doc/GPWPlan1.md item 3).  IBZ needs the density symmetrized over each star to be exact; the
// caller is responsible for that (this only reduces which blocks are BUILT + carries the star weights).
struct KBlock { ivec3_t ik; double weight; };
// The crystal point ops from ONE space-group detection (or {} when NOT folding -- fold + density star-average are
// one package, so a non-folded run carries no ops = trivial {E}).  THREE faces, two distinct roles:
//  - recipFold: reciprocal LINEAR ops U (+time reversal, τ=0 W-only guard) -- folds the k-mesh.  TR supplies the
//    inversion a non-symmorphic Td subgroup lacks, so the fold is already maximal (Oh) without the glide phase.
//  - recipDensity {U|τ}: the FULL crystal point group with fractional translations -- star-averages ρ̃ in G-space
//    with the glide phase e^{+2πi(Um)·τ} (SymmetrizeGMap).  NO time reversal: a real density has the crystal's
//    own point symmetry, not an imposed inversion.
//  - directDensity {W|τ}: the real-space partner (raster star-average, voxel shift N·τ).
// For a SYMMORPHIC crystal every τ=0, so recipDensity/directDensity collapse to the old W-only ops (no phase, no
// shift) and this is bit-compatible with the symmorphic IBZ path.  The reciprocal↔direct + time-reversal + glide
// knowledge lives in SpaceGroup, not here (doc/GPWPlan1.md items 3 + 5).
struct CrystalPointOps
{
    Symmetry::Lattice_3D::SymmetryPolicy                 policy;         //!< the ONE run-level imposition decision (§3)
    std::vector<Matrix3D<double>>                        recipFold;      //!< linear reciprocal ops (+TR) for the k-fold ({} unless imposed)
    std::vector<Symmetry::Lattice_3D::ReciprocalOp>      recipDensity;   //!< {U|τ} for the G-space density star-average ({} unless imposed)
    std::vector<Symmetry::Lattice_3D::DirectOp>          directDensity;  //!< {W|τ} for the real-space raster star-average ({} unless imposed)
    std::vector<Symmetry::Lattice_3D::ReciprocalOp>      recipDetected;  //!< the FULL detected {U|τ} set, ALWAYS (§3 diagnostic reference)
    //! The SHUBNIKOV group when the run imposes a MAGNETIC decoration (S3, {} otherwise): direct-frame
    //! \f$\{W|\tau,\sigma\}\f$ ops.  The legacy faces above are then its SPATIAL projections (σ dropped)
    //! -- COSET-COMPLETE, i.e. they include the anti-translation of a magnetically doubled cell, which
    //! the detected one-coset-per-W group misses -- and serve every EVEN (total-density) consumer
    //! unchanged; σ itself is consumed only by the channel-pair (ρ,m) machinery in the XC quadrature.
    std::vector<Symmetry::Lattice_3D::SymOp>             magneticDirect;
};
//! ONE place resolves the run's symmetry policy (doc/SymmetryUpgradePlan.md §3): detection ALWAYS runs
//! (cheap, and a free run needs the full group as the order-parameter diagnostic's reference), but the
//! IMPOSED faces (k-fold + density star-average -- one package) fill only when the caller asserted
//! imposition (p.imposeSymmetry).  Downstream consumers read the policy from here, never p.imposeSymmetry again.
//!
//! S3 (Shubnikov): imposition on a run carrying a MAGNETIC decoration (p.siteSpins with a genuine
//! staggering) resolves to the magnetic group -- imposing the grey group on a staggered run would
//! star-average +m onto -m sites and erase the order (the reason MnO ran free all campaign).  Grey
//! imposition of such a cell remains reachable only through the §3 release-audit machinery.
static CrystalPointOps DetectPointOps(const ::qchem::Lattice_3D& lat, const GPWParams& p)
{
    namespace SL = qchem::Symmetry::Lattice_3D;
    // Detection + the UnitCell -> {A, sites} adapter now live ON the lattice (Lattice_3D::GetSpaceGroup,
    // 2026-08-01) -- common to every lattice basis flavour; this factory only resolves the §3 POLICY.
    const SL::SpaceGroup& sg = lat.GetSpaceGroup();

    CrystalPointOps ops;
    ops.policy.impose = p.imposeSymmetry ? SL::SymmetryPolicy::Impose::FullGroup
                                   : SL::SymmetryPolicy::Impose::None;
    ops.recipDetected = sg.ReciprocalOps();                                        // full {U|τ}: the diagnostic reference
    if (!ops.policy.Imposes()) return ops;

    const bool magnetic = [&]{ for (int s : p.siteSpins) if (s!=0) return true; return false; }();
    if (magnetic)
    {
        ops.magneticDirect = lat.ShubnikovOps(p.siteSpins);
        size_t nFlip=0; for (const auto& op : ops.magneticDirect) if (op.sigma==SL::SpinAction::Flip) nFlip++;
        std::cout << "[symmetry] IMPOSING the SHUBNIKOV group of the declared ordering: "
                  << ops.magneticDirect.size() << " ops (" << ops.magneticDirect.size()-nFlip
                  << " spatial + " << nFlip << " spin-flip; the grey group would erase the order)" << std::endl;
        // The LEGACY faces get ONLY the σ=None (sublattice-preserving) subgroup.  MEASURED (MnO run 37,
        // 2026-08-11): they must NOT get the flip ops' spatial parts -- the composite star-average is
        // applied PER CHANNEL (each channel's GetFourierDensity/GetRhoOnGrid symmetrizes ITSELF, and the
        // ρ̃ MIXER consumes the channels separately through exactly those faces), and a sublattice-swap op
        // is not a symmetry of ONE channel: averaging under it equalizes the channels and annihilates m
        // MACHINE-EXACTLY at iteration 1 (the seed reads m1=+0.717, the first mixed density m1=6e-14).
        // Each channel IS invariant under the σ=None subgroup, so that is the correct per-channel
        // projector; the flip content is a PAIR property and is enforced where the pair is in hand -- the
        // XC engine's (ρ,m) star-average under the full σ-carrying set (itsMagneticOps).  The total's
        // anti-translation component is then driven by the SCF rather than projected -- a PARTIAL
        // projector is legal (folding merely reduced, never wrong).
        // ★ THE IMPOSITION GATE (doc/SymmetryUpgradePlan.md "SUPERCELLS", 2026-09-07): a star-average is a
        // PROJECTOR only if its op set is a GROUP, and nothing checked that until a Si 2x2x2 supercell
        // diverged from its first Fock on a set with 864 of 2304 products outside it.  Checked HERE
        // because this is the one place a run decides what it imposes.  Closure, NOT completeness -- the
        // sigma=None subgroup extracted just below is a legal partial projector (see the note there).
        SL::RequireClosedGroup(ops.magneticDirect, "GPW imposition (Shubnikov group)");
        for (const auto& op : ops.magneticDirect)
            if (op.sigma==SL::SpinAction::None)
            {
                ops.directDensity.push_back({op.W, op.tau});
                ops.recipDensity .push_back({Transpose(op.W), op.tau});   // the G-index scatter U = W^T
            }
        // k-fold: Γ-only for now -- the magnetic little-group bookkeeping at k≠Γ is the T3.4b-adjacent
        // increment.  Leaving recipFold empty means a multi-k imposed magnetic run keeps its FULL mesh
        // (folding is merely absent, never wrong).
        return ops;
    }
    ops.recipFold     = sg.ReciprocalPointOps(/*timeReversal*/true, /*symmorphicOnly*/true);  // fold: linear + TR
    ops.recipDensity  = sg.ReciprocalOps();                                    // density G-space: {U|τ}
    ops.directDensity = sg.DirectOps();                                        // density raster: {W|τ}
    // The k-sum only respects the mesh-symmetry SUBGROUP (a shifted fcc MP mesh is not Oh-invariant), and the
    // fold above already uses just that subgroup.  Star-averaging rho over ops the mesh does not respect would
    // sample a LARGER k-set than the one whose bands are summed -- so impose the same subgroup on rho.
    {
        const ivec3_t N=lat.GetLimits();
        auto bad=[&](const Matrix3D<double>& U){ return !SL::MapsMeshOntoItself(U,N,p.kShift); };
        std::erase_if(ops.recipDensity,  [&](const auto& o){ return bad(o.U); });
        std::erase_if(ops.directDensity, [&](const auto& o){ return bad(Transpose(o.W)); });
    }
    // ★ THE IMPOSITION GATE -- see the note in the magnetic branch above.  The reciprocal set is the
    // transpose of this one, so closure of {W|τ} settles both.
    // DirectOp carries no σ (a spatial group), so it lifts to SymOp with σ=None -- the closure test is
    // then exactly the spatial one.
    {
        std::vector<SL::SymOp> spatial;
        spatial.reserve(ops.directDensity.size());
        for (const auto& op : ops.directDensity) spatial.push_back({op.W, op.tau});
        SL::RequireClosedGroup(spatial, "GPW imposition (spatial space group)");
    }
    return ops;
}

static std::vector<KBlock> BuildKBlocks(const ::qchem::Lattice_3D& lat, const GPWParams& p,
                                        const CrystalPointOps& cops)
{
    namespace SL = qchem::Symmetry::Lattice_3D;
    const std::vector<Matrix3D<double>>& ops = cops.recipFold;
    const rvec3_t kShift=p.kShift;
    const ivec3_t N=lat.GetLimits();
    std::vector<KBlock> blocks;
    if (!cops.policy.Imposes())
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
{
    const rvec3_t kShift=p.kShift;
    const ivec3_t N=lat.GetLimits();
    const Gaussian::GPW_Evaluator* first=nullptr;   // the shared evaluator base -- serves ReportGrids for EITHER block scalar
    CrystalPointOps ops = DetectPointOps(lat, p);   // ONE detection + the §3 policy: fold ops + {U|τ} (ρ̃) + {W|τ} (raster)
    itsReciprocalOps = ops.recipDensity;            // exposed via GetReciprocalPointOps for the composite G-space density
    itsDetectedOps   = ops.recipDetected;           // the FULL detected group, imposed or not (§3 diagnostic reference)
    const std::vector<KBlock> blocks = BuildKBlocks(lat, p, ops);
    // T3 route (b) STREAM FOLD (doc/SymmetryUpgradePlan.md §6b/T3.2): on an IMPOSED Γ-ONLY run, fold the
    // shared molecular evaluator's (pair, R) collocation streams under the crystal ops -- reduced build +
    // replay (demand and per-iteration scatter/gather ÷ ~orbit factor); the existing per-iteration
    // SymmetrizeGMap/SymmetrizeRaster sites complete the group-average, so no raster commensurability is
    // needed.  Γ-only because a general-k block needs the little-group restriction + edge phases (T3.4);
    // multi-k IBZ runs keep full streams.
    // ARMED BY DEFAULT since 2026-08-19 (T3.5; OpenWork Step 2).  It was opt-in from 2026-08-03 because the
    // fold then imposed STRICTLY MORE than imposeSymmetry itself: the ρ star-average tolerates a
    // symmetry-BROKEN iterate D (it projects it pointwise), while the reduced replay read only
    // orbit-REPRESENTATIVE D elements -- i.e. it SAMPLED the orbit and thereby asserted D symmetric, which a
    // DEGENERATE OPEN SHELL at integer fill breaks permanently (measured: the imposed Si pseudo-atom-in-a-box
    // p² run flipped from the benign rotating-ρ mode into charge-transfer sloshing, ~0.26 Ha off).  The
    // replay now reads the ORBIT-PROJECTED D instead (NR_Evaluator::FoldProjectedD), and since collocation is
    // linear and equivariant that makes P ρ_red[D] == ρ[P D] == P ρ_full[D] for ANY iterate: the folded and
    // the unfolded imposed run compute the SAME density, so arming is a pure COST decision and needs no
    // auto-arm criterion.  GPW_STREAM_FOLD=0 is the opt-OUT, resolved through qchem::theRunPolicy() with the
    // other CP2K deviations (N5) rather than read here -- a gate that A/Bs it in one process calls
    // ReresolveRunPolicy() between the arms.
    if (ops.policy.Imposes() && blocks.size()==1
        && blocks[0].ik.x==0 && blocks[0].ik.y==0 && blocks[0].ik.z==0
        && kShift.x==0.0 && kShift.y==0.0 && kShift.z==0.0)
    {
        if (theRunPolicy().StreamFold())
        {
            // S3 guard: the stream fold asserts D ITSELF symmetric under its ops PER CHANNEL, and a
            // flip op relates D_up to D_dn -- so a magnetic imposition may fold streams only under the
            // σ=None (sublattice-preserving) subgroup.  Folding is merely reduced, never wrong.
            namespace SL = qchem::Symmetry::Lattice_3D;
            std::vector<SL::DirectOp> streamOps;
            if (ops.magneticDirect.empty()) streamOps = ops.directDensity;
            else for (const auto& op : ops.magneticDirect)
                if (op.sigma==SL::SpinAction::None) streamOps.push_back({op.W, op.tau});
            const BasisSet::Real_OIBS* orb=nullptr;
            for (auto ibs : const_cast<BasisSet::Real_BS&>(*mol).Iterate<BasisSet::Real_OIBS>()) { orb=ibs; break; }
            // THE NARROW FACE (ISP split 2026-09-08): arming the fold is all this does, so it names the
            // two-method StreamFoldable capability rather than the whole periodic aggregate.
            if (const auto* sf=dynamic_cast<const Gaussian::StreamFoldable*>(orb))
            {
                const size_t used=sf->SetStreamSymmetryOps(streamOps, lat.GetUnitCell());
                std::cout << "[stream fold] imposed Γ run: " << used << "/" << streamOps.size()
                          << " ops folded into the collocation streams" << std::endl;
            }
        }
    }
    for (const auto& kb : blocks)                                // the k-fold uses the linear reciprocal ops (+TR)
    {
        // Build the Bloch irrep WITH its BZ weight (star weight under IBZ) and the primary sym_t ctor -- the
        // weight carries the Sum_k w_k so the BZ-summed charge/energy are per-cell, not xNk.
        const sym_t irrep = Symmetry::BlochFactory(N, kb.ik, kb.weight, kShift);
        // THE FACTORY DECISION (doc/RealComplexPlan.md §1, Step 3c-3): a block is real ⇔ its irrep is
        // (TRIM, exact integer arithmetic) ∧ every Hamiltonian term preserves realness (the composition
        // root's fact, threaded in as GPWParams::hamPreservesReal).  One scalar-generic build serves
        // both alternatives; the typed Insert (Step 3b) files each under its own child slot.
        auto build=[&]<class U>()
        {
            auto* b=new Gaussian::tGPW_IBS<U>(lat.GetUnitCell(), irrep,
                                    mol, p.densityEcut, p.images, p.cutoffFactor, p.raster, p.ladderFactor,
                                    ops.directDensity,    // mol shared across k-blocks; {W|τ} = the IBZ raster star ops
                                    p.rasterFields,       // field-sharpness routing (HartreeOnly = the Becke-XC partner)
                                    ops.magneticDirect);  // Shubnikov {W|τ,σ} on a magnetic imposition (S3; {} = grey)
            if (!first) first=b;
            Insert(b);
        };
        if (p.hamPreservesReal && irrep->IsReal()) build.template operator()<double>();
        else                                       build.template operator()<dcmplx>();
    }
    // GRID DIAGNOSTIC (doc/GPWPlan §0e): every stored grid, once per run.  When a run report is open the
    // orchestrator now emits grids EXPLICITLY (EmitGpwGrids) AFTER its conditioning pre-flight, so a singular
    // basis aborts before the grid ladder is ever built (fail-fast -- the report surfaced that building grids
    // here forced EnsureLevels before the eigen-analysis could veto the basis).  So only the legacy unbracketed
    // console path stays here; ReportGrids also forces the grid build, which is exactly what we must NOT do early.
    if (first && report::Depth() == 0) first->ReportGrids(std::cout);
}

// Compact irrep label for the report (symmetry symbol + spin arrow), via the Streamable Write().


void EmitGpwGrids(const Complex_BS& bs)
{
    // The grids are k-independent, so the first GPW block speaks for all (and forces EnsureLevels HERE --
    // after the pre-flight, so a vetoed basis never reaches this).  MIXED-AWARE (Step 3c-2): walk the
    // typed child slot -- the evaluator face is shared by both block scalars, so one generic visit serves.
    auto* imp=dynamic_cast<const BasisSet::BasisSetImp<dcmplx>*>(&bs);
    assert(imp && "EmitGpwGrids: the periodic basis set must be a BasisSetImp");
    for (size_t i=0;i<bs.GetNumIBS();++i)
    {
        const bool done=std::visit([&](const auto& b)
        {
            if (auto* ev=dynamic_cast<const Gaussian::GPW_Evaluator*>(b.get())) { ev->EmitGridsReport(); return true; }
            return false;
        }, imp->GetChild(i));
        if (done) return;
    }
}

size_t VetGpwConditioning(const Complex_BS& bs)
{
    namespace rpt = qchem::report;
    rpt::Set("nFunctions", (long)bs.GetNumFunctions());   // basis-level scalar (cursor at the open "basis" section)
    size_t total = 0;
    // MIXED-AWARE (Step 3c-2): one generic visit over the typed child slot serves both block scalars --
    // the eigen-analysis and the Cholesky drop detector are scalar-generic already.
    auto* imp=dynamic_cast<const BasisSet::BasisSetImp<dcmplx>*>(&bs);
    assert(imp && "VetGpwConditioning: the periodic basis set must be a BasisSetImp");
    for (size_t i=0;i<bs.GetNumIBS();++i)
        std::visit([&](const auto& b)
        {
            const auto& S = b->Overlap();                 // ANALYTIC Bloch overlap -- no grids needed
            using U=typename std::decay_t<decltype(S)>::ElementType;
            const std::string label = b->GetSymmetry().GetLabel();   // the block's SPATIAL label (k-point)
            {
                rpt::Row row("perIrrep");
                rpt::Set("irrep",      label);
                rpt::Set("nFunctions", (long)b->GetNumFunctions());
                rpt::Set("real",       b->IsReal());   // TRIM fact from the irrep (doc/RealComplexPlan.md Step 1)
                rpt::Set("runsReal",   std::is_same_v<U,double>);   // the block's CONSTRUCTED scalar (3c-3 evidence)
                rvec_t d; mat_t<U> Uv; blazem::eigen(S, d, Uv);      // ascending eigenvalues of the Hermitian S
                const double mn = d[0], mx = d[d.size()-1];
                double msv = std::fabs(d[0]); for (double v : d) msv = std::min(msv, std::fabs(v));
                rpt::Set("lambdaMin", mn);
                rpt::Set("lambdaMax", mx);
                rpt::Set("cond", msv > 0 ? mx/msv : 0.0);
            }
            for (size_t idx : qchem::PivotedCholeskyDrops<U>(S))     // redundant AOs (empty on a healthy basis)
            {
                rpt::Row r("removed");
                rpt::Set("irrep", label);
                rpt::Set("index", (long)idx);
                ++total;
            }
        }, imp->GetChild(i));
    return total;
}

void PrintGpwConditioning(const Complex_BS& bs)
{
    auto* imp=dynamic_cast<const BasisSet::BasisSetImp<dcmplx>*>(&bs);
    assert(imp && "PrintGpwConditioning: the periodic basis set must be a BasisSetImp");
    double lminWorst=1e300, condWorst=0.0, lmaxBest=0.0;
    std::string lminAt, condAt;
    for (size_t i=0;i<bs.GetNumIBS();++i)
        std::visit([&](const auto& b)
        {
            const auto& S = b->Overlap();
            using U=typename std::decay_t<decltype(S)>::ElementType;
            rvec_t d; mat_t<U> Uv; blazem::eigen(S, d, Uv);
            const double mn=d[0], mx=d[d.size()-1];
            double msv=std::fabs(d[0]); for (double v : d) msv=std::min(msv, std::fabs(v));
            const double cond=msv>0 ? mx/msv : 1e300;
            const std::string label=b->GetSymmetry().GetLabel();
            if (mn<lminWorst)    { lminWorst=mn;    lminAt=label; }
            if (cond>condWorst)  { condWorst=cond;  condAt=label; }
            lmaxBest=std::max(lmaxBest, mx);
        }, imp->GetChild(i));
    std::cout<<"[basis cond] "<<bs.GetNumIBS()<<" k-block(s): min eig S = "<<lminWorst<<" ("<<lminAt
             <<"), max cond(S) = "<<condWorst<<" ("<<condAt<<"), max eig S = "<<lmaxBest<<std::endl;
}

// ------------------------------------------------------------------ the VET-STAGE TRIM (pin 22)

VetTrimResult VetStageTrim(const ::qchem::Lattice_3D& lat,
                           const std::function<std::shared_ptr<const BasisSet::Real_BS>(const Gaussian::ShellTrim&)>& make,
                           const GPWParams& p, double orthoTol)
{
    // The cell's sites, to name a shell's element from its centre.
    std::vector<std::pair<int,rvec3_t>> sites;
    lat.GetUnitCell().ForEachSite([&](int Z, const rvec3_t& R, bool){ sites.push_back({Z,R}); });
    auto elementAt=[&](const rvec3_t& c)
    {
        for (const auto& [Z,R] : sites)
        {
            const rvec3_t d=c-R;
            if (d.x*d.x+d.y*d.y+d.z*d.z < 1e-12) return Z;
        }
        throw std::logic_error("VetStageTrim: a basis shell whose centre is no site of the cell");
    };
    GPWParams full=p;
    full.imposeSymmetry=false;                 // the rank decision is the BASIS's: every k of the mesh, no fold
    full.siteSpins.clear();
    full.densityEcut=0.0;                      // the vet needs ANALYTIC overlaps only: no DFT tier, no grid ladder
    VetTrimResult r;
    // The smallest overlap eigenvalue over the k-mesh, and where -- the conditioning number a reader judges a
    // basis by (user 2026-09-28: the trim must show it BEFORE and AFTER).  One per pass; [0] = the untrimmed basis.
    struct LambdaMin { double value=1e300; std::string k; };
    std::vector<LambdaMin> lambdas;
    auto fmtLambda=[](const LambdaMin& l){ std::ostringstream os; os << l.value << " at k=" << l.k; return os.str(); };
    for (int pass=0; pass<64; ++pass)
    {
        r.mol=make(r.trim);
        GPW_BasisSet bs(lat, r.mol, full);
        const Complex_BS& cbs=bs;                  // the public face (GetNumIBS is protected on the Imp)
        // Per pass: every implicated (element, l), and the most diffuse shell of each (its exponents).
        struct Candidate { int Z; int l; rvec_t exponents; double alphaMin; };
        std::vector<Candidate> hits;
        size_t nDrops=0;
        LambdaMin lmin;
        for (size_t i=0;i<cbs.GetNumIBS();++i)
            std::visit([&](const auto& b)
            {
                const auto& S=b->Overlap();                                  // ANALYTIC, grid-free
                using U=typename std::decay_t<decltype(S)>::ElementType;
                {
                    rvec_t w; mat_t<U> Uw; blazem::eigen(S, w, Uw);          // ascending
                    if (w[0]<lmin.value) {lmin.value=w[0]; std::ostringstream os; os << b->GetSymmetry(); lmin.k=os.str();}
                }
                const std::vector<size_t> drops=qchem::PivotedCholeskyDrops<U>(S, orthoTol);
                if (drops.empty()) return;
                nDrops+=drops.size();
                const auto* src=dynamic_cast<const BasisSet::AoShellSource*>(b.get());   // abstract -> abstract
                if (!src) throw std::logic_error("VetStageTrim: a Bloch block that cannot name its shells (no AoShellSource)");
                const auto shells=src->GetAoShells();
                for (size_t idx : drops)
                    for (const auto& sh : shells)
                    {
                        if (idx<sh.offset || idx>=sh.offset+sh.nComponents()) continue;
                        const int l=sh.rep->L();
                        if (l>=2 && sh.nComponents()!=size_t(2*l+1))
                            throw std::runtime_error("VetStageTrim: a CARTESIAN l>=2 shell is implicated -- its l-2 "
                                "contaminants make 'the most diffuse shell of that l' the wrong cut (doc/Pins.md pin 22 (a)); "
                                "run the spherical view (GPW_SPHERICAL) for a vet-stage trim");
                        const int Z=elementAt(sh.center);
                        // The most diffuse shell of this (element, l) among ALL the block's shells.
                        Candidate best{Z,l,{},1e300};
                        for (const auto& t : shells)
                        {
                            if (t.rep->L()!=l || elementAt(t.center)!=Z || t.exponents.size()==0) continue;
                            double amin=t.exponents[0]; for (double a : t.exponents) amin=std::min(amin,a);
                            if (amin<best.alphaMin) best={Z,l,t.exponents,amin};
                        }
                        if (best.exponents.size()) hits.push_back(best);
                        break;
                    }
            }, bs.GetChild(i));
        lambdas.push_back(lmin);
        if (nDrops==0)
        {
            std::cout << "[basis trim] vet stage: ";
            if (r.trim.empty()) std::cout << "nothing to trim";
            else { std::cout << "trimmed "; r.trim.Write(std::cout); }
            std::cout << " -- no k-block of the full mesh drops an AO at orthoTol=" << orthoTol
                      << ";  min eig S over the k-mesh: ";
            if (lambdas.size()==1) std::cout << fmtLambda(lambdas[0]);
            else std::cout << "BEFORE " << fmtLambda(lambdas.front()) << "  ->  AFTER " << fmtLambda(lambdas.back());
            std::cout << std::endl;
            return r;
        }
        if (hits.empty()) throw std::logic_error("VetStageTrim: drops that map to no nameable shell");
        // THE cut of this pass: the globally most diffuse implicated shell.
        const Candidate* cut=&hits[0];
        for (const auto& h : hits) if (h.alphaMin<cut->alphaMin) cut=&h;
        if (r.trim.Removes(cut->Z, cut->l, cut->exponents))
            throw std::logic_error("VetStageTrim: the trimmed shell is still in the basis (the reader ignored the trim)");
        r.trim.shells.push_back({cut->Z, cut->l, cut->exponents});
        std::cout << "[basis trim] vet stage pass " << pass << ": min eig S " << fmtLambda(lmin) << "; " << nDrops
                  << " AO drop(s) over the k-mesh at orthoTol=" << orthoTol << " -> removing ";
        Gaussian::ShellTrim one; one.shells.push_back(r.trim.shells.back()); one.Write(std::cout);
        std::cout << " from EVERY site of that element, every k" << std::endl;
    }
    throw std::runtime_error("VetStageTrim: 64 passes without a clean basis -- something other than diffuse shells is near-dependent");
}

Complex_BS* GPWFactory(const ::qchem::Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                       const GPWParams& p)
{
    return new GPW_BasisSet(lat, std::move(mol), p);
}
Complex_BS* GPWFactory(const ::qchem::Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                       double densityEcut, rvec3_t kShift, Gaussian::CellImages images, double cutoffFactor)
{
    return new GPW_BasisSet(lat, std::move(mol),
                            GPWParams{densityEcut, cutoffFactor, PlaneWave::RasterPolicy::BallOnly, images, kShift});
}

} //namespace
