// File: Structure/Lattice_3D/Imp/Imposition.C  ResolveImposition: the run's symmetry policy, moved verbatim from the GPW
// basis factory's DetectPointOps (V-BS-SYMMETRY S1).
module;
#include <cassert>
#include <iostream>   // the [symmetry] imposition line; the [XC mesh] / [Becke grid] lines
#include <format>     // the k shift on that line, independent of the stream's precision
#include <memory>
#include <vector>
module qchem.Lattice_3D.Imposition;
import qchem.Symmetry.Lattice_3D.Fold;   // RequireClosedGroup (the imposition gate), MapsMeshOntoItself
import qchem.Matrix3D;                   // Transpose (the G-index scatter U = W^T)
import qchem.SymmetrizeMesh;             // FoldMesh, kMeshMatchTol (the XC mesh's orbit fold)
import qchem.Reporting;                  // report::Timed -- the XC-mesh buckets (ParallelAndOraclePlan 1.1(b))

namespace qchem
{

CrystalImposition ResolveImposition(const Lattice_3D& lat, const ImpositionOptions& o)
{
    namespace SL = qchem::Symmetry::Lattice_3D;
    const SL::SpaceGroup& sg = lat.GetSpaceGroup();

    CrystalImposition ops;
    ops.options = o;
    ops.cell    = std::make_shared<const UnitCell>(lat.GetUnitCell());
    ops.policy.impose = o.imposeSymmetry ? SL::SymmetryPolicy::Impose::FullGroup
                                         : SL::SymmetryPolicy::Impose::None;
    if (!ops.policy.Imposes())
    {
        // KAxisPlan ruling 1 (user, 2026-10-27): a free run gets symmetric NUMERICS, not a symmetric answer.
        if (o.freeMeshInvariant)
            for (const auto& op : sg.DirectOps()) ops.meshInvariance.push_back({op.W, op.tau});
        return ops;
    }

    // THE INTERSECTION RULE (R1.0r, user ruling 2026-10-29): the imposed group is the detected group (spatial or Shubnikov)
    // INTERSECTED with the ops that map the k mesh onto itself.  The k-sum only samples the mesh, so star-averaging ρ over an op
    // the mesh does not respect symmetrizes it beyond its sampling: the density is then not the variational one, and its forces
    // are not dE/dR (a shifted-MP Si run read 1.02 mHa low before ecc9fe23).  The ops that map a mesh onto itself form a group,
    // so the intersection of two groups is a group and the imposition gate below still holds.  When it drops ops it reports
    // itself, one level into the caller's scope.
    const ivec3_t N=lat.GetLimits();
    auto offMesh=[&](const Matrix3D<double>& W){ return !SL::MapsMeshOntoItself(Transpose(W),N,o.kShift); };   // U = W^T
    auto reportDropped=[&](const char* group, size_t detected, size_t dropped)
    {
        if (dropped==0) return;
        std::cout << "[symmetry] the k mesh " << N.x << "x" << N.y << "x" << N.z << " (shift "
                  << std::format("{:g},{:g},{:g}",o.kShift.x,o.kShift.y,o.kShift.z) << ") keeps " << detected-dropped << " of the " << detected << " " << group
                  << " ops: the other " << dropped << " do not map it onto itself, and rho is symmetrized only under ops its sampling has"
                  << std::endl;
        report::Set("kMeshIntersection", report::json{{"group",group},{"detected",detected},{"kept",detected-dropped},{"dropped",dropped},
            {"kMesh",{N.x,N.y,N.z}},{"kShift",{o.kShift.x,o.kShift.y,o.kShift.z}},
            {"why","the imposed group is intersected with the ops that map the k mesh onto itself (R1.0r): a density symmetrized beyond its sampling is not variational"}});
    };

    const bool magnetic = [&]{ for (int s : o.siteSpins) if (s!=0) return true; return false; }();
    if (magnetic)
    {
        ops.magneticDirect = lat.ShubnikovOps(o.siteSpins);
        const size_t nShubnikov=ops.magneticDirect.size();
        const size_t nOff=std::erase_if(ops.magneticDirect, [&](const SL::SymOp& op){ return offMesh(op.W); });
        size_t nFlip=0; for (const auto& op : ops.magneticDirect) if (op.sigma==SL::SpinAction::Flip) nFlip++;
        std::cout << "[symmetry] IMPOSING the SHUBNIKOV group of the declared ordering: "
                  << ops.magneticDirect.size() << " ops (" << ops.magneticDirect.size()-nFlip
                  << " spatial + " << nFlip << " spin-flip; the grey group would erase the order)" << std::endl;
        reportDropped("Shubnikov", nShubnikov, nOff);
        // The density faces get ONLY the σ=None (sublattice-preserving) subgroup.  MEASURED (MnO run 37, 2026-08-11): they must
        // NOT get the flip ops' spatial parts -- the composite star-average is applied PER CHANNEL (each channel's
        // GetFourierDensity/GetRhoOnGrid symmetrizes ITSELF, and the ρ̃ MIXER consumes the channels separately through exactly
        // those faces), and a sublattice-swap op is not a symmetry of ONE channel: averaging under it equalizes the channels
        // and annihilates m MACHINE-EXACTLY at iteration 1 (the seed reads m1=+0.717, the first mixed density m1=6e-14).
        // Each channel IS invariant under the σ=None subgroup, so that is the correct per-channel projector; the flip content
        // is a PAIR property and is enforced where the pair is in hand (the spin-pair projection, the XC engine's (ρ,m)
        // star-average under the full σ-carrying set).  A PARTIAL projector is legal (folding merely reduced, never wrong).
        // ★ THE IMPOSITION GATE (doc/SymmetryUpgradePlan.md "SUPERCELLS", 2026-09-07): a star-average is a PROJECTOR only if
        // its op set is a GROUP, and nothing checked that until a Si 2x2x2 supercell diverged from its first Fock on a set with
        // 864 of 2304 products outside it.  Closure, NOT completeness -- the σ=None subgroup below is a legal partial projector.
        SL::RequireClosedGroup(ops.magneticDirect, "crystal imposition (Shubnikov group)");
        for (const auto& op : ops.magneticDirect)
            if (op.sigma==SL::SpinAction::None)
            {
                ops.directDensity.push_back({op.W, op.tau});
                ops.recipDensity .push_back({Transpose(op.W), op.tau});   // the G-index scatter U = W^T
            }
        // The PAIR projection needs every op to map the one block onto itself with no Bloch phase: a single unshifted Γ block.
        // A multi-k magnetic run keeps the partial (σ=None) projection the density faces carry.
        if (lat.GetLimits()==ivec3_t(1,1,1) && o.kShift==rvec3_t(0,0,0)) ops.spinPairOps=ops.magneticDirect;
        // k-fold: Γ-only for now -- the magnetic little-group bookkeeping at k≠Γ is the T3.4b-adjacent increment.  Leaving
        // recipFold empty means a multi-k imposed magnetic run keeps its FULL mesh (folding is merely absent, never wrong).
        return ops;
    }
    ops.recipFold     = sg.ReciprocalPointOps(/*timeReversal*/true, /*symmorphicOnly*/true);  // fold: linear + TR
    ops.recipDensity  = sg.ReciprocalOps();                                    // density G-space: {U|τ}
    ops.directDensity = sg.DirectOps();                                        // density raster: {W|τ}
    // The k-sum only respects the mesh-symmetry SUBGROUP (a shifted fcc MP mesh is not Oh-invariant), and the fold above
    // already uses just that subgroup.  Star-averaging rho over ops the mesh does not respect would sample a LARGER k-set
    // than the one whose bands are summed -- so impose the same subgroup on rho.
    {
        const size_t nDetected=ops.directDensity.size();
        std::erase_if(ops.recipDensity,  [&](const auto& op){ return offMesh(Transpose(op.U)); });
        const size_t nOff=std::erase_if(ops.directDensity, [&](const auto& op){ return offMesh(op.W); });
        assert(ops.recipDensity.size()==ops.directDensity.size());
        reportDropped("spatial", nDetected, nOff);
    }
    // ★ THE IMPOSITION GATE -- see the magnetic branch.  The reciprocal set is the transpose of this one, so closure of {W|τ}
    // settles both.  DirectOp carries no σ (a spatial group), so it lifts to SymOp with σ=None.
    {
        std::vector<SL::SymOp> spatial;
        spatial.reserve(ops.directDensity.size());
        for (const auto& op : ops.directDensity) spatial.push_back({op.W, op.tau});
        SL::RequireClosedGroup(spatial, "crystal imposition (spatial space group)");
    }
    return ops;
}

// Moved verbatim from tGPW_IBS::CreateXCQuadrature (V-BS-SYMMETRY S2): the δ-fit XC quadrature (doc/SymmetryUpgradePlan.md
// §6a).  On a §3-imposed run:
//  - Becke grid: the SITE-ADAPTED mesh (W2b) -- rep atoms carry site-group-invariant angular sets, partners the op-rotated
//    copies -- invariant BY CONSTRUCTION, no group-average growth pass;
//  - uniform grid (the cross cell): group-average it invariant (W1's MakeInvariant -- a no-op dedup when the grid is already
//    commensurate).
// Either way the geometry-fixed orbit fold feeds the engine's per-iteration rho star-average.  The CALLER opens the report
// scope ("xcQuadrature"), so every mesh built below reports one level into it.
qcMesh::FoldedMesh CrystalImposition::MakeXCMesh(const qcMesh::MeshParams& mp) const
{
    assert(cell && "CrystalImposition::MakeXCMesh: an imposition not made by ResolveImposition has no cell");
    const auto& dops = directDensity;
    if (dops.empty() && meshInvariance.empty())
        return {std::make_shared<const qcMesh::Mesh>(cell->CreateIntegrationMesh(mp)), {}};
    // FREE run: the same invariant mesh an imposed run gets, but NO fold -- nothing averages the density.
    if (dops.empty())
    {
        qcMesh::Mesh mesh = cell->CreateIntegrationMesh(mp, meshInvariance);
        std::cout << "[XC mesh] free run: mesh invariant under the " << meshInvariance.size() << " detected ops, "
                  << mesh.size() << " points; the density is NOT averaged" << std::endl;
        return {std::make_shared<const qcMesh::Mesh>(std::move(mesh)), {}};
    }

    // On a MAGNETICALLY imposed run (S3) the σ-carrying Shubnikov set drives the quadrature; its SPATIAL parts equal the
    // density raster ops (ResolveImposition filled both from one group), so the mesh/fold machinery is unchanged -- σ only ADDS
    // the per-op spin actions and the odd-field flags the engine's (ρ,m) channel-pair star-average consumes.
    std::vector<Symmetry::Lattice_3D::SymOp> ops;
    if (!magneticDirect.empty()) ops = magneticDirect;
    else { ops.reserve(dops.size()); for (const auto& op : dops) ops.push_back({op.W, op.tau}); }
    const double tol = qchem::kMeshMatchTol;   // the ONE mesh-point coincidence tolerance (R2.12)
    const Matrix3D<double>& A = cell->GetCellMatrix();

    // Ask the CELL for a mesh invariant under the ops.  WHICH strategy realises that (site-adapted by construction for Becke,
    // group-average for uniform) is the cell's business -- the imposition only needs the invariance, which is the T2
    // precondition for the fold below.
    qcMesh::Mesh mesh = cell->CreateIntegrationMesh(mp, ops);
    // ★ 1.1(b): the XC quadrature build is 23.75 s/call EXCLUSIVE of the mesh build and the Φ tables, and it is the largest
    // non-threading block in a GPW run.  These two buckets are the only candidates left inside it -- both are point-matching
    // passes over ~97k mesh points against 24 ops.
    auto fold = [&]{ qchem::report::Timed timed("setup: XC mesh orbit fold (FoldMesh)");
                     return FoldMesh(mesh, A, ops, tol); }();
    std::cout << "[Becke grid] imposed symmetry (" << ops.size() << " ops): invariant mesh "
              << mesh.size() << " points in " << fold.Reps()
              << " orbits; rho star-averaged each iteration (doc/SymmetryUpgradePlan.md 6a)" << std::endl;
    if (magneticDirect.empty())
        return {std::make_shared<const qcMesh::Mesh>(std::move(mesh)), std::move(fold)};

    // The σ tags (parallel to ops -- the fold's edge opIndex indexes them) and the odd-field zero flags: FRACTIONAL mesh
    // coordinates for the torus-metric fixed-point test.
    std::vector<Symmetry::SpinAction> sigmas;
    sigmas.reserve(ops.size());
    for (const auto& op : ops) sigmas.push_back(op.sigma);
    auto flags = [&]{
        qchem::report::Timed timed("setup: XC mesh Shubnikov flip-fixed test");
        std::vector<rvec3_t> frac;
        frac.reserve(mesh.size());
        const Matrix3D<double> Ainv = Invert(A);
        for (size_t i=0;i<mesh.size();++i) frac.push_back(Ainv*mesh.Points()[i]);
        return Symmetry::Lattice_3D::FlipFixedPointsPeriodic(frac, ops, tol); }();
    size_t nFixed=0; for (char c : flags) if (c) nFixed++;
    std::cout << "[Becke grid] Shubnikov: " << nFixed << " flip-fixed mesh points (m==0 there exactly)"
              << std::endl;
    return {std::make_shared<const qcMesh::Mesh>(std::move(mesh)), std::move(fold),
            std::move(sigmas), std::move(flags)};
}

} // namespace qchem
