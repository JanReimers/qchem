// File: Structure/Lattice_3D/Imp/Imposition.C  ResolveImposition: the run's symmetry policy, moved verbatim from the GPW
// basis factory's DetectPointOps (V-BS-SYMMETRY S1).
module;
#include <iostream>   // the [symmetry] imposition line
#include <vector>
module qchem.Lattice_3D.Imposition;
import qchem.Symmetry.Lattice_3D.Fold;   // RequireClosedGroup (the imposition gate), MapsMeshOntoItself
import qchem.Matrix3D;                   // Transpose (the G-index scatter U = W^T)

namespace qchem
{

CrystalImposition ResolveImposition(const Lattice_3D& lat, const ImpositionOptions& o)
{
    namespace SL = qchem::Symmetry::Lattice_3D;
    const SL::SpaceGroup& sg = lat.GetSpaceGroup();

    CrystalImposition ops;
    ops.options = o;
    ops.policy.impose = o.imposeSymmetry ? SL::SymmetryPolicy::Impose::FullGroup
                                         : SL::SymmetryPolicy::Impose::None;
    if (!ops.policy.Imposes())
    {
        // KAxisPlan ruling 1 (user, 2026-10-27): a free run gets symmetric NUMERICS, not a symmetric answer.
        if (o.freeMeshInvariant)
            for (const auto& op : sg.DirectOps()) ops.meshInvariance.push_back({op.W, op.tau});
        return ops;
    }

    const bool magnetic = [&]{ for (int s : o.siteSpins) if (s!=0) return true; return false; }();
    if (magnetic)
    {
        ops.magneticDirect = lat.ShubnikovOps(o.siteSpins);
        size_t nFlip=0; for (const auto& op : ops.magneticDirect) if (op.sigma==SL::SpinAction::Flip) nFlip++;
        std::cout << "[symmetry] IMPOSING the SHUBNIKOV group of the declared ordering: "
                  << ops.magneticDirect.size() << " ops (" << ops.magneticDirect.size()-nFlip
                  << " spatial + " << nFlip << " spin-flip; the grey group would erase the order)" << std::endl;
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
        const ivec3_t N=lat.GetLimits();
        auto bad=[&](const Matrix3D<double>& U){ return !SL::MapsMeshOntoItself(U,N,o.kShift); };
        std::erase_if(ops.recipDensity,  [&](const auto& op){ return bad(op.U); });
        std::erase_if(ops.directDensity, [&](const auto& op){ return bad(Transpose(op.W)); });
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

} // namespace qchem
