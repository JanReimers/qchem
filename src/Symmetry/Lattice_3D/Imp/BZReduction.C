// File: Symmetry/Lattice_3D/Imp/BZReduction.C  IBZ reduction implementation.
//
// Re-expressed on the shared orbit-fold primitive (doc/SymmetryUpgradePlan.md §7 step 2):
// FoldGrid scans in the same KMesh order with the same orbit closure, so the resulting
// IBZMesh is bit-identical to the original standalone reduction.
module;
#include <string>
#include <stdexcept>
#include <cmath>
#include <vector>
module qchem.Symmetry.Lattice_3D.BZReduction;
import qchem.Symmetry.Lattice_3D.Fold;

namespace qchem::Symmetry::Lattice_3D
{

double IBZMesh::WeightSum() const
{
    double s = 0.0;
    for (const auto& p : points) s += p.weight;
    return s;
}

//---------------------------------------------------------------------------------------
IBZMesh ReduceToIBZ(const ivec3_t& N, const rvec3_t& shift,
                    const std::vector<Matrix3D<double>>& ops)
{
    std::vector<SymOp> sops;
    sops.reserve(ops.size());
    for (const auto& U : ops) sops.push_back({U, rvec3_t(0,0,0)});   // tau is phase-only on the k-side
    Fold f = FoldGrid(N, shift, sops);

    IBZMesh mesh;
    mesh.N = N;
    mesh.shift = shift;
    mesh.ownerOfGrid = std::move(f.owner);

    const size_t Ntot = size_t(N.x) * N.y * N.z;
    mesh.points.reserve(f.repRaw.size());
    for (size_t r = 0; r < f.repRaw.size(); ++r)
    {
        int lin = f.repRaw[r];                                       // KMesh linear order (ix outer, iz inner)
        ivec3_t idx(lin / (N.y * N.z), (lin / N.z) % N.y, lin % N.z);
        IBZPoint p;
        p.index    = idx;
        p.k        = rvec3_t((idx.x + shift.x) / N.x,
                             (idx.y + shift.y) / N.y,
                             (idx.z + shift.z) / N.z);
        p.starSize = f.starSize[r];
        p.weight   = double(f.starSize[r]) / double(Ntot);
        mesh.points.push_back(p);
    }

    // ★ THE PARTITION INVARIANT, ENFORCED (2026-09-09).  The stars must PARTITION the full k-grid, so
    // Sum_k w_k = Sum_r starSize_r / Ntot == 1 exactly (up to roundoff).  IBZMesh::WeightSum's own doc has
    // always said "should be 1" -- it just said it to nobody: the value was printed on the [IBZ] banner of
    // every multi-k run and never checked.
    //
    // WHY IT MATTERS MORE THAN IT LOOKS.  Every BZ-summed quantity carries these weights, so a weight sum
    // of W scales the electron count by W: the disabled GPW_SCF.SiliconMultiKPlumbing gate reports
    // charge=12 against 8 valence electrons, and its banner reads "Sum(w)=1.5" -- 8 x 1.5, exactly.  The
    // charge is not a stale anchor and cannot become one; 8 valence electrons is physics.  Caught here,
    // the diagnosis is one line at the source instead of a wrong energy several layers downstream.
    {
        const double W=mesh.WeightSum();
        if (std::abs(W-1.0) > 1e-12*double(Ntot))
            throw std::runtime_error("IBZMesh: the irreducible k-point weights sum to "+std::to_string(W)
                +", not 1.  The stars must PARTITION the "+std::to_string(Ntot)+"-point grid, so a sum "
                "above 1 means orbits OVERLAP (a point counted in more than one star) and a sum below 1 "
                "means points were dropped.  Every BZ-summed quantity scales by this factor -- the "
                "electron count first -- so continuing would produce a plausible-looking wrong answer.");
    }
    return mesh;
}

} // namespace
