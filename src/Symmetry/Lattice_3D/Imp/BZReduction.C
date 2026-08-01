// File: Symmetry/Lattice_3D/Imp/BZReduction.C  IBZ reduction implementation.
//
// Re-expressed on the shared orbit-fold primitive (doc/SymmetryUpgradePlan.md §7 step 2):
// FoldGrid scans in the same KMesh order with the same orbit closure, so the resulting
// IBZMesh is bit-identical to the original standalone reduction.
module;
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
    return mesh;
}

} // namespace
