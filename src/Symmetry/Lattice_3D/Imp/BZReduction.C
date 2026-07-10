// File: Symmetry/Lattice_3D/Imp/BZReduction.C  IBZ reduction implementation.
module;
#include <vector>
module qchem.Symmetry.Lattice_3D.BZReduction;
import qchem.Math;   // fabs, lround

namespace qchem::Symmetry::Lattice_3D
{

double IBZMesh::WeightSum() const
{
    double s = 0.0;
    for (const auto& p : points) s += p.weight;
    return s;
}

//---------------------------------------------------------------------------------------
//  Grid index <-> linear index in KMesh loop order (ix outer, iz inner).
//
static int LinearIndex(const ivec3_t& i, const ivec3_t& N)
{
    return (i.x * N.y + i.y) * N.z + i.z;
}

//! Apply reciprocal op \a U to grid point \a i (with grid \a shift).  A symmetry maps
//! \f$k=(i+s)/N\f$ to \f$(i'+s)/N = U(i+s)/N\f$, so \f$i' = U(i+s)-s\f$ reduced mod \a N.
//! Returns false if \f$i'\f$ is not integral (the op does not map this grid onto itself).
static bool ApplyToGrid(const Matrix3D<double>& U, const ivec3_t& i,
                        const ivec3_t& N, const rvec3_t& shift, ivec3_t& out)
{
    rvec3_t ks(double(i.x) + shift.x, double(i.y) + shift.y, double(i.z) + shift.z);
    rvec3_t im = U * ks;                         // U(i+shift)
    double  t[3] = {im.x - shift.x, im.y - shift.y, im.z - shift.z};   // i' before mod
    int     n[3] = {N.x, N.y, N.z};
    int     r[3];
    for (int c = 0; c < 3; ++c)
    {
        long ii = lround(t[c]);
        if (fabs(t[c] - double(ii)) > 1e-6) return false;   // not grid-closed under this op
        r[c] = int(((ii % n[c]) + n[c]) % n[c]);            // reduce into [0,N)
    }
    out = ivec3_t(r[0], r[1], r[2]);
    return true;
}

//---------------------------------------------------------------------------------------
IBZMesh ReduceToIBZ(const ivec3_t& N, const rvec3_t& shift,
                    const std::vector<Matrix3D<double>>& ops)
{
    IBZMesh mesh;
    mesh.N = N;
    mesh.shift = shift;
    const size_t Ntot = size_t(N.x) * N.y * N.z;
    mesh.ownerOfGrid.assign(Ntot, -1);

    // Scan in KMesh order; the first unowned point of each star is its representative (and,
    // by scan order, the lowest-index member of the star).
    for (int ix = 0; ix < N.x; ++ix)
    for (int iy = 0; iy < N.y; ++iy)
    for (int iz = 0; iz < N.z; ++iz)
    {
        ivec3_t seed(ix, iy, iz);
        int seedLin = LinearIndex(seed, N);
        if (mesh.ownerOfGrid[seedLin] != -1) continue;

        // Build the star by orbit closure (BFS): the grid-closed subset of ops need not be
        // a subgroup, so iterate to a fixpoint rather than a single pass.
        std::vector<ivec3_t> orbit{seed};
        std::vector<int>     orbitLin{seedLin};
        for (size_t q = 0; q < orbit.size(); ++q)
            for (const auto& U : ops)
            {
                ivec3_t img;
                if (!ApplyToGrid(U, orbit[q], N, shift, img)) continue;
                int lin = LinearIndex(img, N);
                bool seen = false;
                for (int l : orbitLin) if (l == lin) { seen = true; break; }
                if (!seen) { orbit.push_back(img); orbitLin.push_back(lin); }
            }

        int repIdx = int(mesh.points.size());
        for (int l : orbitLin) mesh.ownerOfGrid[l] = repIdx;

        IBZPoint p;
        p.index    = seed;
        p.k        = rvec3_t((seed.x + shift.x) / N.x,
                             (seed.y + shift.y) / N.y,
                             (seed.z + shift.z) / N.z);
        p.starSize = int(orbit.size());
        p.weight   = double(orbit.size()) / double(Ntot);
        mesh.points.push_back(p);
    }
    return mesh;
}

} // namespace
