// File: Symmetry/Lattice_3D/Imp/Fold.C  Orbit-fold implementation.
module;
#include <vector>
#include <utility>
#include <map>
#include <array>
#include <algorithm>
module qchem.Symmetry.Lattice_3D.Fold;
import qchem.Math;   // fabs, lround

namespace qchem::Symmetry::Lattice_3D
{

//---------------------------------------------------------------------------------------
//  The shared orbit-closure engine.  apply(o,i) returns the raw index of ops[o] acting on
//  point i, or -1 when the image is not in the set (that op is skipped for that point).
//  BFS to a fixpoint because the applicable op subset need not be a subgroup.  Seeds scan
//  in raw order, so each representative is the lowest raw index of its star -- the
//  ReduceToIBZ convention, preserved bit-identically.
//
template <class F> static Fold FoldByAction(int n, int nops, const F& apply)
{
    Fold f;
    f.owner.assign(n, -1);
    for (int seed = 0; seed < n; ++seed)
    {
        if (f.owner[seed] != -1) continue;

        // Orbit closure; a member discovered DIRECTLY from the seed (q==0) records its edge op on the
        // spot -- for a group action that is every member (the o-ascending scan gives the minimal op
        // index, identical to the fallback search below), so the per-member op search only runs for
        // members reachable solely by composition (non-closed hand-made op subsets).
        std::vector<int> orbit{seed};
        std::vector<int> edge{-1};
        for (size_t q = 0; q < orbit.size(); ++q)
            for (int o = 0; o < nops; ++o)
            {
                int img = apply(o, orbit[q]);
                if (img < 0) continue;
                bool seen = false;
                for (int l : orbit) if (l == img) { seen = true; break; }
                if (!seen) { orbit.push_back(img); edge.push_back(q == 0 ? o : -1); }
            }

        int rep = int(f.repRaw.size());
        f.repRaw.push_back(seed);
        f.starSize.push_back(int(orbit.size()));

        std::vector<std::pair<int,int>> edges;
        edges.reserve(orbit.size());
        for (size_t q = 0; q < orbit.size(); ++q)
        {
            int m = orbit[q];
            f.owner[m] = rep;
            int opIdx = edge[q];
            if (opIdx < 0)    // the rep itself, or a composition-only member: search (may stay -1)
                for (int o = 0; o < nops; ++o)
                    if (apply(o, seed) == m) { opIdx = o; break; }
            edges.emplace_back(m, opIdx);
        }
        f.members.push_back(std::move(edges));
    }
    return f;
}

//---------------------------------------------------------------------------------------
//  Grid path ({k},{G} grids): exact integer permutation of the periodic grid.
//
static int LinearIndex(const ivec3_t& i, const ivec3_t& N)
{
    return (i.x * N.y + i.y) * N.z + i.z;
}

//! Apply linear part \a U to grid point \a i (with grid \a shift): a symmetry maps
//! \f$k=(i+s)/N\f$ to \f$(i'+s)/N = U(i+s)/N\f$, so \f$i' = U(i+s)-s\f$ reduced mod \a N.
//! Returns false if \f$i'\f$ is not integral (the op does not map this grid onto itself).
static bool ApplyToGrid(const Matrix3D<double>& U, const ivec3_t& i,
                        const ivec3_t& N, const rvec3_t& shift, ivec3_t& out)
{
    rvec3_t ks(double(i.x) + shift.x, double(i.y) + shift.y, double(i.z) + shift.z);
    rvec3_t im = U * ks;                                                // U(i+shift)
    double  t[3] = {im.x - shift.x, im.y - shift.y, im.z - shift.z};    // i' before mod
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

Fold FoldGrid(const ivec3_t& N, const rvec3_t& shift, const std::vector<SymOp>& ops)
{
    const int Ntot = N.x * N.y * N.z;
    std::vector<ivec3_t> pts;
    pts.reserve(Ntot);
    for (int ix = 0; ix < N.x; ++ix)          // KMesh raw order: ix outer, iz inner
    for (int iy = 0; iy < N.y; ++iy)
    for (int iz = 0; iz < N.z; ++iz)
        pts.push_back(ivec3_t(ix, iy, iz));

    auto apply = [&](int o, int i) -> int
    {
        ivec3_t img;
        if (!ApplyToGrid(ops[o].W, pts[i], N, shift, img)) return -1;
        return LinearIndex(img, N);
    };
    return FoldByAction(Ntot, int(ops.size()), apply);
}

//---------------------------------------------------------------------------------------
//  Explicit G-list path: exact integer arithmetic, images matched by index lookup.
//
Fold FoldGVectors(const std::vector<ivec3_t>& ms, const std::vector<SymOp>& ops)
{
    const int n = int(ms.size());
    std::map<std::array<int,3>, int> index;
    for (int i = 0; i < n; ++i) index[{ms[i].x, ms[i].y, ms[i].z}] = i;

    auto apply = [&](int o, int i) -> int
    {
        rvec3_t im = ops[o].W * rvec3_t(double(ms[i].x), double(ms[i].y), double(ms[i].z));
        std::array<int,3> key;
        double t[3] = {im.x, im.y, im.z};
        for (int c = 0; c < 3; ++c)
        {
            long ii = lround(t[c]);
            if (fabs(t[c] - double(ii)) > 1e-6) return -1;   // non-integer image: not a G-index map
            key[c] = int(ii);
        }
        auto it = index.find(key);
        return it == index.end() ? -1 : it->second;
    };
    return FoldByAction(n, int(ops.size()), apply);
}

//---------------------------------------------------------------------------------------
//  Tolerance path ({r}/Becke): r -> W r + tau, images matched within Euclidean tol via a
//  sorted-by-x window search.
//
Fold FoldPoints(const std::vector<rvec3_t>& pts, const std::vector<SymOp>& ops, double tol)
{
    const int n = int(pts.size());
    std::vector<int> byx(n);
    for (int i = 0; i < n; ++i) byx[i] = i;
    std::sort(byx.begin(), byx.end(), [&](int a, int b) {return pts[a].x < pts[b].x;});

    auto lookup = [&](const rvec3_t& r) -> int
    {
        auto lo = std::lower_bound(byx.begin(), byx.end(), r.x - tol,
                                   [&](int a, double v) {return pts[a].x < v;});
        for (auto it = lo; it != byx.end() && pts[*it].x <= r.x + tol; ++it)
        {
            const rvec3_t& p = pts[*it];
            double dx = p.x - r.x, dy = p.y - r.y, dz = p.z - r.z;
            if (dx*dx + dy*dy + dz*dz <= tol*tol) return *it;
        }
        return -1;
    };
    auto apply = [&](int o, int i) -> int {return lookup(ops[o].W * pts[i] + ops[o].tau);};
    return FoldByAction(n, int(ops.size()), apply);
}

} // namespace
