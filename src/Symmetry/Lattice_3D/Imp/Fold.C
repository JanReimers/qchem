// File: Symmetry/Lattice_3D/Imp/Fold.C  Orbit-fold implementation.
module;
#include <vector>
#include <utility>
#include <map>
#include <unordered_map>
#include <array>
#include <algorithm>
module qchem.Symmetry.Lattice_3D.Fold;
import qchem.Math;   // fabs, lround, floor, sqrt

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

//---------------------------------------------------------------------------------------
//  Torus path (periodic {r}/Becke meshes): fractional coordinates matched modulo the lattice,
//  via a uniform bucket index (real meshes are 1e4-1e5 points x 48 ops -- lookups must be
//  sublinear; buckets wrap, so the boundary needs no special casing).
//
static double Wrap01(double x) {double r = x - floor(x); return r < 1.0 ? r : 0.0;}

static double TorusDelta(double a, double b)   // distance of two wrapped coords on the unit circle
{
    double d = fabs(a - b);
    return d <= 0.5 ? d : 1.0 - d;
}

class TorusIndex
{
public:
    TorusIndex(const std::vector<rvec3_t>& wrapped, double tol)
        : itsPts(wrapped), itsTol(tol)
    {
        // Bucket edge must exceed tol so a query only needs the 3x3x3 neighbourhood.
        itsNb = 64;
        while (itsNb > 1 && 1.0/itsNb <= 2.0*tol) itsNb /= 2;
        for (int i = 0; i < int(itsPts.size()); ++i)
            itsBuckets[Key(itsPts[i])].push_back(i);
    }

    //! Index of the mesh point within tol of \a f (torus metric), or -1.
    int Find(const rvec3_t& f) const
    {
        int bx = Bucket(f.x), by = Bucket(f.y), bz = Bucket(f.z);
        for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
        for (int dz = -1; dz <= 1; ++dz)
        {
            auto it = itsBuckets.find(Key((bx+dx+itsNb)%itsNb, (by+dy+itsNb)%itsNb, (bz+dz+itsNb)%itsNb));
            if (it == itsBuckets.end()) continue;
            for (int i : it->second)
            {
                const rvec3_t& p = itsPts[i];
                double ex = TorusDelta(p.x, f.x), ey = TorusDelta(p.y, f.y), ez = TorusDelta(p.z, f.z);
                if (ex*ex + ey*ey + ez*ez <= itsTol*itsTol) return i;
            }
        }
        return -1;
    }

private:
    int       Bucket(double x) const {return int(Wrap01(x)*itsNb) % itsNb;}
    long long Key(int bx, int by, int bz) const {return (static_cast<long long>(bx)*itsNb + by)*itsNb + bz;}
    long long Key(const rvec3_t& p) const {return Key(Bucket(p.x), Bucket(p.y), Bucket(p.z));}

    const std::vector<rvec3_t>& itsPts;
    double itsTol;
    int    itsNb;
    std::unordered_map<long long, std::vector<int>> itsBuckets;
};

Fold FoldPointsPeriodic(const std::vector<rvec3_t>& pts, const std::vector<SymOp>& ops, double tol)
{
    const int n = int(pts.size());
    std::vector<rvec3_t> wrapped;
    wrapped.reserve(n);
    for (const auto& p : pts) wrapped.push_back(rvec3_t(Wrap01(p.x), Wrap01(p.y), Wrap01(p.z)));
    TorusIndex index(wrapped, tol);

    auto apply = [&](int o, int i) -> int
    {
        rvec3_t img = ops[o].W * wrapped[i] + ops[o].tau;
        return index.Find(rvec3_t(Wrap01(img.x), Wrap01(img.y), Wrap01(img.z)));
    };
    return FoldByAction(n, int(ops.size()), apply);
}

std::vector<int> CountUnmappedPeriodic(const std::vector<rvec3_t>& pts,
                                       const std::vector<SymOp>& ops, double tol)
{
    std::vector<rvec3_t> wrapped;
    wrapped.reserve(pts.size());
    for (const auto& p : pts) wrapped.push_back(rvec3_t(Wrap01(p.x), Wrap01(p.y), Wrap01(p.z)));
    TorusIndex index(wrapped, tol);

    std::vector<int> bad(ops.size(), 0);
    for (size_t o = 0; o < ops.size(); ++o)
        for (const auto& f : wrapped)
        {
            rvec3_t img = ops[o].W * f + ops[o].tau;
            if (index.Find(rvec3_t(Wrap01(img.x), Wrap01(img.y), Wrap01(img.z))) < 0) ++bad[o];
        }
    return bad;
}

// The odd-field fixed-point audit (S2): a point some Flip op maps onto itself.  No point-set index
// needed -- the question is per point, against its OWN image (torus metric).
std::vector<char> FlipFixedPointsPeriodic(const std::vector<rvec3_t>& pts,
                                          const std::vector<SymOp>& ops, double tol)
{
    std::vector<char> fixed(pts.size(), 0);
    for (size_t i = 0; i < pts.size(); ++i)
    {
        const rvec3_t f(Wrap01(pts[i].x), Wrap01(pts[i].y), Wrap01(pts[i].z));
        for (const auto& op : ops)
        {
            if (op.sigma != SpinAction::Flip) continue;
            rvec3_t g = op.W * f + op.tau;
            double ex = TorusDelta(Wrap01(g.x), f.x), ey = TorusDelta(Wrap01(g.y), f.y),
                   ez = TorusDelta(Wrap01(g.z), f.z);
            if (ex*ex + ey*ey + ez*ez <= tol*tol) { fixed[i] = 1; break; }
        }
    }
    return fixed;
}

} // namespace
