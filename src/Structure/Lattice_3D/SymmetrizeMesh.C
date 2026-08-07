// File: Structure/Lattice_3D/SymmetrizeMesh.C  Symmetry folding of quadrature meshes (T2 groundwork).
//
// doc/SymmetryUpgradePlan.md §2a: qcSymmetry owns the FOLD (a pure partition over primitive vectors);
// the Mesh-typed convenience lives CLIENT-side.  It sits HERE, in qcStructure, where Mesh and Symmetry
// MEET WITHOUT DEPENDING ON EACH OTHER (user decision 2026-08-01): qcMesh stays a pure quadrature leaf,
// qcSymmetry stays the pure group-theory leaf, and qcStructure -- which links both, builds the periodic
// Becke mesh (UnitCell), and owns the crystal's detected SpaceGroup (Lattice_3D) -- is where quadrature
// meets group action.  (The G-space sibling SymmetrizeGMap lives with the G-currency in qcBasisSet.)
//
// The {r}-quadrature toolkit for the T2 site-adapted Becke grid, BEFORE any XC wiring:
//   - FoldMesh       : the partition itself (consumers needing per-member data, e.g. SymmetrizeValues);
//   - Symmetrize     : the compact folded mesh -- one representative per star, weight = the SUM of the
//                      member weights (exact for a totally symmetric integrand whatever the weight
//                      distribution, == w*n_star when the mesh is op-invariant);
//   - UnmatchedCounts: the INVARIANCE checker -- per op, how many points fail to land on a mesh point.
//                      A pointwise star-average / fold is exact only on an invariant mesh (§5/§6 T2);
//   - MakeInvariant  : the CONSTRUCTIVE half -- group-average any mesh into an invariant one (the
//                      union of all op images at weight w/|ops|, coincident points merged).  Each
//                      rotated copy is a quadrature of the same polynomial exactness, so the average
//                      keeps the original exactness degree AND gains exact op-invariance; folding it
//                      back (Symmetrize) returns to ~the original point count.  The generic route to
//                      the site-group-adapted Becke grid before any bespoke minimal-direction design.
//
// Mesh points are CARTESIAN (Mesh's currency); the ops act on FRACTIONAL coordinates on the unit
// torus, so every entry point takes the cell matrix A (columns = lattice vectors, r = A f) -- or a
// SpaceGroup, which carries both its ops and its cell (§2c: no ops at the call site).
module;
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include <vector>
export module qchem.SymmetrizeMesh;
export import qchem.Mesh;                            // qcMesh::Mesh
export import qchem.Symmetry.Lattice_3D.Fold;        // Fold, SymOp, FoldPointsPeriodic, SymmetrizeValues
export import qchem.Symmetry.Lattice_3D.SpaceGroup;  // SpaceGroup (the ops-encapsulated overloads)
import qchem.Math;                                   // floor, fabs

export namespace qchem
{

namespace SL = qchem::Symmetry::Lattice_3D;

//! \brief Default coincidence tolerance for mesh-point matching, in FRACTIONAL coordinates (the torus
//! metric these routines work in), so it is cell-size independent by construction.
//!
//! Named once rather than repeated as a literal at each entry point: every routine here decides the SAME
//! question -- "are these two mesh points the same point?" -- so they must agree, and a reader should not
//! have to diff four default arguments to confirm that they do (R2.12).  1e-8 is far below any real mesh
//! spacing and far above the round-off of a fractional coordinate that has been through one W matrix.
inline constexpr double kMeshMatchTol = 1e-8;

//! The fold partition of \a m's points under \f$\{W|\tau\}\f$ ops acting on fractional coordinates
//! (torus metric, tolerance \a tol in fractional units).
inline SL::Fold FoldMesh(const qcMesh::Mesh& m, const Matrix3D<double>& A,
                         const std::vector<SL::SymOp>& ops, double tol);

//! The compact folded mesh: representatives only, weight = the sum over the star's members --
//! exact for a totally symmetric integrand:
//! \f$\sum_i w_i f(r_i) = \sum_\mathrm{stars} f(r_\mathrm{rep})\sum_\mathrm{members} w\f$.
inline qcMesh::Mesh Symmetrize(const qcMesh::Mesh& m, const Matrix3D<double>& A,
                               const std::vector<SL::SymOp>& ops, double tol);

//! Per-op invariance defect: how many of \a m's points an op maps OFF the mesh (all zeros ==
//! the mesh is invariant, the §6 T2 precondition for pointwise star-averaging and folding).
inline std::vector<int> UnmatchedCounts(const qcMesh::Mesh& m, const Matrix3D<double>& A,
                                        const std::vector<SL::SymOp>& ops, double tol);

//! Group-average \a m into an op-invariant mesh: the union of all op images, each carrying weight
//! \f$w/|\mathrm{ops}|\f$, coincident points (torus metric, \a tol) merged.  Pass the FULL op set
//! including the identity -- the average is a projector only over a group.  Preserves \f$\sum w\f$
//! and the quadrature's polynomial exactness; the point count grows toward \f$|\mathrm{ops}|\times\f$
//! for a generic mesh (\c Symmetrize folds it back to ~the original count).
inline qcMesh::Mesh MakeInvariant(const qcMesh::Mesh& m, const Matrix3D<double>& A,
                                  const std::vector<SL::SymOp>& ops, double tol);

// ---- The ops-encapsulated overloads (§2c): the SpaceGroup supplies its own {W|τ} set + cell. ----
inline SL::Fold         FoldMesh       (const qcMesh::Mesh& m, const SL::SpaceGroup& sg, double tol=kMeshMatchTol);
inline qcMesh::Mesh     Symmetrize     (const qcMesh::Mesh& m, const SL::SpaceGroup& sg, double tol=kMeshMatchTol);
inline std::vector<int> UnmatchedCounts(const qcMesh::Mesh& m, const SL::SpaceGroup& sg, double tol=kMeshMatchTol);
inline qcMesh::Mesh     MakeInvariant  (const qcMesh::Mesh& m, const SL::SpaceGroup& sg, double tol=kMeshMatchTol);

//! \brief SITE-GROUP-INVARIANT angular quadrature of polynomial degree \a L (doc/SymmetryUpgradePlan.md
//! §6a W2 -- "the site group fixes both how many and which directions").  Directions = unions of
//! orbits of directions under the site ops, so the grid is op-invariant BY CONSTRUCTION with equal
//! weights within each orbit.  THE MIXED RULE (§6a W2 refinement): the candidate pool holds the
//! group's SPECIAL orbits first (op axes + mirror-plane directions -- the small Lebedev-like orbits),
//! filtered against \a avoid (the ATOM's actual bond directions -- only some special families are
//! poison, e.g. diamond T_d's bonded <111> tetrahedron, while <100>/<110>/the anti-bond tetrahedron
//! point between bonds and are safe), then generic Fibonacci-sphere seeds.  Orbit weights solve the
//! spherical-harmonic moment conditions \f$Q(\bar Y_{lm}) = \sqrt{4\pi}\,\delta_{l0}\,\delta_{m0},\
//! l\le L\f$ by NNLS (Lawson-Hanson, per-DIRECTION cost-biased column selection), which guarantees
//! nonnegative weights and picks a SPARSE orbit support -- plain least squares was measured to leave
//! a ~-5e-4 negative-weight floor at L=29 however many orbits were added.  Only the positive-weight
//! support is emitted; \f$\sum w = 4\pi\f$ (the \f$Y_{00}\f$ condition).  Deterministic.
//! \param ops   the site group as CARTESIAN orthogonal matrices (rotations + reflections; include E).
//!              Fractional \f$W\f$ ops convert via \f$A W A^{-1}\f$ upstream.
//! \param L     the polynomial exactness degree (GL-\a L quality class).
//! \param avoid unit directions (e.g. the atom's bonds) no quadrature point may lie on.
qcMesh::Mesh MakeInvariantAngularMesh(const std::vector<Matrix3D<double>>& ops, int L,
                                      const std::vector<rvec3_t>& avoid = {});

} // export namespace

//-----------------------------------------------------------------------------------------------
//  Implementation (module-private helpers + the inline definitions).
//
namespace qchem
{

static double Wrap01(double x) {double r = x - floor(x); return r < 1.0 ? r : 0.0;}

//! Cartesian mesh points -> wrapped fractional coordinates (f = A^{-1} r).
static std::vector<rvec3_t> Fractional(const qcMesh::Mesh& m, const Matrix3D<double>& A)
{
    const Matrix3D<double> Ainv = Invert(A);
    std::vector<rvec3_t> f;
    f.reserve(m.size());
    for (size_t i = 0; i < m.size(); ++i)   // indexed: Blaze iterator ops are not visible cross-module
    {
        rvec3_t x = Ainv * m.Points()[i];
        f.push_back(rvec3_t(Wrap01(x.x), Wrap01(x.y), Wrap01(x.z)));
    }
    return f;
}

static std::vector<SL::SymOp> OpsOf(const SL::SpaceGroup& sg)
{
    std::vector<SL::SymOp> ops;
    ops.reserve(sg.Ops().size());
    for (const auto& op : sg.Ops()) ops.push_back({op.W, op.tau});
    return ops;
}

inline SL::Fold FoldMesh(const qcMesh::Mesh& m, const Matrix3D<double>& A,
                         const std::vector<SL::SymOp>& ops, double tol)
{
    return SL::FoldPointsPeriodic(Fractional(m, A), ops, tol);
}

inline qcMesh::Mesh Symmetrize(const qcMesh::Mesh& m, const Matrix3D<double>& A,
                               const std::vector<SL::SymOp>& ops, double tol)
{
    SL::Fold f = FoldMesh(m, A, ops, tol);
    rvec3vec_t pts(f.Reps());
    rvec_t     w(f.Reps());
    for (size_t r = 0; r < f.Reps(); ++r)
    {
        pts[r] = m.Points()[f.repRaw[r]];            // the representative keeps its ORIGINAL Cartesian point
        double s = 0.0;
        for (auto [mi, o] : f.members[r]) s += m.Weights()[mi];
        w[r] = s;                                    // summed member weights: exact for a symmetric integrand
    }
    return qcMesh::Mesh(std::move(pts), std::move(w));
}

inline std::vector<int> UnmatchedCounts(const qcMesh::Mesh& m, const Matrix3D<double>& A,
                                        const std::vector<SL::SymOp>& ops, double tol)
{
    return SL::CountUnmappedPeriodic(Fractional(m, A), ops, tol);
}

//  MakeInvariant's growing point set needs a torus-metric merge with sublinear lookup: a uniform
//  bucket hash on the wrapped fractional cube (buckets wrap, edge > 2*tol so a query only visits
//  the 3x3x3 neighbourhood).
class MergeIndex
{
public:
    MergeIndex(double tol) : itsTol(tol)
    {
        itsNb = 64;
        while (itsNb > 1 && 1.0/itsNb <= 2.0*tol) itsNb /= 2;
    }
    int Find(const std::vector<rvec3_t>& pts, const rvec3_t& q) const
    {
        int bx = Bucket(q.x), by = Bucket(q.y), bz = Bucket(q.z);
        for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
        for (int dz = -1; dz <= 1; ++dz)
        {
            auto it = itsBuckets.find(Key((bx+dx+itsNb)%itsNb, (by+dy+itsNb)%itsNb, (bz+dz+itsNb)%itsNb));
            if (it == itsBuckets.end()) continue;
            for (int i : it->second)
            {
                double ex = Delta(pts[i].x, q.x), ey = Delta(pts[i].y, q.y), ez = Delta(pts[i].z, q.z);
                if (ex*ex + ey*ey + ez*ez <= itsTol*itsTol) return i;
            }
        }
        return -1;
    }
    void Insert(const rvec3_t& p, int idx) {itsBuckets[Key(Bucket(p.x), Bucket(p.y), Bucket(p.z))].push_back(idx);}

private:
    static double Delta(double a, double b) {double d = fabs(a - b); return d <= 0.5 ? d : 1.0 - d;}
    int       Bucket(double x) const {return int(Wrap01(x)*itsNb) % itsNb;}
    long long Key(int bx, int by, int bz) const {return (static_cast<long long>(bx)*itsNb + by)*itsNb + bz;}

    double itsTol;
    int    itsNb;
    std::unordered_map<long long, std::vector<int>> itsBuckets;
};

inline qcMesh::Mesh MakeInvariant(const qcMesh::Mesh& m, const Matrix3D<double>& A,
                                  const std::vector<SL::SymOp>& ops, double tol)
{
    const double wScale = 1.0/double(ops.size());

    std::vector<rvec3_t> outF;      // wrapped fractional points of the growing invariant mesh
    std::vector<double>  outW;
    MergeIndex           index(tol);

    std::vector<rvec3_t> f = Fractional(m, A);
    for (size_t i = 0; i < f.size(); ++i)
        for (const auto& op : ops)
        {
            rvec3_t img = op.W * f[i] + op.tau;
            img = rvec3_t(Wrap01(img.x), Wrap01(img.y), Wrap01(img.z));
            int j = index.Find(outF, img);
            if (j >= 0) outW[j] += wScale*m.Weights()[i];
            else
            {
                index.Insert(img, int(outF.size()));
                outF.push_back(img);
                outW.push_back(wScale*m.Weights()[i]);
            }
        }

    rvec3vec_t pts(outF.size());
    rvec_t     w(outW.size());
    for (size_t i = 0; i < outF.size(); ++i) { pts[i] = A * outF[i]; w[i] = outW[i]; }
    return qcMesh::Mesh(std::move(pts), std::move(w));
}

inline SL::Fold FoldMesh(const qcMesh::Mesh& m, const SL::SpaceGroup& sg, double tol)
    {return FoldMesh(m, sg.CellMatrix(), OpsOf(sg), tol);}
inline qcMesh::Mesh Symmetrize(const qcMesh::Mesh& m, const SL::SpaceGroup& sg, double tol)
    {return Symmetrize(m, sg.CellMatrix(), OpsOf(sg), tol);}
inline std::vector<int> UnmatchedCounts(const qcMesh::Mesh& m, const SL::SpaceGroup& sg, double tol)
    {return UnmatchedCounts(m, sg.CellMatrix(), OpsOf(sg), tol);}
inline qcMesh::Mesh MakeInvariant(const qcMesh::Mesh& m, const SL::SpaceGroup& sg, double tol)
    {return MakeInvariant(m, sg.CellMatrix(), OpsOf(sg), tol);}

//-----------------------------------------------------------------------------------------------
//  MakeInvariantAngularMesh (§6a W2): generic-orbit invariant angular quadrature.
//

//! All (L+1)^2 ORTHONORMAL real spherical harmonics at unit vector \a n, packed l-major
//! (l=0..L; within l: m=0, then (cos,sin) pairs m=1..l).  Standard stable normalized
//! associated-Legendre recurrences (good far beyond L=30).
static void EvalRealYlm(const rvec3_t& n, int L, std::vector<double>& out)
{
    const double ct = n.z;                                // cos(theta)
    const double st = sqrt(n.x*n.x + n.y*n.y);            // sin(theta) >= 0
    const double phi = atan2(n.y, n.x);

    // Pbar[m][l] built column-by-column in m; store flat as Pbar(l,m) -> p[l*(L+1)+m].
    std::vector<double> p((L+1)*(L+1), 0.0);
    p[0] = sqrt(1.0/(4.0*Pi));                            // Pbar_0^0
    for (int m = 1; m <= L; ++m)                          // diagonal: Pbar_m^m
        p[m*(L+1)+m] = -sqrt((2.0*m+1.0)/(2.0*m)) * st * p[(m-1)*(L+1)+(m-1)];
    for (int m = 0; m < L; ++m)                           // first off-diagonal: Pbar_{m+1}^m
        p[(m+1)*(L+1)+m] = sqrt(2.0*m+3.0) * ct * p[m*(L+1)+m];
    for (int m = 0; m <= L; ++m)                          // upward in l
        for (int l = m+2; l <= L; ++l)
        {
            const double a = sqrt((4.0*l*l-1.0)/(double(l)*l - double(m)*m));
            const double b = sqrt((double(l-1)*(l-1) - double(m)*m)/(4.0*double(l-1)*(l-1) - 1.0));
            p[l*(L+1)+m] = a*(ct*p[(l-1)*(L+1)+m] - b*p[(l-2)*(L+1)+m]);
        }

    out.assign((L+1)*(L+1), 0.0);
    size_t j = 0;
    const double r2 = sqrt(2.0);
    for (int l = 0; l <= L; ++l)
    {
        out[j++] = p[l*(L+1)];                            // m = 0
        for (int m = 1; m <= l; ++m)
        {
            out[j++] = r2 * p[l*(L+1)+m] * cos(m*phi);
            out[j++] = r2 * p[l*(L+1)+m] * sin(m*phi);
        }
    }
}

//! Solve the (small, SPD after a tiny ridge) normal equations G w = g by Cholesky in place.
static std::vector<double> SolveSPD(std::vector<double> G, std::vector<double> g)
{
    const int K = int(g.size());
    double tr = 0.0;
    for (int i = 0; i < K; ++i) tr += G[i*K+i];
    for (int i = 0; i < K; ++i) G[i*K+i] += 1e-14*tr;     // ridge: rank-deficiency guard
    for (int i = 0; i < K; ++i)                            // Cholesky G = L L^T
    {
        for (int k = 0; k < i; ++k) G[i*K+i] -= G[i*K+k]*G[i*K+k];
        G[i*K+i] = sqrt(G[i*K+i]);
        for (int jj = i+1; jj < K; ++jj)
        {
            for (int k = 0; k < i; ++k) G[jj*K+i] -= G[jj*K+k]*G[i*K+k];
            G[jj*K+i] /= G[i*K+i];
        }
    }
    for (int i = 0; i < K; ++i)                            // forward
    {
        for (int k = 0; k < i; ++k) g[i] -= G[i*K+k]*g[k];
        g[i] /= G[i*K+i];
    }
    for (int i = K-1; i >= 0; --i)                         // backward
    {
        for (int k = i+1; k < K; ++k) g[i] -= G[k*K+i]*g[k];
        g[i] /= G[i*K+i];
    }
    return g;
}

//! The axis of a (non-identity, non-inversion) orthogonal op: the eigenvector for eigenvalue +1
//! (proper rotation) / -1 (mirror or rotoreflection), i.e. a null vector of \f$R \mp I\f$ -- taken
//! as the largest cross product of two of its rows.  False for E and the inversion (no axis).
static bool OpAxis(const Matrix3D<double>& R, rvec3_t& axis)
{
    const double s = (Determinant(R) > 0) ? -1.0 : +1.0;
    Matrix3D<double> M = R;
    M(1,1) += s; M(2,2) += s; M(3,3) += s;
    rvec3_t r1(M(1,1),M(1,2),M(1,3)), r2(M(2,1),M(2,2),M(2,3)), r3(M(3,1),M(3,2),M(3,3));
    rvec3_t best = Cross(r1,r2);
    rvec3_t c    = Cross(r2,r3); if (c*c > best*best) best = c;
    c            = Cross(r3,r1); if (c*c > best*best) best = c;
    if (best*best < 1e-16) return false;
    axis = best*(1.0/sqrt(best*best));
    return true;
}

//! Lawson-Hanson NNLS on the angular moment system \f$A w = \sqrt{4\pi}\,e_{Y_{00}},\ w \ge 0\f$
//! (\a A row-major nY x K).  The entering column maximizes the gradient PER DIRECTION (\a cost =
//! orbit size), so cheap special orbits win over generic full-size orbits when comparable -- the
//! Lebedev-efficiency bias of the §6a mixed rule.  True when the sup-norm residual reaches \a tol;
//! the returned \a w is nonnegative with a sparse (<= #independent moment conditions) support.
static bool SolveNNLS(const std::vector<double>& A, int nY, int K, const std::vector<double>& cost,
                      std::vector<double>& w, double tol)
{
    const double b0 = sqrt(4.0*Pi);
    w.assign(K, 0.0);
    std::vector<int>    P;                                // passive (support) set
    std::vector<char>   inP(K, 0);
    std::vector<double> r(nY);

    auto solveP = [&]() -> std::vector<double>            // unconstrained LS on the passive set
    {
        const int m = int(P.size());
        std::vector<double> G(size_t(m)*m, 0.0), g(m, 0.0);
        for (int a = 0; a < m; ++a)
        {
            g[a] = A[size_t(0)*K + P[a]]*b0;              // b is nonzero only on the Y00 row
            for (int bb = a; bb < m; ++bb)
            {
                double s = 0.0;
                for (int j = 0; j < nY; ++j) s += A[size_t(j)*K+P[a]]*A[size_t(j)*K+P[bb]];
                G[size_t(a)*m+bb] = G[size_t(bb)*m+a] = s;
            }
        }
        return SolveSPD(std::move(G), std::move(g));
    };

    for (int iter = 0; iter < 3*K + 100; ++iter)
    {
        for (int j = 0; j < nY; ++j)                      // residual r = b - A w (drift-free recompute)
        {
            double s = (j == 0) ? b0 : 0.0;
            for (int k : P) s -= A[size_t(j)*K+k]*w[k];
            r[j] = s;
        }
        double rmax = 0.0;
        for (double x : r) rmax = std::max(rmax, fabs(x));
        if (rmax < tol) return true;

        int jin = -1; double best = 0.0;
        for (int k = 0; k < K; ++k)
        {
            if (inP[k]) continue;
            double g = 0.0;
            for (int j = 0; j < nY; ++j) g += A[size_t(j)*K+k]*r[j];
            if (g/cost[k] > best) {best = g/cost[k]; jin = k;}
        }
        if (jin < 0) return false;                        // gradient exhausted: the pool is too small
        P.push_back(jin); inP[jin] = 1;

        for (;;)                                          // inner: restore strict positivity on P
        {
            std::vector<double> z = solveP();
            double zmin = 1e300;
            for (double x : z) zmin = std::min(zmin, x);
            if (zmin > 0.0) {for (size_t a = 0; a < P.size(); ++a) w[P[a]] = z[a]; break;}
            double alpha = 1e300;                         // step to the first variable hitting zero
            for (size_t a = 0; a < P.size(); ++a)
                if (z[a] <= 0.0) alpha = std::min(alpha, w[P[a]]/(w[P[a]] - z[a]));
            for (size_t a = 0; a < P.size(); ++a) w[P[a]] += alpha*(z[a] - w[P[a]]);
            std::vector<int> Pn;
            for (int k : P) {if (w[k] > 1e-14) Pn.push_back(k); else {w[k] = 0.0; inP[k] = 0;}}
            P.swap(Pn);
            if (P.empty()) break;
        }
    }
    return false;
}

qcMesh::Mesh MakeInvariantAngularMesh(const std::vector<Matrix3D<double>>& ops, int L,
                                      const std::vector<rvec3_t>& avoid)
{
    const int    nG   = int(ops.size());
    const double dtol = 1e-9;                             // direction-coincidence tolerance
    const int    nY   = (L+1)*(L+1);

    // Deterministic generic seed pool: Fibonacci sphere (never lands on symmetry axes for
    // generic counts), scanned in order.
    auto seed = [](int i)
    {
        const double z   = 1.0 - 2.0*(i + 0.5)/1024.0;
        const double phi = i * (Pi*(3.0 - sqrt(5.0)));
        const double s   = sqrt(1.0 - z*z);
        return rvec3_t(s*cos(phi), s*sin(phi), z);
    };
    auto dist2 = [](const rvec3_t& a, const rvec3_t& b)
    { double dx=a.x-b.x, dy=a.y-b.y, dz=a.z-b.z; return dx*dx+dy*dy+dz*dz; };

    // The orbit of a candidate direction; empty when a member lands ON an avoided (bond) direction
    // (exact coincidence class: the site ops permute the bonds, so a hit is a hit to full precision).
    auto orbitOf = [&](const rvec3_t& d0) -> std::vector<rvec3_t>
    {
        std::vector<rvec3_t> orb;
        for (const auto& R : ops)
        {
            rvec3_t im = R*d0;
            bool dup = false;
            for (const auto& q : orb) if (dist2(im,q) < dtol*dtol) {dup = true; break;}
            if (!dup) orb.push_back(im);
        }
        for (const auto& m : orb)
            for (const auto& b : avoid)
                if (m*b > 1.0 - 1e-6) return {};
        return orb;
    };

    std::vector<std::vector<rvec3_t>> orbits;             // the accepted candidate pool, cheap-first
    auto tryAdd = [&](std::vector<rvec3_t> orb) -> bool   // reject empties and near-duplicate orbits
    {
        if (orb.empty()) return false;
        for (const auto& o : orbits)
            for (const auto& q : o)
                for (const auto& m : orb)
                    if (dist2(m,q) < 1e-8) return false;
        orbits.push_back(std::move(orb));
        return true;
    };

    // ---- THE MIXED RULE (§6a W2 refinement): special-orbit candidates first, smallest first. ----
    // (i) op axes (+-n: opposite special orbits are distinct without inversion -- e.g. T_d's bonded
    // vs anti-bond <111> tetrahedra, only one of which the avoid list kills);
    // (ii) mirror-plane directions (seeds projected into each mirror plane: |G|/2-size semi-special
    // orbits, steerable within the plane).
    {
        std::vector<std::vector<rvec3_t>> cand;
        for (const auto& R : ops)
        {
            rvec3_t n;
            if (!OpAxis(R, n)) continue;
            auto o = orbitOf(n);                          if (!o.empty()) cand.push_back(std::move(o));
            o      = orbitOf(rvec3_t(-n.x,-n.y,-n.z));    if (!o.empty()) cand.push_back(std::move(o));
            if (Determinant(R) < 0.0 && fabs(R(1,1)+R(2,2)+R(3,3) - 1.0) < 1e-9)
                for (int t = 0; t < 6; ++t)               // a mirror: a few in-plane candidates
                {
                    rvec3_t s = seed((t*401) & 1023);
                    rvec3_t v = s - n*(s*n);
                    if (v*v < 0.1) continue;              // seed too close to the normal
                    o = orbitOf(v*(1.0/sqrt(v*v)));
                    if (!o.empty()) cand.push_back(std::move(o));
                }
        }
        std::stable_sort(cand.begin(), cand.end(),
                         [](const auto& a, const auto& b) {return a.size() < b.size();});
        for (auto& o : cand) tryAdd(std::move(o));
    }
    const int nSpecial = int(orbits.size());

    // ---- Generic full-size orbits appended on demand. ----
    int nextSeed = 0;
    auto addGeneric = [&]() -> bool
    {
        for (; nextSeed < 1024; ++nextSeed)
        {
            // Stride the Fibonacci indices (401 coprime to 1024): consecutive CANDIDATES then jump
            // around the sphere, so early orbits cover it instead of clustering at the pole (which
            // left the moment system near-degenerate however many orbits were added).
            std::vector<rvec3_t> orb = orbitOf(seed((nextSeed*401) & 1023));
            if (int(orb.size()) != nG) continue;          // special direction: not generic, skip
            if (!tryAdd(std::move(orb))) continue;
            ++nextSeed;
            return true;
        }
        return false;
    };

    // Grow the pool from the invariant-subspace estimate and NNLS-solve until the moment conditions
    // are met (the solver guarantees nonnegativity; only the positive support is emitted).
    std::vector<double> w, Y;
    for (int K = std::max(nSpecial, 2*(nY/nG) + 8); ; K = 2*K)
    {
        while (int(orbits.size()) < K && addGeneric()) {}
        const int nOrb = int(orbits.size());

        // Rows: sum of Ylm over each orbit.  b: sqrt(4pi) on the Y00 row.
        std::vector<double> A(size_t(nY)*nOrb, 0.0), cost(nOrb);
        for (int k = 0; k < nOrb; ++k)
        {
            cost[k] = double(orbits[k].size());
            for (const auto& m : orbits[k])
            {
                EvalRealYlm(m, L, Y);
                for (int j = 0; j < nY; ++j) A[size_t(j)*nOrb + k] += Y[j];
            }
        }
        const bool ok = SolveNNLS(A, nY, nOrb, cost, w, 1e-9);
        if (getenv("QCHEM_ANGMESH_DEBUG"))
        {
            int    sup = 0; size_t nd = 0;
            for (int k = 0; k < nOrb; ++k) if (w[k] > 0.0) {++sup; nd += orbits[k].size();}
            printf("[angmesh] pool=%d (special %d) nnls=%s support=%d dirs=%zu\n",
                   nOrb, nSpecial, ok ? "ok" : "FAIL", sup, nd);
        }
        if (ok)
        {
            size_t n = 0;                                 // flatten: every member carries its orbit weight
            for (int k = 0; k < nOrb; ++k) if (w[k] > 0.0) n += orbits[k].size();
            rvec3vec_t P(n);
            rvec_t     W(n);
            size_t i = 0;
            for (int k = 0; k < nOrb; ++k)
            {
                if (w[k] <= 0.0) continue;
                for (const auto& m : orbits[k]) {P[i] = m; W[i] = w[k]; ++i;}
            }
            return qcMesh::Mesh(std::move(P), std::move(W));
        }
        if (nOrb < K) break;                              // seed pool exhausted (assert below)
    }
    assert(false && "MakeInvariantAngularMesh: no positive exact rule found (seed pool exhausted)");
    return qcMesh::Mesh();
}

} // namespace
