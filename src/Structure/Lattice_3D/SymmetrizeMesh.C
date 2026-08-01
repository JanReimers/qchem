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
inline SL::Fold         FoldMesh       (const qcMesh::Mesh& m, const SL::SpaceGroup& sg, double tol=1e-8);
inline qcMesh::Mesh     Symmetrize     (const qcMesh::Mesh& m, const SL::SpaceGroup& sg, double tol=1e-8);
inline std::vector<int> UnmatchedCounts(const qcMesh::Mesh& m, const SL::SpaceGroup& sg, double tol=1e-8);
inline qcMesh::Mesh     MakeInvariant  (const qcMesh::Mesh& m, const SL::SpaceGroup& sg, double tol=1e-8);

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

} // namespace
