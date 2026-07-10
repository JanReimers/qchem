// File: Symmetry/Lattice_3D/Imp/SpaceGroup.C  Space-group detection implementation.
module;
#include <vector>
#include <cassert>
module qchem.Symmetry.Lattice_3D.SpaceGroup;
import qchem.Math;   // sqrt, fabs, floor, lround (via qchem.CMath)

namespace qchem::Symmetry::Lattice_3D
{

//---------------------------------------------------------------------------------------
//  Small local helpers (file-internal).
//
static double dot(const rvec3_t& a, const rvec3_t& b) {return a.x*b.x + a.y*b.y + a.z*b.z;}

//! A 3x3 matrix from its three columns.
static Matrix3D<double> FromColumns(const rvec3_t& c1, const rvec3_t& c2, const rvec3_t& c3)
{
    return Matrix3D<double>(c1.x, c2.x, c3.x,
                            c1.y, c2.y, c3.y,
                            c1.z, c2.z, c3.z);
}

//! Reduce one coordinate into [0,1), snapping values within \a tol of an integer to 0.
static double Frac1(double x, double tol)
{
    double r = x - floor(x);
    if (r > 1.0 - tol) r -= 1.0;   // snap the just-below-1 case down to 0
    if (r < 0.0) r += 1.0;
    return r;
}

//! Do two fractional positions coincide modulo a lattice translation, within \a tol?
static bool SameSiteModLattice(const rvec3_t& a, const rvec3_t& b, double tol)
{
    for (int c = 0; c < 3; ++c)
    {
        double d = ((&a.x)[c]) - ((&b.x)[c]);
        d -= floor(d + 0.5);            // nearest-integer wrap into (-1/2, 1/2]
        if (fabs(d) > tol) return false;
    }
    return true;
}

//! Are two matrices equal to within \a tol (entrywise)?
static bool MatEqual(const Matrix3D<double>& a, const Matrix3D<double>& b, double tol)
{
    for (int i = 1; i <= 3; ++i)
        for (int j = 1; j <= 3; ++j)
            if (fabs(a(i,j) - b(i,j)) > tol) return false;
    return true;
}

//---------------------------------------------------------------------------------------
//  Holohedry: all integer ops W (lattice coordinates) with W^T M W = M, found by searching
//  for images of the three lattice vectors that preserve every pairwise dot product.  The
//  image candidates are enumerated lattice vectors of the matching length.
//
static std::vector<Matrix3D<double>> Holohedry(const Matrix3D<double>& A, double tol)
{
    const Matrix3D<double> M = Transpose(A) * A;   // metric: M(i,j) = a_i . a_j

    // Enumerate lattice vectors A*n for small integer n, tagged with their integer coords.
    // nmax=2 comfortably covers the same-length images for standard (cubic/FCC/BCC/hex) cells,
    // whose symmetry images have coordinates in {-1,0,1}.
    struct LVec { rvec3_t n; rvec3_t v; double len2; };
    std::vector<LVec> lattice;
    const int nmax = 2;
    for (int i = -nmax; i <= nmax; ++i)
        for (int j = -nmax; j <= nmax; ++j)
            for (int k = -nmax; k <= nmax; ++k)
            {
                rvec3_t n{double(i), double(j), double(k)};
                rvec3_t v = A * n;
                lattice.push_back({n, v, dot(v,v)});
            }

    // Length tolerance on |v|^2 (scaled so tol is a fractional length tolerance).
    const double cellScale = M(1,1) + M(2,2) + M(3,3);
    const double len2tol = 2.0 * tol * cellScale + tol * tol;

    // Candidate images for each of the three lattice vectors (same length).
    std::vector<rvec3_t> cand[3];
    for (int col = 0; col < 3; ++col)
        for (const auto& lv : lattice)
            if (fabs(lv.len2 - M(col+1,col+1)) <= len2tol)
                cand[col].push_back(lv.n);

    std::vector<Matrix3D<double>> ops;
    for (const auto& n1 : cand[0])
    for (const auto& n2 : cand[1])
    for (const auto& n3 : cand[2])
    {
        Matrix3D<double> W = FromColumns(n1, n2, n3);   // columns = integer images
        // Metric preservation  W^T M W == M  <=>  (A n_p).(A n_q) == M(p,q) for all p,q.
        if (MatEqual(Transpose(W) * M * W, M, len2tol))
            ops.push_back(W);
    }
    return ops;
}

//---------------------------------------------------------------------------------------
//  Space-group compatibility: does holohedry op W admit a fractional translation tau so
//  that {W|tau} maps the basis onto itself?  Trial tau values come from mapping atom 0 onto
//  each same-species atom; the first that works is returned (primitive-cell assumption).
//
static bool FindTau(const Matrix3D<double>& W, const std::vector<AtomSite>& basis,
                    double tol, rvec3_t& tauOut)
{
    if (basis.empty()) { tauOut = rvec3_t(0,0,0); return true; }

    const AtomSite& a0 = basis.front();
    const rvec3_t   Wf0 = W * a0.f;

    for (const AtomSite& aj : basis)
    {
        if (aj.species != a0.species) continue;
        rvec3_t tau(Frac1(aj.f.x - Wf0.x, tol),
                    Frac1(aj.f.y - Wf0.y, tol),
                    Frac1(aj.f.z - Wf0.z, tol));

        bool ok = true;
        for (const AtomSite& ak : basis)
        {
            rvec3_t g = W * ak.f + tau;
            bool matched = false;
            for (const AtomSite& am : basis)
                if (am.species == ak.species && SameSiteModLattice(g, am.f, tol))
                    { matched = true; break; }
            if (!matched) { ok = false; break; }
        }
        if (ok) { tauOut = tau; return true; }
    }
    return false;
}

//---------------------------------------------------------------------------------------
SpaceGroup SpaceGroup::Detect(const Matrix3D<double>& A,
                              const std::vector<AtomSite>& basis, double tol)
{
    std::vector<SpaceGroupOp> ops;
    for (const Matrix3D<double>& W : Holohedry(A, tol))
    {
        rvec3_t tau;
        if (FindTau(W, basis, tol, tau))
            ops.push_back({W, tau});
    }
    assert(!ops.empty());   // the identity is always a symmetry
    return SpaceGroup(A, std::move(ops));
}

bool SpaceGroup::isSymmorphic() const
{
    for (const auto& op : itsOps)
        if (fabs(op.tau.x) > 1e-9 || fabs(op.tau.y) > 1e-9 || fabs(op.tau.z) > 1e-9)
            return false;
    return true;
}

std::vector<Matrix3D<double>> SpaceGroup::PointGroupOps() const
{
    std::vector<Matrix3D<double>> W;
    W.reserve(itsOps.size());
    for (const auto& op : itsOps) W.push_back(op.W);
    return W;
}

std::vector<Matrix3D<double>> SpaceGroup::ReciprocalPointOps(bool includeTimeReversal) const
{
    // In the reciprocal fractional basis, direct-space op W acts on k as U = (W^{-1})^T.
    std::vector<Matrix3D<double>> U;
    auto add = [&](const Matrix3D<double>& u)
    {
        for (const auto& e : U) if (MatEqual(e, u, 1e-9)) return;   // dedup
        U.push_back(u);
    };
    for (const auto& op : itsOps)
    {
        Matrix3D<double> u = Transpose(Invert(op.W));
        add(u);
        if (includeTimeReversal) add(-1.0 * u);   // time reversal k -> -k
    }
    return U;
}

} // namespace
