// File: Symmetry/Imp/SphericalRep.C  Implementation of the real-spherical shell representation.
module;
#include <vector>
#include <array>
#include <utility>
#include <map>
#include <cmath>
#include <cassert>
module qchem.Symmetry.SphericalRep;

namespace qchem::Symmetry
{

// Tiny dense helpers (shells are <= 7 wide); hand-rolled so this stays independent of which free-function
// Blaze ops the module surface re-exports.
static rmat_t transpose(const rmat_t& A)
{
    rmat_t T(A.columns(), A.rows(), 0.0);
    for (size_t i = 0; i < A.rows(); ++i) for (size_t j = 0; j < A.columns(); ++j) T(j, i) = A(i, j);
    return T;
}
static rmat_t matmul(const rmat_t& A, const rmat_t& B)
{
    assert(A.columns() == B.rows());
    rmat_t C(A.rows(), B.columns(), 0.0);
    for (size_t i = 0; i < A.rows(); ++i)
        for (size_t p = 0; p < A.columns(); ++p)
        {
            const double a = A(i, p); if (a == 0.0) continue;
            for (size_t j = 0; j < B.columns(); ++j) C(i, j) += a * B(p, j);
        }
    return C;
}

// Small dense inverse via Gauss-Jordan with partial pivoting.  qcSymmetry is LAPACK-free; the matrix here
// is the (2l+1)x(2l+1) Gram matrix C C^T (SPD, size <= 7 for s..f), so a plain inverse is fine and cheap.
static rmat_t SmallInverse(rmat_t A)
{
    const size_t n = A.rows();
    rmat_t I(n, n, 0.0);
    for (size_t i = 0; i < n; ++i) I(i, i) = 1.0;
    for (size_t col = 0; col < n; ++col)
    {
        size_t piv = col; double best = std::abs(A(col, col));
        for (size_t r = col + 1; r < n; ++r)
            if (std::abs(A(r, col)) > best) { best = std::abs(A(r, col)); piv = r; }
        assert(best > 1e-300 && "SmallInverse: singular Gram matrix (harmonics not independent?)");
        if (piv != col)
            for (size_t c = 0; c < n; ++c) { std::swap(A(col, c), A(piv, c)); std::swap(I(col, c), I(piv, c)); }
        const double d = A(col, col);
        for (size_t c = 0; c < n; ++c) { A(col, c) /= d; I(col, c) /= d; }
        for (size_t r = 0; r < n; ++r) if (r != col)
        {
            const double f = A(r, col);
            for (size_t c = 0; c < n; ++c) { A(r, c) -= f * A(col, c); I(r, c) -= f * I(col, c); }
        }
    }
    return I;
}

rmat_t SphericalShellRep::Rep(const Matrix3D<double>& R) const
{
    const HarmonicC2S& c2s = itsC2S;
    const size_t nSph = c2s.size();

    // Collect the distinct Cartesian monomials appearing across the harmonics -> the Cartesian shell.
    std::map<IVec3, size_t> index;
    std::vector<IVec3>      exps;
    for (const auto& h : c2s)
        for (const auto& [e, c] : h)
            if (index.find(e) == index.end()) { index[e] = exps.size(); exps.push_back(e); }
    const size_t nCart = exps.size();

    rmat_t C(nSph, nCart, 0.0);                      // harmonics in the monomial basis
    for (size_t m = 0; m < nSph; ++m)
        for (const auto& [e, c] : c2s[m]) C(m, index[e]) += c;

    const rmat_t Dcart = CartesianShellRep(exps).Rep(R); // nCart x nCart, the Cartesian rep on `exps`
    const rmat_t Ct    = transpose(C);
    const rmat_t Ginv  = SmallInverse(matmul(C, Ct));            // (C C^T)^{-1}, nSph x nSph
    return matmul(Ginv, matmul(C, matmul(Dcart, Ct)));          // D_sph = (CC^T)^{-1} C Dcart C^T
}

} //namespace
