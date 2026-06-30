// File: LASolver/Internal/Imp/LASolverLapack.C  Eigen/SVD/Cholesky solver implementations.
module;
#include <cassert>
#include <iostream>
#include <type_traits>
module qchem.LASolver.Internal.Lapack;
import qchem.Blaze;

namespace qchem {

//-----------------------------------------------------------------------------
//  Internal helpers — not exported.  Template on V/Vd types so Blaze can use
//  typed (triangular, etc.) matrix BLAS paths for each solver variant.
//-----------------------------------------------------------------------------

template <class T>
static hmat_t<T> make_hermitian(mat_t<T>& A)
{
    if constexpr (std::is_floating_point_v<T>)
        A = 0.5*(A + blazem::trans(A));
    else
        A = 0.5*(A + blazem::ctrans(A));
    return hmat_t<T>(A);
}

template <class T, class Vd_t, class V_t>
static typename LASolver<T>::Ud_t
solve_impl(const hmat_t<T>& H, const Vd_t& Vd, const V_t& V)
{
    mat_t<T> Hp = Vd * H * V;
    hmat_t<T> Hs = make_hermitian<T>(Hp);
    rvec_t d; mat_t<T> U;
    blazem::eigen(Hs, d, U);
    U = V * U;
    return {U, d};
}

template <class T, class V_t>
static typename LASolver<T>::UUd_t
solve_ortho_impl(const hmat_t<T>& Hprime, const V_t& V)
{
    rvec_t d; mat_t<T> Up;
    blazem::eigen(Hprime, d, Up);
    mat_t<T> U = V * Up;
    return {U, Up, d};
}

template <class T, class Vd_t, class V_t>
static hmat_t<T>
transform_impl(const hmat_t<T>& M, const Vd_t& Vd, const V_t& V)
{
    mat_t<T> Mp = Vd * M * V;
    return make_hermitian<T>(Mp);
}

//=============================================================================
//  LASolverEigen
//=============================================================================

template <class T> void LASolverEigen<T>::SetBasisOverlap(const hmat_t<T>& S)
{
    rvec_t d; mat_t<T> U;
    blazem::eigen(S, d, U);
    Truncate(U, d, itsTruncationTolerance);
    Rescale(U, d);
    itsVd = blazem::ctrans(U);  // U^H: correct for both real and complex
    itsV  = std::move(U);
    itsD  = d;
}

template <class T> typename LASolver<T>::Ud_t
LASolverEigen<T>::Solve(const hmat_t<T>& H) const
{
    return solve_impl<T>(H, itsVd, itsV);
}

template <class T> typename LASolver<T>::UUd_t
LASolverEigen<T>::SolveOrtho(const hmat_t<T>& Hprime) const
{
    return solve_ortho_impl<T>(Hprime, itsV);
}

template <class T> hmat_t<T>
LASolverEigen<T>::Transform(const hmat_t<T>& M) const
{
    return transform_impl<T>(M, itsVd, itsV);
}

template <class T> mat_t<T>
LASolverEigen<T>::BackTransform(const mat_t<T>& Uprime) const
{
    return itsV * Uprime;
}

template <class T> void LASolverEigen<T>::Rescale(mat_t<T>& V, const rvec_t& w)
{
    for (size_t j = 0; j < V.columns(); j++)
        blazem::column(V, j) /= sqrt(w[j]);
}

template <class T> void LASolverEigen<T>::Truncate(mat_t<T>& U, rvec_t& w, double tol)
{
    assert(U.columns() == w.size());
    size_t index = 0;
    for (auto v : w)
        if (v < tol) index++;
        else         break;

    size_t n = w.size();
    assert(w[index] >= tol);
    if (index > 0)
    {
        assert(w[index-1] < tol);
        std::cout << "LASolverEigen truncating " << index << " eigen values."
                  << " Min(w)=" << w[0] << " tol=" << tol << std::endl;
        assert(U.rows() == U.columns());
        size_t nr = U.rows();
        U = blazem::submatrix(U, 0, index, nr, nr-index);
        w = blazem::subvector(w, index, nr-index);
    }
}

//=============================================================================
//  LASolverSVD
//=============================================================================

template <class T> void LASolverSVD<T>::SetBasisOverlap(const hmat_t<T>& S)
{
    rvec_t s; mat_t<T> U, Vt;
    blazem::svd(S, U, s, Vt);
    Truncate(U, s, Vt, itsTruncationTolerance);
    Rescale(U, s, Vt);
    itsV  = std::move(U);
    itsVd = std::move(Vt);
    itsD  = s;
}

template <class T> typename LASolver<T>::Ud_t
LASolverSVD<T>::Solve(const hmat_t<T>& H) const
{
    return solve_impl<T>(H, itsVd, itsV);
}

template <class T> typename LASolver<T>::UUd_t
LASolverSVD<T>::SolveOrtho(const hmat_t<T>& Hprime) const
{
    return solve_ortho_impl<T>(Hprime, itsV);
}

template <class T> hmat_t<T>
LASolverSVD<T>::Transform(const hmat_t<T>& M) const
{
    return transform_impl<T>(M, itsVd, itsV);
}

template <class T> mat_t<T>
LASolverSVD<T>::BackTransform(const mat_t<T>& Uprime) const
{
    return itsV * Uprime;
}

template <class T> void LASolverSVD<T>::Rescale(mat_t<T>& U, const rvec_t& s, mat_t<T>& Vt)
{
    for (size_t j = 0; j < U.columns(); j++)
        blazem::column(U, j) /= sqrt(s[j]);
    for (size_t i = 0; i < Vt.rows(); i++)
        blazem::row(Vt, i) /= sqrt(s[i]);
}

template <class T> void LASolverSVD<T>::Truncate(mat_t<T>& U, rvec_t& s, mat_t<T>& Vt, double tol)
{
    assert(blazem::isSquare(U));
    assert(blazem::isSquare(Vt));
    assert(U.rows()  == s.size());
    assert(Vt.rows() == s.size());
    size_t index = 0;
    for (auto v : s)
        if (v >= tol) index++;
        else          break;
    size_t n = s.size();
    assert(s[index-1] >= tol);
    assert(index == n || s[index] < tol);
    if (n - index > 0)
        std::cout << "LASolverSVD truncating " << n-index << " singular values."
                  << " Min(s)=" << s[n-1] << " tol=" << tol << std::endl;
    if (index < n)
    {
        U  = blazem::submatrix(U,  0, 0, n,     index);
        s  = blazem::subvector(s,  0,    index);
        Vt = blazem::submatrix(Vt, 0, 0, index, n    );
    }
}

//=============================================================================
//  LASolverCholesky
//=============================================================================

template <class T> void LASolverCholesky<T>::SetBasisOverlap(const hmat_t<T>& S)
{
    mat_t<T> Sm(S);
    blazem::potrf(Sm, 'U');
    // Cholesky diagonal is always real; extract accordingly.
    if constexpr (std::is_floating_point_v<T>)
        itsD = blazem::diagonal(Sm);
    else
        itsD = blazem::real(blazem::diagonal(Sm));
    blazem::trtri(Sm, 'U', 'N');
    // zero out lower triangle (trtri leaves it undefined)
    size_t N = Sm.rows();
    for (size_t i = 0; i < N; i++)
        for (size_t j = i+1; j < N; j++)
            Sm(j,i) = T(0);
    itsV  = umat_t(Sm);
    itsVd = lmat_t(blazem::ctrans(Sm));  // (U^{-1})^H for complex, same as trans for real
}

template <class T> typename LASolver<T>::Ud_t
LASolverCholesky<T>::Solve(const hmat_t<T>& H) const
{
    return solve_impl<T>(H, itsVd, itsV);
}

template <class T> typename LASolver<T>::UUd_t
LASolverCholesky<T>::SolveOrtho(const hmat_t<T>& Hprime) const
{
    return solve_ortho_impl<T>(Hprime, itsV);
}

template <class T> hmat_t<T>
LASolverCholesky<T>::Transform(const hmat_t<T>& M) const
{
    return transform_impl<T>(M, itsVd, itsV);
}

template <class T> mat_t<T>
LASolverCholesky<T>::BackTransform(const mat_t<T>& Uprime) const
{
    return itsV * Uprime;
}

//=============================================================================
//  Explicit instantiations
//=============================================================================

template class LASolverEigen  <double>;
template class LASolverSVD    <double>;
template class LASolverCholesky<double>;
template class LASolverEigen  <dcmplx>;
template class LASolverSVD    <dcmplx>;
template class LASolverCholesky<dcmplx>;

} // namespace qchem