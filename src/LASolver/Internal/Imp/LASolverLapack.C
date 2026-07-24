// File: LASolver/Internal/Imp/LASolverLapack.C  Eigen/SVD/Cholesky solver implementations.
module;
#include <cassert>
#include <iostream>
#include <type_traits>
#include <cmath>
#include <algorithm>
#include <memory>
module qchem.LASolver.Internal.Lapack;
import qchem.Blaze;

namespace qchem {

// One-line conditioning report of the basis overlap S (min/max eigenvalue, min singular value = min|eig| for
// Hermitian S, condition number), emitted at SetBasisOverlap when the process-wide toggle is on -- so a
// near-singular (small min eig) or indefinite (negative min eig -> Cholesky then fails) overlap is visible.
template <class T> static void report_conditioning(const hmat_t<T>& S)
{
    if (!ReportOverlapConditioning()) return;
    rvec_t d; mat_t<T> U;
    blazem::eigen(S, d, U);                          // ascending eigenvalues of the Hermitian overlap
    double mn = d[0], mx = d[d.size()-1], msv = std::fabs(d[0]);
    for (double v : d) msv = std::min(msv, std::fabs(v));
    std::cout << "[overlap S] n=" << S.rows() << "  min eig=" << mn << "  min sv=" << msv
              << "  max eig=" << mx << "  cond=" << (msv > 0 ? mx/msv : 0.0) << std::endl;
}

// Gap-detection auto-tolerance (doc/GPWPlan §1 step 2).  Given the ASCENDING overlap eigenvalues w, decide
// how many leading NEAR-NULL modes to drop -- without the user guessing a tolerance:
//   (1) force-drop w<=0 (indefinite from roundoff -- Cholesky would fail, canonical ortho must not sqrt them);
//   (2) in the LOW spectrum (below kLowFrac*max, where an over-complete null CLUSTER lives) find the largest
//       consecutive ratio w[i+1]/w[i]; a ratio > kGapRatio is a CLEAN spectral gap -> cut there;
//   (3) else (a smoothly over-complete basis: no clean gap, a geometric ramp) fall back to the ABSOLUTE
//       singularity floor kEpsFloor -- drop every mode below it.  Absolute, NOT wmax-relative: an
//       over-complete basis can have a huge wmax (measured 8.9e7), which would make a relative floor cut
//       arbitrarily high; and a normal basis with a 1e-9 near-null needs the floor to actually reach it.
// Returns {nDrop, gap ratio, clean?, the tolerance ACTUALLY APPLIED} (for a never-silent, self-explaining WARN).
struct NullGap { size_t nDrop; double gapRatio; bool clean; double tol; };
static NullGap detect_null_gap(const rvec_t& w)
{
    static constexpr double kGapRatio = 30.0;   //!< a >30x consecutive jump in the low spectrum = a clean gap
    static constexpr double kLowFrac  = 1e-3;   //!< only hunt the null/physical gap below max*kLowFrac
    static constexpr double kEpsFloor = 1e-8;   //!< ABSOLUTE singularity floor (no-clean-gap fallback)
    const size_t n = w.size();
    if (n==0) return {0, 1.0, true, 0.0};
    const double wmax = w[n-1];
    size_t k = 0;                                // (1) leading non-positive
    while (k<n && w[k] <= 0.0) ++k;
    size_t cut = k; double bestR = 1.0;          // (2) largest low-spectrum ratio
    for (size_t i=k; i+1<n && w[i] < wmax*kLowFrac; ++i)
    {
        const double r = w[i+1]/w[i];
        if (r > bestR) { bestR = r; cut = i+1; }
    }
    if (bestR > kGapRatio) return {cut, bestR, true, w[cut>0?cut-1:0]};     // clean gap: tol = last dropped eig
    size_t m = k;                                                          // (3) ambiguous -> absolute floor
    while (m<n && w[m] < kEpsFloor) ++m;
    return {m, bestR, m==0, kEpsFloor};
}

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
    report_conditioning<T>(S);
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
    const size_t n = w.size();
    size_t index = 0;
    bool ambiguous = false; double gapRatio = 0.0, appliedTol = tol;
    if (tol < 0.0)                                  // AUTO tolerance: gap detection (doc/GPWPlan §1 step 2)
    {
        const NullGap g = detect_null_gap(w);
        index = g.nDrop; gapRatio = g.gapRatio; ambiguous = (index>0 && !g.clean); appliedTol = g.tol;
    }
    else                                            // EXPLICIT floor (tol=0 keeps everything -- PSD w>0)
        for (auto v : w) { if (v < tol) index++; else break; }

    if (index == 0) return;
    // NEVER SILENT (doc/GPWPlan §1 step 3): every drop is reported on cerr with the evidence, INCLUDING the
    // tolerance actually applied (a clean gap's boundary eig, or the absolute floor for a smooth over-complete
    // ramp -- e.g. an even-tempered basis, whose overlap is a geometric ramp with no clean spectral gap).
    std::cerr << "[ortho] dropped " << index << " near-null overlap mode(s) of " << n
              << " (min eig=" << w[0] << ", first kept=" << (index<n ? w[index] : 0.0) << "; ";
    if (tol < 0.0 && !ambiguous) std::cerr << "AUTO clean gap: ratio " << gapRatio << " at eig " << appliedTol;
    else if (tol < 0.0)          std::cerr << "AUTO no clean gap (max ratio " << gapRatio
                                           << ") -> dropped below floor " << appliedTol;
    else                         std::cerr << "explicit tol " << tol;
    std::cerr << ")" << std::endl;

    assert(U.rows() == U.columns());
    const size_t nr = U.rows();
    U = blazem::submatrix(U, 0, index, nr, nr-index);
    w = blazem::subvector(w, index, nr-index);
}

//=============================================================================
//  LASolverSVD
//=============================================================================

template <class T> void LASolverSVD<T>::SetBasisOverlap(const hmat_t<T>& S)
{
    report_conditioning<T>(S);
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
    report_conditioning<T>(S);   // reported BEFORE potrf -- an indefinite S (min eig<0) is why Cholesky then fails
    mat_t<T> Sm(S);
    blazem::potrf(Sm, 'U');
    // Cholesky diagonal is always real; extract accordingly.
    if constexpr (std::is_floating_point_v<T>) itsD = blazem::diagonal(Sm);
    else                                       itsD = blazem::real(blazem::diagonal(Sm));
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
//  LASolverAuto -- eigen-analyse S, then pick Cholesky (well-conditioned) or Eigen(auto tol) (near-null).
//=============================================================================

template <class T> void LASolverAuto<T>::SetBasisOverlap(const hmat_t<T>& S)
{
    report_conditioning<T>(S);
    rvec_t w; mat_t<T> U;
    blazem::eigen(S, w, U);                        // ascending eigenvalues -- the "eigen analysis"
    const NullGap g = detect_null_gap(w);          // 0 droppable modes <=> S is safe for Cholesky
    if (g.nDrop == 0)
        itsImpl = std::make_unique<LASolverCholesky<T>>();
    else
    {
        std::cerr << "[ortho] Auto: overlap has " << g.nDrop << " near-null mode(s) (min eig=" << w[0]
                  << ") -> canonical Eigen orthogonalisation instead of Cholesky." << std::endl;
        itsImpl = std::make_unique<LASolverEigen<T>>(-1.0);   // <0 = AUTO gap-detection tol (re-WARNs the drop)
    }
    itsImpl->SetBasisOverlap(S);
}

//=============================================================================
//  Explicit instantiations
//=============================================================================

template class LASolverEigen  <double>;
template class LASolverSVD    <double>;
template class LASolverCholesky<double>;
template class LASolverAuto    <double>;
template class LASolverEigen  <dcmplx>;
template class LASolverSVD    <dcmplx>;
template class LASolverCholesky<dcmplx>;
template class LASolverAuto    <dcmplx>;

} // namespace qchem