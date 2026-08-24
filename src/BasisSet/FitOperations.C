// File: BasisSet/FitOperations.C  Free ALGORITHMS over the scalar fit-basis faces.
//
// Neither of these is an interface, so neither belongs in Fit_IBS.C (user, 2026-08-23) -- they are generic
// operations WRITTEN IN TERMS OF one, and each exists because two consumers in DIFFERENT libraries need
// the same few lines and neither owns the other:
//
//   OrthogonalFit   -- the scalar fitter (qcFitting) and a matrix-free density (qcChargeDensity)
//   Overlap3CFace   -- the scalar fitter (adjoint direction) and a density block (forward direction)
//
// Kept as free templates rather than members precisely BECAUSE they are not part of any face's contract: a
// fit basis promises the projection, the metric diagonal and its 3-centre overlap; how a caller combines
// those is the caller's algorithm, shared here only to have one copy of it.
module;
#include <cassert>
#include <type_traits>   // std::is_same_v -- Overlap3CFace's same-scalar/mixed arm
export module qchem.BasisSet.FitOperations;
export import qchem.BasisSet.Orbital_DFT_IBS;   // FIT_SF_ABS<T> / Integrals_Overlap3C<U> -- what these operate on
import qchem.Blaze;                     // blazem::real -- these are TEMPLATES, so they need it HERE

export namespace qchem::BasisSet
{

//! \brief THE ORTHOGONAL FIT: \f$c_a=\langle f_a|f\rangle/\langle f_a|f_a\rangle\f$ -- the projection
//! divided by the metric diagonal, which is the whole of a fit when \c isOrtho() is true.
//!
//! REAL by return type, and honestly so: \a f is a real field and the representations that reach here have
//! a real metric diagonal (\f$\delta\f$: \f$w_g\f$) and hence a real projection.  Divided COMPONENT-WISE on
//! the real parts rather than as a complex quotient: \c std::complex division of \f$(x,0)/(y,0)\f$ goes
//! through \f$(xy)/(y^2)\f$, which is not \f$x/y\f$ to the last bit -- and for \f$\delta\f$ the whole point
//! is that \f$w_g f_g/w_g\f$ cancels.
template <class T> inline rvec_t OrthogonalFit(const FIT_SF_ABS<T>& b, const ScalarFunction<double>& f)
{
    assert(b.isOrtho() && "OrthogonalFit: the fit is projection/diagonal only for an ORTHOGONAL basis");
    const vec_t<T> p=b.Overlap(f), d=b.OverlapDiagonal();
    assert(p.size()==d.size() && "OrthogonalFit: one projection and one metric entry per fit function");
    rvec_t c(p.size());
    for (size_t a=0; a<p.size(); a++) c[a]=blazem::real(p[a])/blazem::real(d[a]);
    return c;
}

//! \brief A fit basis's 3-centre face FOR A GIVEN ORBITAL-BLOCK SCALAR \a U.
//!
//! Same-scalar is a base and simply converts; the other one is the MIXED-run (3c-3) capability, so it is a
//! genuine "can you?" question -- a REFERENCE cross-cast, which throws rather than being release-mode UB.
//! \c if \c constexpr, not a run-time branch: which arm applies is fixed by the block's type.
template <class U, class TFit>
inline const Integrals_Overlap3C<U,TFit>& Overlap3CFace(const FIT_SF_ABS<TFit>& b)
{
    if constexpr (std::is_same_v<U,TFit>) return b;
    else return dynamic_cast<const Integrals_Overlap3C<U,TFit>&>(b);
}

}//namespace
