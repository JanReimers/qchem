// File: Fitting/FitOperations.C  A free ALGORITHM over the scalar fit-basis face.
//
// OrthogonalFit is not an interface, so it does not belong in Fit_IBS.C (user, 2026-08-23) -- it is a
// generic operation WRITTEN IN TERMS OF one.  And since 2026-08-24 it is not in qcBasisSet either (user:
// "OrthogonalFit belongs in the qcFitting library.  Client code needs to use that framework, not dodge
// around it"): PERFORMING A FIT is what this library is for, and a fit algorithm sitting in the BASIS
// library was an invitation to use it instead of a fitter.
//
// Kept as a free template rather than a member precisely BECAUSE it is not part of any face's contract: a
// fit basis promises the projection, the metric diagonal and its 3-centre overlap; how a caller combines
// those is the caller's algorithm, and this library owns that algorithm.
//
// ⚠ WHAT THE MOVE DOES NOT YET FIX, and it is the second half of the user's sentence.  Two callers still
// use this INSTEAD of a FunctionFitter_Scalar -- XC_SinglesQuadrature's matrix-free branches and
// tDM_CD::ProjectOnto's default -- because both want the COEFFICIENT VECTOR and the fitter face
// deliberately has no accessor for one (a Coefficients() getter is the smell deleted in increment 6).
// Being in this library at least makes that visible for what it is: framework algorithm, no framework
// object.  The remaining question is not how to hide the coefficients -- v_xc is pointwise NONLINEAR, so
// for a delta basis they genuinely must reach the term as rho(r_g) -- but WHAT TYPE should carry them.
// See doc/OpenWork.md, item 2.
//
// GONE 2026-08-24 (user): Overlap3CFace<U>(b), which returned the fit basis's 3-centre face for a given
// block scalar.  Once the face itself took BOTH axes (Integrals_Overlap3C<U,TFit>, increment 7) the helper
// was one reference cross-cast wrapped in an if constexpr, saving two lines at three call sites -- and the
// if constexpr arm was buying nothing: the same-scalar case is a derived-to-virtual-base cast, which never
// fails, and every call site runs once per Fock block, not in a loop.  Each site now spells its own cast.
module;
#include <cassert>
export module qchem.Fitting.FitOperations;
export import qchem.BasisSet.Orbital_DFT_IBS;   // FIT_SF_ABS<T> -- the face this operates on
import qchem.Blaze;                     // blazem::real -- these are TEMPLATES, so they need it HERE

export namespace qchem::Fitting
{

//! \brief THE ORTHOGONAL FIT: \f$c_a=\langle f_a|f\rangle/\langle f_a|f_a\rangle\f$ -- the projection
//! divided by the metric diagonal, which is the whole of a fit when \c isOrtho() is true.
//!
//! REAL by return type, and honestly so: \a f is a real field and the representations that reach here have
//! a real metric diagonal (\f$\delta\f$: \f$w_g\f$) and hence a real projection.  Divided COMPONENT-WISE on
//! the real parts rather than as a complex quotient: \c std::complex division of \f$(x,0)/(y,0)\f$ goes
//! through \f$(xy)/(y^2)\f$, which is not \f$x/y\f$ to the last bit -- and for \f$\delta\f$ the whole point
//! is that \f$w_g f_g/w_g\f$ cancels.
template <class T> inline rvec_t OrthogonalFit(const BasisSet::FIT_SF_ABS<T>& b, const ScalarFunction<double>& f)
{
    assert(b.isOrtho() && "OrthogonalFit: the fit is projection/diagonal only for an ORTHOGONAL basis");
    const vec_t<T> p=b.Overlap(f), d=b.OverlapDiagonal();
    assert(p.size()==d.size() && "OrthogonalFit: one projection and one metric entry per fit function");
    rvec_t c(p.size());
    for (size_t a=0; a<p.size(); a++) c[a]=blazem::real(p[a])/blazem::real(d[a]);
    return c;
}


}//namespace
