// File: BasisSet/G_FieldEvaluator.C  Real-space evaluation of a G-space coefficient map.
//
// The dependency-inversion (SOLID DIP) seam between an ORTHONORMAL (plane-wave) fitter and its {G} fit
// basis.  The ortho fitter's fitted field is a set of Hermitian G-space coefficients c(dm) (a ΔG_Map); to
// plot it -- rho_fit(r) or v_xc,fit(r) -- it must inverse-transform them:  f(r) = Re Σ_dm c(dm) e^{i(B·dm)·r}.
// That needs the reciprocal lattice B (dm -> Cartesian G), which the STRUCTURE-NEUTRAL fitting interface
// must not know.  So the fitter (high-level policy) depends on THIS abstraction and asks; the concrete
// plane-wave fit basis (which owns B via its grid engine) implements it.  Both depend on the abstraction,
// neither on the other's concretion -- and the fitter reaches it by an "I want more" cross-cast of the fit
// face it already holds, never a cast into a concrete PW_Evaluator.
module;
export module qchem.BasisSet.G_FieldEvaluator;
import qchem.Math.GMap;   // ΔG_Map (the G-space coefficient map to evaluate)
import qchem.Types;       // rvec3_t

export namespace qchem::BasisSet
{

//! \brief Evaluate a G-space coefficient map \f$c(\Delta m)\f$ (Hermitian) as a REAL field in real space,
//! \f$f(\vec r)=\mathrm{Re}\sum_{\Delta m} c(\Delta m)\,e^{i(B\,\Delta m)\cdot\vec r}\f$, and its gradient.
//! The capability an orthonormal (plane-wave) fitter cross-casts its held fit basis to, so it can plot its
//! fitted field without knowing the reciprocal lattice (DIP -- see the file note).  Implemented by the
//! concrete plane-wave fit basis (which owns \f$B\f$).
class G_FieldEvaluator
{
public:
    virtual ~G_FieldEvaluator() = default;
    //! \f$f(\vec r)=\mathrm{Re}\sum_{\Delta m} c(\Delta m)\,e^{i(B\Delta m)\cdot\vec r}\f$.
    virtual double  EvalField        (const ΔG_Map& c, const rvec3_t& r) const=0;
    //! \f$\nabla f(\vec r)=\sum_{\Delta m}(B\Delta m)\,(-\mathrm{Im}[c(\Delta m)e^{i(B\Delta m)\cdot\vec r}])\f$.
    virtual rvec3_t EvalFieldGradient(const ΔG_Map& c, const rvec3_t& r) const=0;
};

} //namespace
