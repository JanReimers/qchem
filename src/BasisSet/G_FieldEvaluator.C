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
#include <functional>
export module qchem.BasisSet.G_FieldEvaluator;
import qchem.BasisSet.Internal.GMap;   // ΔG_Map (the G-space coefficient map to evaluate)
import qchem.Types;       // rvec3_t, rvec_t, cvec_t, rvec3vec_t, ivec3_t, dcmplx
import qchem.Structure;   // Structure, Atom (MakeFourierDensity's structure-factor sum)

export namespace qchem::BasisSet
{

//! \brief The plane-wave GRID-ENGINE dependency-inversion (SOLID DIP) seam: the structure-neutral fitting
//! (and Hamiltonian XC) layer asks a plane-wave \f$\{G\}\f$ basis to do the reciprocal-lattice-aware FFT
//! quadrature it CANNOT do itself (it must not know \f$B\f$), and the concrete plane-wave basis -- which owns
//! the grid engine -- answers.  Two responsibilities, both needing \f$B\f$:
//!  - **evaluate** a Hermitian coefficient map \f$c(\Delta m)\f$ as a real field \f$f(r)=\mathrm{Re}\sum c\,e^{iG\cdot r}\f$ (GUI / fit-residual);
//!  - **quadrature** a field on the FFT grid: sample points, inverse/forward transforms, coefficient lookup.
//! An orthonormal scalar fitter reaches these by an "I want more" cross-cast of its OWN fit basis (the fit face
//! it already holds) -- never a cast into a concrete \c PW_Evaluator.  This is now a PURE density/potential grid
//! engine: the \f$\langle G_i|V|G_j\rangle=\tilde V(m_i-m_j)\f$ potential->orbital-matrix assembly (the ONE
//! method that assumed the ORBITALS are plane waves) has moved OFF here onto the orbital face
//! (\c Band_FT_IBS::MakeOverlap), so a Gaussian-orbital (GPW) density grid can reuse this engine wholesale.
//! Implemented by \c PW_Evaluator, so both the orbital and the auxiliary fit basis carry it.
class G_FieldEvaluator
{
public:
    virtual ~G_FieldEvaluator() = default;

    // --- evaluate a coefficient map as a real field (GUI / fit-residual diagnostic) ---
    //! \f$f(\vec r)=\mathrm{Re}\sum_{\Delta m} c(\Delta m)\,e^{i(B\Delta m)\cdot\vec r}\f$.
    virtual double  EvalField        (const ΔG_Map& c, const rvec3_t& r) const=0;
    //! \f$\nabla f(\vec r)=\sum_{\Delta m}(B\Delta m)\,(-\mathrm{Im}[c(\Delta m)e^{i(B\Delta m)\cdot\vec r}])\f$.
    virtual rvec3_t EvalFieldGradient(const ΔG_Map& c, const rvec3_t& r) const=0;

    // --- the FFT quadrature grid engine (a fit basis quadratures v_xc on its OWN, possibly denser, grid) ---
    //! Cartesian points of the FFT grid (raster order), the quadrature mesh a scalar fitter samples a field on.
    virtual const rvec3vec_t& GridPoints() const=0;
    //! Inverse-FFT a G-space coefficient map (keyed by reciprocal-index difference) to \f$\rho(r)\f$ on the grid.
    virtual rvec_t     RhoOnGrid  (const ΔG_Map& rhoTilde) const=0;
    //! Forward-FFT a real-space grid field to the FULL normalised (\f$/N_{pts}\f$) G-space grid (raster order).
    virtual cvec_t     ForwardFFT (const rvec_t& V) const=0;
    //! Look up \f$\tilde V(\Delta m)\f$ in a \c ForwardFFT grid (wrapping \f$\Delta m\f$ into THIS grid).
    virtual dcmplx     GridCoeff  (const cvec_t& Vt, const ivec3_t& dm) const=0;
    //! The evaluatable fitted-field coefficients \f$c(G)=\tilde V(G)\f$ over THIS basis's own \f$\{G\}\f$ (a
    //! \c GridCoeff gather).  \c EvalField plots \f$\sum_G c(G)e^{iG\cdot r}\f$ -- the scalar fitter's op(r).
    virtual ΔG_Map     FieldCoeffs(const cvec_t& Vt) const=0;
    //! \f$\int f\,d^3r\f$ on the FFT grid (weight \f$\Omega/N_{pts}\f$) -- the XC energy quadrature on the fit grid.
    virtual double     Integral(const rvec_t& f) const=0;

    //! \brief Analytic structure-factor DENSITY over THIS engine's own \f$\{G\}\f$:
    //! \f$\tilde\rho(G)=\frac1\Omega\sum_{\text{atoms}} f(Z,|G|^2)\,e^{-iG\cdot R}\f$, with \a f the per-species
    //! radial form factor (an atomic density's 1-D Fourier transform).  This is the ANALYTIC assembly (no 3-D
    //! grid, so no aliasing of a peaked density) that builds a SAD seed's \f$\tilde\rho\f$; the seed reaches it
    //! through its OWN density-fit basis (which is a \c G_FieldEvaluator), never the orbital basis.  Keeps
    //! \f$G=0\f$ (the total charge).  The density analogue of the pseudopotential's MakeLocalPotential.
    virtual ΔG_Map MakeFourierDensity(const Structure* atoms,
                          const std::function<double(int Z, double g2)>& formFactor) const=0;
};

} //namespace
