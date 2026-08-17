// File: BasisSet/G_FieldEvaluator.C  The plane-wave grid-engine DIP seams, segregated by client (V1.5).
//
// The dependency-inversion (SOLID DIP) seams between the structure-neutral fitting / Hamiltonian /
// ChargeDensity layers and a plane-wave {G} fit basis: those layers must not know the reciprocal lattice B,
// so they depend on these abstractions and ask; the concrete plane-wave basis (which owns B via its grid
// engine) answers.  A consumer reaches its face by an "I want more" cross-cast of the fit face it already
// holds -- never a cast into a concrete PW_Evaluator.
//
// V1.5 (2026-08-16): ONE interface became FOUR, cut along what each CLIENT consumes (the V1.27 lesson --
// name the capability for the consumption, not the mechanism).  The old union forced every consumer to see
// eleven methods when it used one to three, and forced the seed's structure-factor assembly and the mixer's
// Kerker filter to ride the quadrature engine:
//   - G_FieldEvaluator   : evaluate a fitted coefficient map as a real field (the ortho fitters' op(r)).
//   - G_Quadrature       : the FFT quadrature grid engine (fit sampling, XC energy, rho on grid).
//   - G_StructureFactor  : the analytic species-superposition density (the SAD seed's rho-tilde).
//   - G_SpectralFilter   : an isotropic spectral multiplier on the raster (the mixer's Kerker step).
// PW_Grid_Evaluator implements all four; each consumer casts to exactly the face it consumes.
module;
#include <functional>
export module qchem.BasisSet.G_FieldEvaluator;
import qchem.BasisSet.Internal.GMap;   // ΔG_Map (the G-space coefficient map to evaluate)
import qchem.Types;       // rvec3_t, rvec_t, cvec_t, rvec3vec_t, ivec3_t, dcmplx
import qchem.Structure;   // Structure, Atom (MakeFourierDensity's structure-factor sum)

export namespace qchem::BasisSet
{

//! \brief Evaluate a Hermitian G-space coefficient map \f$c(\Delta m)\f$ as a REAL field -- the ortho
//! fitters' \c op(r)/Gradient (what the GUI plots; the seed of the \f$v_{xc}-v_{xc,fit}\f$ fit-residual
//! diagnostic).  Needs \f$B\f$ (\f$\Delta m\to\f$ Cartesian \f$G\f$), which is exactly why it is a seam.
class G_FieldEvaluator
{
public:
    virtual ~G_FieldEvaluator() = default;
    //! \f$f(\vec r)=\mathrm{Re}\sum_{\Delta m} c(\Delta m)\,e^{i(B\Delta m)\cdot\vec r}\f$.
    virtual double  EvalField        (const ΔG_Map& c, const rvec3_t& r) const=0;
    //! \f$\nabla f(\vec r)=\sum_{\Delta m}(B\Delta m)\,(-\mathrm{Im}[c(\Delta m)e^{i(B\Delta m)\cdot\vec r}])\f$.
    virtual rvec3_t EvalFieldGradient(const ΔG_Map& c, const rvec3_t& r) const=0;
};

//! \brief The FFT QUADRATURE grid engine: sample points, forward/inverse transforms, coefficient lookup,
//! and the grid integral.  A scalar fitter samples \f$v_{xc}(r)\f$ on THIS grid (the FIT basis's own,
//! possibly denser than the orbital's -- relCutoff / GridCutoffFactor control it) and the XC term
//! integrates \f$E_{xc}\f$ on the same grid, borrowed through the fitter (\c GriddedScalarFitter::Grid(),
//! one owner).  Pure density/potential machinery: nothing here assumes plane-wave ORBITALS, so a
//! Gaussian-orbital (GPW) density grid implements it wholesale.
class G_Quadrature
{
public:
    virtual ~G_Quadrature() = default;
    //! Cartesian points of the FFT grid (raster order), the quadrature mesh a scalar fitter samples a field on.
    virtual const rvec3vec_t& GridPoints() const=0;
    //! Inverse-FFT a G-space coefficient map (keyed by reciprocal-index difference) to \f$\rho(r)\f$ on the grid.
    virtual rvec_t     RhoOnGrid  (const ΔG_Map& rhoTilde) const=0;
    //! Forward-FFT a real-space grid field to the FULL normalised (\f$/N_{pts}\f$) G-space grid (raster order).
    virtual cvec_t     ForwardFFT (const rvec_t& V) const=0;
    //! Look up \f$\tilde V(\Delta m)\f$ in a \c ForwardFFT grid (wrapping \f$\Delta m\f$ into THIS grid).
    virtual dcmplx     GridCoeff  (const cvec_t& Vt, const ivec3_t& dm) const=0;
    //! The evaluatable fitted-field coefficients \f$c(G)=\tilde V(G)\f$ over THIS basis's own \f$\{G\}\f$ (a
    //! \c GridCoeff gather).  \c G_FieldEvaluator::EvalField then plots \f$\sum_G c(G)e^{iG\cdot r}\f$.
    virtual ΔG_Map     FieldCoeffs(const cvec_t& Vt) const=0;
    //! \f$\int f\,d^3r\f$ on the FFT grid (weight \f$\Omega/N_{pts}\f$) -- the XC energy quadrature on the fit grid.
    virtual double     Integral(const rvec_t& f) const=0;
    // NB: grid REPORTING is deliberately not here: the owning fit basis self-reports at construction,
    // role-labeled (see PlaneWaveFit_IBS; user ruling 2026-08-16).
};

//! \brief The analytic structure-factor DENSITY over this engine's own \f$\{G\}\f$:
//! \f$\tilde\rho(G)=\frac1\Omega\sum_{\text{atoms}} f(Z,|G|^2)\,e^{-iG\cdot R}\f$, with \a f the per-species
//! radial form factor (an atomic density's 1-D Fourier transform).  The ANALYTIC assembly (no 3-D grid, so
//! no aliasing of a peaked density) that builds a SAD seed's \f$\tilde\rho\f$; the seed reaches it through
//! its OWN density-fit basis, never the orbital basis.  Keeps \f$G=0\f$ (the total charge).  The density
//! analogue of the pseudopotential's MakeLocalPotential.
class G_StructureFactor
{
public:
    virtual ~G_StructureFactor() = default;
    virtual ΔG_Map MakeFourierDensity(const Structure* atoms,
                          const std::function<double(int Z, double g2)>& formFactor) const=0;
};

//! \brief Apply an ISOTROPIC spectral multiplier to a real grid field over the FULL FFT box:
//! \f$f\mapsto\mathcal F^{-1}[k(|G|^2)\,\mathcal F f]\f$.  A SMOOTH \a k truncates nothing, so no Gibbs
//! ringing is introduced -- the raster-space Kerker preconditioner \f$k=G^2/(G^2+G_0^2)\f$ of the raw-XC
//! mixing pipeline (doc/GPWPlan 0.5(f2)); \f$k(0)\f$ scales the mean (Kerker's \f$k(0)=0\f$ conserves
//! charge).  The density MIXER's face -- its one and only ask of the grid engine.
class G_SpectralFilter
{
public:
    virtual ~G_SpectralFilter() = default;
    virtual rvec_t ApplySpectralFilter(const rvec_t& f,
                                       const std::function<double(double g2)>& k) const=0;
};

} //namespace
