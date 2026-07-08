// File: BasisSet/Lattice_3D/Evaluators/PeriodicGridEvaluator.C  The reciprocal-cell FFT/Poisson grid engine.
//
// The k-INDEPENDENT density/potential grid machinery shared by every periodic evaluator (plane waves today,
// GPW next).  Given the reciprocal cell B, the direct-cell volume Omega, and an FFT grid resolution N, it
// answers the pure {r} <-> {G} FFT-quadrature questions -- sample points, forward/inverse transforms,
// coefficient lookup, real-space integration, and evaluating a G-space coefficient map as a real field.
//
// It knows NOTHING about the orbitals: no {G} cutoff set, no k, no Bloch phase.  A concrete evaluator
// (PW_Evaluator, later GPW_Evaluator) HOLDS one of these and forwards its G_FieldEvaluator grid virtuals to
// it, while supplying the orbital-specific part (Eval, overlap/kinetic, the potential->matrix bridge) itself
// -- the same "share the data + code, add only the functional form" split as the atom side's
// ExponentialEvaluator under Slater/Gaussian.  Because the engine is k-independent, ONE instance can be
// shared across the Brillouin-zone k-blocks of a calculation (the container hands each block a shared_ptr).
module;
#include <vector>
export module qchem.BasisSet.Lattice_3D.Evaluators.PeriodicGridEvaluator;
export import qchem.ReciprocalLattice;        // ReciprocalLattice / UnitCell (the B cell; source of G)
export import qchem.BasisSet.Internal.GMap;   // ΔG_Map (the G-space coefficient map EvalField/RhoOnGrid speak)
import qchem.Types;                           // ivec3_t, rvec3_t, rvec_t, rvec3vec_t, cvec_t, dcmplx

export namespace qchem::BasisSet::Lattice_3D
{

//! \brief The k-independent FFT/Poisson grid engine of a periodic basis: the reciprocal cell \f$B\f$, the
//! direct-cell volume \f$\Omega\f$, and the real-space FFT grid \f$N\f$, plus the pure \f$\{r\}\leftrightarrow
//! \{G\}\f$ quadrature it needs.  A concrete evaluator holds one and forwards the grid virtuals to it; the
//! \f$\{G\}\f$-cutoff set, the crystal momentum \f$k\f$, and the orbital matrices are NOT here (they are the
//! per-k, orbital-specific part on the concrete evaluator).  Held by \c shared_ptr so one engine serves many
//! k-blocks.
class PeriodicGridEvaluator
{
public:
    //! \a recip = reciprocal cell \f$B\f$; \a volume = direct-cell volume \f$\Omega\f$; \a fftGrid = the FFT
    //! divisions \f$N\f$ (a concrete evaluator sizes \f$N\f$ from its own \f$\{G\}\f$ extent and hands it in).
    PeriodicGridEvaluator(const ReciprocalLattice& recip, double volume, const ivec3_t& fftGrid)
        : itsRecip(recip), itsVolume(volume), itsN(fftGrid) {}

    double                   Volume()  const {return itsVolume;}   //!< direct cell volume \f$\Omega\f$
    const ReciprocalLattice& Recip()   const {return itsRecip;}    //!< reciprocal cell (matrix \f$B\f$)
    const ivec3_t&           FFTGrid() const {return itsN;}        //!< FFT grid divisions \f$N\f$
    rvec3_t GetGCartesian(const ivec3_t& m) const;                 //!< \f$G=B\,m\f$ (Cartesian a.u.)

    //! Fractional \f$(i/n)\f$ grid (weight \f$\Omega/\prod n\f$), the raster order \c RhoOnGrid / \c ForwardFFT use.
    std::vector<rvec3_t> UniformGrid(const ivec3_t& n) const;
    //! Cartesian points of the FFT grid (raster order): \f$r=A\,(i/N)\f$ -- the quadrature mesh a fitter samples on.
    const rvec3vec_t&    GridPoints() const;
    //! \f$\rho(r)\f$ on the FFT grid = inverse-FFT of a G-space map keyed by \f$\Delta m\f$ (coefficients are physical, \f$/\Omega\f$).
    rvec_t   RhoOnGrid  (const ΔG_Map& rhoTilde) const;
    //! Forward-FFT a real-space grid field to the FULL, normalised (\f$/N_{pts}\f$) G-space grid (raster order).
    cvec_t   ForwardFFT (const rvec_t& V) const;
    //! Look up \f$\tilde V(\Delta m)\f$ in a \c ForwardFFT grid, wrapping \f$\Delta m\f$ into \f$N\f$ (alias-free for \f$|\Delta m|<N/2\f$).
    dcmplx   GridCoeff  (const cvec_t& Vt, const ivec3_t& dm) const;
    //! \f$\int f\,d^3r\f$ on the FFT grid: uniform quadrature, weight \f$\Omega/N_{pts}\f$.
    double   Integral   (const rvec_t& f) const;
    //! \f$f(\vec r)=\mathrm{Re}\sum_{\Delta m} c(\Delta m)\,e^{i(B\Delta m)\cdot\vec r}\f$ (a sparse point eval, no FFT).
    double   EvalField        (const ΔG_Map& c, const rvec3_t& r) const;
    //! \f$\nabla f(\vec r)=\sum_{\Delta m}(B\Delta m)\,(-\mathrm{Im}[c(\Delta m)e^{i(B\Delta m)\cdot\vec r}])\f$.
    rvec3_t  EvalFieldGradient(const ΔG_Map& c, const rvec3_t& r) const;

private:
    ReciprocalLattice  itsRecip;      //!< reciprocal cell \f$B\f$; source of \f$G=B\,m\f$ and the direct cell \f$A\f$
    double             itsVolume;     //!< direct cell volume \f$\Omega\f$
    ivec3_t            itsN;          //!< FFT grid divisions \f$N\f$
    mutable rvec3vec_t itsGridPoints; //!< cached Cartesian grid points (built once)
};

} //namespace
