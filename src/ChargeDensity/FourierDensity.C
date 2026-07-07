// File: ChargeDensity/FourierDensity.C  A density's reciprocal-space (G-space) coefficients.
//
// For a PERIODIC system rho-tilde(G) is the density's NATIVE representation -- you diagonalize in
// G-space and store rho-tilde; the real-space rho(r) is the FFT-derived view, not the reverse.  So this
// is core charge-density functionality for solids, not an escape hatch: it is the dual of
// BasisSet::Band_FT_IBS, and lets the plane-wave basis assemble the Hartree (FFT-free) and XC (FFT)
// matrices from rho-tilde directly instead of O(Npts*n^2) pointwise sampling.  A FINITE (molecular)
// density has no reciprocal-lattice Fourier series, so it does not provide one (cf. Structure::isFinite);
// a DFT term reaches a periodic density's rho-tilde by dynamic_cast (abstract->abstract).  The composite
// density sums the per-block (BZ-weighted) rho-tilde over the k-mesh.
module;
#include <type_traits>
export module qchem.ChargeDensity.FourierDensity;
export import qchem.BasisSet.Internal.GMap;
import qchem.BasisSet.Fit_IBS;   // cFIT_CD_ABS (the CD fit basis GetRepulsion3C keys by)
import qchem.Types;   // dcmplx

export namespace qchem::ChargeDensity
{

class FourierDensity
{
public:
    virtual ~FourierDensity() {}
    //! \brief The density's METRIC-FREE reciprocal-space coefficients \f$\tilde\rho(\Delta m)\f$ (BZ-weighted;
    //! the composite sums \f$\sum_k w_k\tilde\rho_k\f$).  Feeds the XC \f$\rho(r)\f$ path (inverse FFT).
    virtual ΔG_Map GetFourierDensity() const=0;

    //! \brief The density's COULOMB projection \f$V_H(\Delta m)=4\pi\tilde\rho/|G|^2\f$ for CD fit basis \a c
    //! -- the reciprocal analogue of the molecular \c IrrepCD::GetRepulsion3C.  A matrix-carrying density
    //! contracts \f$D\f$ against \c Band_FT_IBS::Repulsion3C (kernel baked); a matrix-free seed applies
    //! \c CoulombKernel to its \f$\tilde\rho\f$.  The Hartree term assembles \f$\langle i|V_H|j\rangle\f$ from it.
    virtual ΔG_Map GetRepulsion3C(const BasisSet::cFIT_CD_ABS& c) const=0;
};

//! Empty (non-polymorphic) stand-in for a FINITE density, which has no reciprocal-space representation.
//! A density template inherits FourierDensity only on the periodic (dcmplx) path; the finite (double)
//! path gets this empty base so its object layout is UNCHANGED -- a polymorphic virtual base would shift
//! the object size/allocation and perturb the (basin-sensitive) molecular SCF.
struct NoFourierDensity {};

//! FourierDensity for the periodic path (T=dcmplx), the empty base for the finite path.
template <class T> using FourierDensityBase =
    std::conditional_t<std::is_same_v<T,dcmplx>, FourierDensity, NoFourierDensity>;

} //namespace
