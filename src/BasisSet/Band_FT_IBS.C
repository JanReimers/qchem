// File: BasisSet/Band_FT_IBS.C  Abstract reciprocal-space (Fourier / plane-wave) DFT-assembly capability.
//
// The dual of Band_DFT_IBS: instead of integrating real-space ScalarFunctions on a mesh, this is the
// RECIPROCAL-SPACE lineage -- a plane-wave basis and a plane-wave density share G-space coefficients (a
// FourierMap, rho-tilde / V-tilde).  The questions here are the DENSITY-driven KS assembly: rho-tilde ->
// the FFT-free Hartree matrix, and the FFT XC route -- replacing the O(Npts*n^2) pointwise density
// sampling with one O(n^2) projection D->rho-tilde.  (The EXTERNAL pseudopotential assembly is a separate
// capability, Pseudopotential::Integrals_Pseudo in qcPseudopotential, so qcBasisSet names no pseudopotential type.)
// A term reaches this the sanctioned way: holding the abstract orbital basis and dynamic_cast-ing UP.
module;
#include <functional>
export module qchem.BasisSet.Band_FT_IBS;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.FourierMap;
import qchem.Types;       // hmat_t<dcmplx>

export namespace BasisSet
{

//! \brief A basis that assembles matrices in reciprocal (G-)space: the density-driven KS matrices from a
//! density's Fourier coefficients, and the external (pseudo)potential matrices from a per-species form
//! factor.  Only the complex (plane-wave) path realizes this; a density exposes its rho-tilde via the
//! matching ChargeDensity FourierDensity sub-interface, and the BZ sum is the per-block rho-tilde summed.
class Band_FT_IBS
    : public virtual Orbital_1E_IBS<dcmplx>
{
public:
    using Orbital_1E_IBS<dcmplx>::Overlap;   // keep the cached no-arg Overlap() visible beside Overlap(rvec_t)

    //! \brief \f$\tilde\rho(\Delta m)=\frac1\Omega\sum_{G_i-G_j=\Delta m}D_{ij}\f$ for a density matrix
    //! \a D in THIS plane-wave block (one \f$O(n^2)\f$ accumulation over the difference set).
    virtual FourierMap MakeFourierDensity(const hmat_t<dcmplx>& D) const=0;

    //! \brief \f$\tilde\rho(\Delta m)=\frac1\Omega\sum_{\text{atoms}} f(Z,|B\Delta m|^2)\,e^{-iG\cdot R}\f$ over
    //! the difference set -- the structure-factor assembly of a per-species radial form factor \a formFactor
    //! (e.g. an atomic valence density for a SAD seed).  The density analogue of the pseudopotential's
    //! MakeLocalPotential, but it KEEPS \f$\Delta m=0\f$ (the total charge), and returns a density not a matrix.
    virtual FourierMap MakeFourierDensity(const Structure*,
                          const std::function<double(int Z, double g2)>& formFactor) const=0;

    //! \brief Hartree matrix + energy directly from the density's G-space coefficients \a rho:
    //! \f$V_H(\Delta m)=4\pi\tilde\rho/|B\Delta m|^2\f$, \f$E_H=\tfrac\Omega2\sum 4\pi|\tilde\rho|^2/G^2\f$.
    virtual hmat_t<dcmplx> Repulsion(const FourierMap& rho, double& Eh) const=0;

    // --- XC route: the basis owns the FFTs and the real-space grid; the term owns the functional. ---
    //! \brief Real-space density \f$\rho(r)\f$ on the basis's FFT grid (inverse FFT of \a rho).  The
    //! values are returned in the grid's internal order; the matching Overlap / Integral
    //! consume values in the SAME order, so the term never needs the grid points -- it just maps the
    //! functional over them (\f$v_{xc}(\rho)\f$ for the matrix, \f$\epsilon_{xc}(\rho)\rho\f$ for the energy).
    virtual rvec_t  RhoOnGrid(const FourierMap& rho) const=0;
    //! \brief Forward-FFT a real-space field \a gridValues on the FFT grid to its G-space coefficients
    //! \f$\tilde V(\Delta m)\f$ over the difference set -- the potential analogue of MakeFourierDensity (a
    //! density's rho-tilde).  So v_xc, like rho, can be carried as a FourierMap (the projected potential).
    virtual FourierMap ForwardGrid(const rvec_t& gridValues) const=0;
    //! \brief Matrix \f$\langle i|V|j\rangle = \tilde V(\Delta m)\f$ from G-space coefficients \a Vtilde
    //! (no kernel -- the overlap 3-centre is the delta).  The XC sibling of Repulsion (which carries 4pi/G^2).
    virtual hmat_t<dcmplx> Overlap(const FourierMap& Vtilde) const=0;
    //! \brief Matrix from potential values \a Vgrid on the FFT grid: ForwardGrid then Overlap(FourierMap).
    virtual hmat_t<dcmplx> Overlap(const rvec_t& Vgrid) const=0;
    //! \brief Scalar integral \f$\int f\,d^3r\f$ from values \a fgrid on the FFT grid (uniform quadrature).
    virtual double  Integral(const rvec_t& fgrid) const=0;
};

} //namespace
