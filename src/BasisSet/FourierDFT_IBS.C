// File: BasisSet/FourierDFT_IBS.C  Abstract G-space (plane-wave) DFT-assembly capability.
//
// The dual of DFTPotential_IBS: instead of integrating real-space ScalarFunctions on a mesh, a plane-wave
// basis and a plane-wave density share G-space coefficients (a FourierMap, rho-tilde / V-tilde).  This
// replaces the O(Npts*n^2) pointwise density sampling with one O(n^2) projection D->rho-tilde plus an
// FFT-free Hartree solve (and, later, FFT-based XC).  A DFT term reaches it the same sanctioned way it
// reaches DFTPotential_IBS -- holding the abstract orbital basis and dynamic_cast-ing UP to this capability.
module;
export module qchem.BasisSet.FourierDFT_IBS;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.FourierMap;
import qchem.Types;   // hmat_t<dcmplx>

export namespace BasisSet
{

//! \brief A basis that assembles DFT matrices in G-space from a density's Fourier coefficients.  Only
//! the complex (plane-wave) path realizes this; a density exposes its rho-tilde via the matching
//! ChargeDensity FourierDensity sub-interface, and the BZ sum is the per-block rho-tilde summed.
class FourierDFT_IBS
    : public virtual Orbital_1E_IBS<dcmplx>
{
public:
    //! \brief \f$\tilde\rho(\Delta m)=\frac1\Omega\sum_{G_i-G_j=\Delta m}D_{ij}\f$ for a density matrix
    //! \a D in THIS plane-wave block (one \f$O(n^2)\f$ accumulation over the difference set).
    virtual FourierMap MakeFourierDensity(const hmat_t<dcmplx>& D) const=0;

    //! \brief Hartree matrix + energy directly from the density's G-space coefficients \a rho:
    //! \f$V_H(\Delta m)=4\pi\tilde\rho/|B\Delta m|^2\f$, \f$E_H=\tfrac\Omega2\sum 4\pi|\tilde\rho|^2/G^2\f$.
    virtual hmat_t<dcmplx> IntegralHartree(const FourierMap& rho, double& Eh) const=0;

    // --- XC route: the basis owns the FFTs and the real-space grid; the term owns the functional. ---
    //! \brief Real-space density \f$\rho(r)\f$ on the basis's FFT grid (inverse FFT of \a rho).  The
    //! values are returned in the grid's internal order; the matching IntegralPotentialGrid / IntegralGrid
    //! consume values in the SAME order, so the term never needs the grid points -- it just maps the
    //! functional over them (\f$v_{xc}(\rho)\f$ for the matrix, \f$\epsilon_{xc}(\rho)\rho\f$ for the energy).
    virtual rvec_t  RhoOnGrid(const FourierMap& rho) const=0;
    //! \brief Matrix \f$\langle i|V|j\rangle\f$ from potential values \a Vgrid on the FFT grid (forward
    //! FFT to \f$\tilde V(\Delta m)\f$, then assemble).
    virtual hmat_t<dcmplx> IntegralPotentialGrid(const rvec_t& Vgrid) const=0;
    //! \brief Scalar integral \f$\int f\,d^3r\f$ from values \a fgrid on the FFT grid (uniform quadrature).
    virtual double  IntegralGrid(const rvec_t& fgrid) const=0;
};

} //namespace
