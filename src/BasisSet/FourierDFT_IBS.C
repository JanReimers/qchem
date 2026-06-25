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
};

} //namespace
