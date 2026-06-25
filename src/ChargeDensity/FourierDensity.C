// File: ChargeDensity/FourierDensity.C  A density's reciprocal-space (G-space) coefficients.
//
// For a PERIODIC system rho-tilde(G) is the density's NATIVE representation -- you diagonalize in
// G-space and store rho-tilde; the real-space rho(r) is the FFT-derived view, not the reverse.  So this
// is core charge-density functionality for solids, not an escape hatch: it is the dual of
// BasisSet::FourierDFT_IBS, and lets the plane-wave basis assemble the Hartree (FFT-free) and XC (FFT)
// matrices from rho-tilde directly instead of O(Npts*n^2) pointwise sampling.  A FINITE (molecular)
// density has no reciprocal-lattice Fourier series, so it does not provide one (cf. Structure::isFinite);
// a DFT term reaches a periodic density's rho-tilde by dynamic_cast (abstract->abstract).  The composite
// density sums the per-block (BZ-weighted) rho-tilde over the k-mesh.
module;
export module qchem.ChargeDensity.FourierDensity;
export import qchem.FourierMap;

export namespace qchem::ChargeDensity
{

class FourierDensity
{
public:
    virtual ~FourierDensity() {}
    //! \brief The density's reciprocal-space coefficients \f$\tilde\rho(\Delta m)\f$ (BZ-weighted; the
    //! composite sums \f$\sum_k w_k\tilde\rho_k\f$ over the k-mesh).  Keyed by reciprocal-index difference.
    virtual FourierMap GetFourierDensity() const=0;
};

} //namespace
