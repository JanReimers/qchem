// File: Hamiltonian/Internal/RealSpaceResidual.C  The capability of a term that owns a real-space quadrature.
module;
export module qchem.Hamiltonian.Internal.RealSpaceResidual;
export import qchem.ChargeDensity;

export namespace qchem::Hamiltonian
{

//! \brief A term that owns the run's real-space quadrature measures two densities on it.  The Hamiltonian holds the term that
//! answers, and \c tHamiltonian::DeltaRho is the only client: nothing above the Hamiltonian sees that a term did the measuring.
template <class T> class tRealSpaceResidual
{
public:
    virtual ~tRealSpaceResidual() = default;
    //! \f$\int(|\rho^a_\uparrow-\rho^b_\uparrow|+|\rho^a_\downarrow-\rho^b_\downarrow|)\,d^3r\f$ (electrons) on this term's quadrature.
    virtual double L1Distance(const ChargeDensity::tChargeDensity<T>* a, const ChargeDensity::tChargeDensity<T>* b) const=0;
};

} //namespace
