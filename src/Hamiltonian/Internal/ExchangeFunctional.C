// File: ExchangeFunctional.C   Exchange potential for DFT.
export module qchem.Hamiltonian.Internal.ExFunctional;
import qchem.ChargeDensity;
import qchem.ScalarFunction;
import qchem.Streamable;
import qchem.Symmetry.Spin;   // Spin -- the requested channel of the spin-native correlation face

export namespace qchem::Hamiltonian
{

using ChargeDensity::rDM_CD;
using ChargeDensity::rChargeDensity;

class ExFunctional
    : public virtual Streamable
    , public virtual ScalarFunction<double>
{
public:
    ExFunctional(               );

    virtual void   InsertChargeDensity(const rChargeDensity*);
    virtual rvec_t GetVxcs(const rvec_t& ChargeDensities) const;
    virtual double GetVxc(                double ChargeDensity) const=0;
    //! \brief Energy density per particle \f$\varepsilon_{xc}(\rho)\f$, so \f$E_{xc}=\int\varepsilon_{xc}\rho\,d^3r\f$.
    //! Default is the EXCHANGE virial \f$\varepsilon_x=\tfrac34 v_x\f$ (exact for Dirac/Slater exchange).
    //! CORRELATION functionals MUST override: \f$\varepsilon_c\neq\tfrac34 v_c\f$ (differs ~15%).
    virtual double GetEpsXc(              double ChargeDensity) const {return 0.75*GetVxc(ChargeDensity);}
    virtual void   SetPolarized(bool p) {isPolarized=p;}

protected:

    const rChargeDensity* itsChargeDensity;
    bool            isPolarized;
};

//! \brief Spin-native correlation face (no data; an abstract capability mixin).
//!
//! Correlation does NOT separate by spin channel the way Slater exchange does: \f$v_c^\sigma\f$ and
//! \f$\varepsilon_c\f$ COUPLE both densities (through \f$r_s(\rho_\uparrow+\rho_\downarrow)\f$ and
//! \f$\zeta\f$), so they cannot be expressed through the single-density \c ExFunctional::GetVxc face that
//! channel-separable exchange uses.  A correlation functional that supports polarized DFT implements this
//! two-channel face; \c FittedVcorrPol consumes it.  Unpolarized is the \f$\rho_\uparrow=\rho_\downarrow\f$
//! collapse (so \c GetVc(h,h,s)==\c ExFunctional::GetVxc(2h)).
class SpinCorrelation
{
public:
    virtual ~SpinCorrelation() {}
    //! Correlation energy density per particle \f$\varepsilon_c(\rho_\uparrow,\rho_\downarrow)\f$.
    virtual double GetEpsC(double rhoUp, double rhoDown) const=0;
    //! Channel correlation potential \f$v_c^\sigma=\varepsilon_c+\rho\,\partial\varepsilon_c/\partial\rho_\sigma\f$.
    virtual double GetVc  (double rhoUp, double rhoDown, const Spin&) const=0;
};

} //namespace