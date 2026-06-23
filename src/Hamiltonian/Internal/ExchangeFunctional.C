// File: ExchangeFunctional.C   Exchange potential for DFT.
export module qchem.Hamiltonian.Internal.ExFunctional;
import qchem.ChargeDensity;
import qchem.ScalarFunction;
import qchem.Streamable;

export namespace qchem::Hamiltonian
{

using ChargeDensity::DM_CD;

class ExFunctional
    : public virtual Streamable
    , public virtual ScalarFunction<double>
{
public:
    ExFunctional(               );

    virtual void   InsertChargeDensity(const DM_CD*);
    virtual rvec_t GetVxcs(const rvec_t& ChargeDensities) const;
    virtual double GetVxc(                double ChargeDensity) const=0;
    //! \brief Energy density per particle \f$\varepsilon_{xc}(\rho)\f$, so \f$E_{xc}=\int\varepsilon_{xc}\rho\,d^3r\f$.
    //! Default is the EXCHANGE virial \f$\varepsilon_x=\tfrac34 v_x\f$ (exact for Dirac/Slater exchange).
    //! CORRELATION functionals MUST override: \f$\varepsilon_c\neq\tfrac34 v_c\f$ (differs ~15%).
    virtual double GetEpsXc(              double ChargeDensity) const {return 0.75*GetVxc(ChargeDensity);}
    virtual void   SetPolarized(bool p) {isPolarized=p;}

protected:

    const DM_CD* itsChargeDensity;
    bool            isPolarized;
};

} //namespace 