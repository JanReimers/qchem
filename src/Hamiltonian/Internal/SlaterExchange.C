// File: SlaterExchange.C Slater exchange potential.
module;
#include <iosfwd>
export module qchem.Hamiltonian.Internal.SlaterExchange;
import qchem.Hamiltonian.Internal.ExFunctional;

export namespace qchem::Hamiltonian
{

class SlaterExchange
    : public  ExFunctional
{
public:
    SlaterExchange(               );
    SlaterExchange(double theAlpha);
    //! The closed-shell (\f$\zeta=0\f$) potential of the TOTAL density: \f$v_x(\rho)=-3\alpha(3(\rho/2)/4\pi)^{1/3}\f$
    //! -- each channel holds \f$\rho/2\f$.  The per-channel \f$v_x^\sigma(\rho_\sigma)\f$ is the base
    //! class's channel-separable default, \f$v_x(2\rho_\sigma)\f$; the Spin-tagged ctor that used to select
    //! between the two is gone (V1.37 step 3) -- a functional does not know which subgroup is imposed.
    using ExFunctional::GetVxc;     // keep the two-channel face visible beside the scalar override
    using ExFunctional::GetEpsXc;
    virtual double GetVxc(double ChargeDensity) const;


    virtual std::ostream& Write(std::ostream&) const;

private:
    double itsAlpha;
};

} //namespace
