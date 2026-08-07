// File: SlaterExchange.C Slater exchange potential.
module;
#include <iosfwd>
export module qchem.Hamiltonian.Internal.SlaterExchange;
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Symmetry.Spin;

export namespace qchem::Hamiltonian
{

class SlaterExchange
    : public  ExFunctional
{
public:
    SlaterExchange(               );
    SlaterExchange(double theAlpha);
    SlaterExchange(double theAlpha, const Spin&);

    //! \f$v_x(\rho)=-3\alpha(3\rho/4\pi)^{1/3}\f$.  \c itsSpin is what expresses polarization here:
    //! Spin::None halves rho first (the closed-shell collapse).  (V1.13 removed the FIELD face --
    //! op(r)/Gradient(r) -- and with it the density pointer and the never-set isPolarized bool.)
    virtual double GetVxc(double ChargeDensity) const;


    virtual std::ostream& Write(std::ostream&) const;

private:
    double itsAlpha;
    Spin   itsSpin;
};

} //namespace
