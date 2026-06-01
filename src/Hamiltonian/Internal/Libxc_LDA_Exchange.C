// File: Hamiltonian/Internal/Libxc_LDA_Exchange.C Any LDA exchange potential defined in libxc.
module;
#include <iosfwd>
#include <src/xc.h>
export module qchem.Hamiltonian.Internal.Libxc_LDA_Exchange;
export import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Symmetry.Spin;
import qchem.Vector3D;

export namespace qchem::Hamiltonian
{

class Libxc_LDA_Exchange
    : public  ExFunctional
{
public:
    Libxc_LDA_Exchange(int id,const Spin&, double Ne); // See https://libxc.gitlab.io/functionals/libxc-7.0.0/

    virtual double  operator()(const rvec3_t&) const;
    virtual rvec3_t Gradient  (const rvec3_t&) const;
    virtual double  GetVxc(double ChargeDensity) const;


    virtual std::ostream& Write(std::ostream&) const;

private:
    double Ne; //# of electrons
    Spin   spin;
    xc_func_type func;
};

} //namespace
