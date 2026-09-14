// File: src/Hamiltonian/Internal/Imp/Libxc_LDA.C  One LDA functional (exchange OR correlation) from libxc.
module;
#include <stdexcept>
#include <ostream>
#include <cassert>
#include <src/xc.h>
module qchem.Hamiltonian.Internal.Libxc_LDA;
import qchem.Streamable;

namespace qchem::Hamiltonian
{

// Unpolarized-only by construction (see the header): the functional inits XC_UNPOLARIZED, so the scalar
// single-density GetVxc/GetEpsXc below are always the correct libxc contract.
Libxc_LDA::Libxc_LDA(int id)
{
    int ok = xc_func_init(&itsFunc, id, XC_UNPOLARIZED);
    assert(ok==0);
}

Libxc_LDA::~Libxc_LDA()
{
    xc_func_end(&itsFunc);
}

double Libxc_LDA::GetVxc(double rho) const
{
    double v;
    xc_lda_vxc(&itsFunc, 1, &rho, &v);
    return v;
}

double Libxc_LDA::GetEpsXc(double rho) const
{
    double eps;
    xc_lda_exc(&itsFunc, 1, &rho, &eps);   // energy density per particle; E = integral eps rho
    return eps;
}

// zeta=0 only (see the header): an unpolarized run's folded doublet hands over rho/2 per channel, which is
// exactly equal; anything else is a polarized density this wrapper cannot evaluate.
static void RequireUnpolarized(double up, double dn)
{
    if (up!=dn)
        throw std::runtime_error("Libxc_LDA: this wrapper is scalar (XC_UNPOLARIZED) by construction and was "
                                 "handed a spin-polarized density (rho_up != rho_down).  Use XC::DiracVWN for "
                                 "polarized LDA.");
}
double Libxc_LDA::GetVxc  (double up, double dn, const Spin&) const {RequireUnpolarized(up,dn); return GetVxc  (up+dn);}
double Libxc_LDA::GetEpsXc(double up, double dn, const Spin&) const {RequireUnpolarized(up,dn); return GetEpsXc(up+dn);}

std::ostream& Libxc_LDA::Write(std::ostream& os) const
{
    return os << itsFunc.info->name;
}

} //namespace
