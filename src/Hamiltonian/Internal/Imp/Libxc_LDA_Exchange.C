// File: src/Hamiltonian/Internal/Imp/Libxc_LDA_Exchange.C  Any LDA exchange functional from libxc
module;
#include <iostream>
#include <cassert>
#include <memory>
#include <src/xc.h>
module qchem.Hamiltonian.Internal.Libxc_LDA_Exchange;
import qchem.ChargeDensity;
import qchem.Streamable;

namespace qchem::Hamiltonian
{

// Unpolarized-only by construction (see the header): both functionals init XC_UNPOLARIZED, so the
// scalar single-density GetVxc below is always the correct libxc contract.
Libxc_LDA_Exchange::Libxc_LDA_Exchange(int id, double _Ne)
: Ne(_Ne)
{
    int ok = xc_func_init(&corr,     id, XC_UNPOLARIZED);   // correlation functional named by id (e.g. 7=VWN)
    assert(ok==0);
    ok     = xc_func_init(&exchange,  1, XC_UNPOLARIZED);   // Dirac exchange (LDA_X)
    assert(ok==0);
};

double Libxc_LDA_Exchange::operator()(const rvec3_t& r) const
{
    double ro = (*itsChargeDensity)(r);
    return GetVxc(ro);
}

double Libxc_LDA_Exchange::GetVxc(double rho) const
{
    // if (spin==Spin::None) rho*=0.5;
    // rho/=Ne;
    double vcorr,vexchange;
    xc_lda_vxc(&corr    , 1, &rho, &vcorr);
    xc_lda_vxc(&exchange, 1, &rho, &vexchange);
    return vcorr + vexchange;
}

rvec3_t Libxc_LDA_Exchange::Gradient(const rvec3_t& r) const
{
    assert(false);
    return rvec3_t(0,0,0);
}

std::ostream& Libxc_LDA_Exchange::Write(std::ostream& os) const
{
    os << corr.info->name << " " << exchange.info->name;
    return os;
}
} //namespace
