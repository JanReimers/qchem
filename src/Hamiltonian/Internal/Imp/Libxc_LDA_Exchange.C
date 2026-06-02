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

Libxc_LDA_Exchange::Libxc_LDA_Exchange(int id, const Spin& s, double _Ne)
: Ne(_Ne), spin(s)
{
    int ok= spin==Spin::None ? 
    xc_func_init(&corr,id, XC_UNPOLARIZED) 
    : xc_func_init(&corr,id, XC_POLARIZED);
    assert(ok==0);
    ok= spin==Spin::None ? 
    xc_func_init(&exchange,1, XC_UNPOLARIZED) 
    : xc_func_init(&exchange,1, XC_POLARIZED);
    assert(ok==0);
    // std::cout <<"exchange.info->n_ext_params=" << exchange.info->ext_params.n << std::endl;
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
    // std::cout << "Ne,rho,vxc = " << Ne << " " << rho << " " << ret*Ne << std::endl; 
    return (vcorr+vexchange*1.006613); //Fudge factor for Z=36!
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
