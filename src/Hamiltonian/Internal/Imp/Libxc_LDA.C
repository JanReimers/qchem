// File: src/Hamiltonian/Internal/Imp/Libxc_LDA.C  One LDA functional (exchange OR correlation) from libxc.
module;
#include <ostream>
#include <cassert>
#include <src/xc.h>
module qchem.Hamiltonian.Internal.Libxc_LDA;
import qchem.Streamable;

namespace qchem::Hamiltonian
{

// Two inits of the same id: the scalar face keeps the XC_UNPOLARIZED contract, the two-channel face uses
// XC_POLARIZED (see the header).
Libxc_LDA::Libxc_LDA(int id)
{
    int ok = xc_func_init(&itsFunc, id, XC_UNPOLARIZED);
    assert(ok==0);
    ok = xc_func_init(&itsPol, id, XC_POLARIZED);
    assert(ok==0);
    itsSeparable = itsFunc.info->kind==XC_EXCHANGE;
}

Libxc_LDA::~Libxc_LDA()
{
    xc_func_end(&itsFunc);
    xc_func_end(&itsPol);
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

// A negative channel (band-limited ringing) is arithmetic, not physics: clamp EACH channel to >= 0, the rule
// libxc applies itself below its density threshold and the one SlaterExchange's spin face follows (a negative
// channel contributes nothing, the other channel is untouched).  The wrapper cannot tell an exchange id from a
// correlation id, so it does not use VWN5's total-preserving zeta clamp, which would move an exchange channel.
void Libxc_LDA::EvalPolarized(double up, double dn, double& eps, double v[2]) const
{
    double r[2]={up>0.0 ? up : 0.0, dn>0.0 ? dn : 0.0};
    if (r[0]+r[1]<=0.0) { eps=0.0; v[0]=v[1]=0.0; return; }
    xc_lda_exc_vxc(&itsPol, 1, r, &eps, v);
}

// Exact zeta=0 stays on the scalar path (byte-identical to the unpolarized-only wrapper); an EXCHANGE id is the
// exact spin-scaling relation (the base class's channel-separable face over the scalar libxc call); otherwise
// XC_POLARIZED.
double Libxc_LDA::GetVxc(double up, double dn, const Spin& s) const
{
    if (up==dn) return GetVxc(up+dn);
    if (itsSeparable) return ExFunctional::GetVxc(up<0.0 ? 0.0 : up, dn<0.0 ? 0.0 : dn, s);
    double eps, v[2]; EvalPolarized(up, dn, eps, v);
    return s==Spin::Down ? v[1] : v[0];
}
double Libxc_LDA::GetEpsXc(double up, double dn, const Spin& s) const
{
    if (up==dn) return GetEpsXc(up+dn);
    if (itsSeparable) return ExFunctional::GetEpsXc(up<0.0 ? 0.0 : up, dn<0.0 ? 0.0 : dn, s);
    double eps, v[2]; EvalPolarized(up, dn, eps, v);
    return eps;                                         // per particle of the TOTAL: the same for both channels
}

std::ostream& Libxc_LDA::Write(std::ostream& os) const
{
    return os << itsFunc.info->name;
}

} //namespace
