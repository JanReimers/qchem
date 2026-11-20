// File: src/Hamiltonian/Internal/Imp/Libxc_GGA.C  One GGA functional (exchange OR correlation) from libxc.
module;
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <src/xc.h>
module qchem.Hamiltonian.Internal.Libxc_GGA;

namespace qchem::Hamiltonian
{

Libxc_GGA::Libxc_GGA(int id)
{
    // THROWS, not an assert: the id is user-reachable (XCFunctional), and an assert is dead in Release.
    if (xc_func_init(&itsPol, id, XC_POLARIZED)!=0)
        throw std::invalid_argument("Libxc_GGA: libxc has no functional with id "+std::to_string(id)+" (or cannot initialise it spin-polarized)");
    if (itsPol.info->family!=XC_FAMILY_GGA)
    {
        const std::string name=itsPol.info->name;
        xc_func_end(&itsPol);
        throw std::invalid_argument("Libxc_GGA: libxc functional id "+std::to_string(id)+" ("+name+") is not a pure GGA -- an LDA, a "
                                    "hybrid (exact exchange) or a meta-GGA needs another evaluator");
    }
}

Libxc_GGA::~Libxc_GGA() {xc_func_end(&itsPol);}

bool Libxc_GGA::IsExchangeCorrelation() const {return itsPol.info->kind==XC_EXCHANGE_CORRELATION;}

[[noreturn]] static void NoPointFace(const xc_func_type& f, const char* face)
{
    throw std::logic_error(std::string("Libxc_GGA (")+f.info->name+"): "+face+" asked without the density gradient -- a GGA answers "
                           "only the batch face with the sigma triple; the route that called this should have been refused by policy");
}
double Libxc_GGA::GetVxc  (double) const {NoPointFace(itsPol, "GetVxc(rho)");}
double Libxc_GGA::GetEpsXc(double) const {NoPointFace(itsPol, "GetEpsXc(rho)");}
double Libxc_GGA::GetVxc  (double, double, const Spin&) const {NoPointFace(itsPol, "GetVxc(up,dn)");}
double Libxc_GGA::GetEpsXc(double, double, const Spin&) const {NoPointFace(itsPol, "GetEpsXc(up,dn)");}
double Libxc_GGA::GetFxc  (double, double, const Spin&, const Spin&) const {NoPointFace(itsPol, "GetFxc");}

// libxc's polarized layout: rho[2g+s], sigma[3g+{uu,ud,dd}], zk[g] (per PARTICLE of the total), vrho[2g+s], vsigma[3g+{uu,ud,dd}].
void Libxc_GGA::Evaluate(const XCInput& in, XCOutput& out) const
{
    if (!in.HasGradient())
        throw std::logic_error(std::string("Libxc_GGA (")+itsPol.info->name+")::Evaluate: the input carries no sigma triple -- the term "
                               "must sample the density gradient for a functional that ConsumesDensityGradient");
    const size_t n=in.size();
    if (n==0) return;
    std::vector<double> rho(2*n), sigma(3*n), zk(n), vrho(2*n), vsigma(3*n);
    for (size_t g=0;g<n;++g)
    {
        rho[2*g  ]=in.up[g]>0.0 ? in.up[g] : 0.0;
        rho[2*g+1]=in.dn[g]>0.0 ? in.dn[g] : 0.0;
        sigma[3*g  ]=in.sigUU[g];
        sigma[3*g+1]=in.sigUD[g];
        sigma[3*g+2]=in.sigDD[g];
    }
    xc_gga_exc_vxc(&itsPol, n, rho.data(), sigma.data(), zk.data(), vrho.data(), vsigma.data());
    for (size_t g=0;g<n;++g)
    {
        out.e   [g]+=zk[g]*(rho[2*g]+rho[2*g+1]);   // per particle -> per unit volume
        out.vUp [g]+=vrho[2*g  ];
        out.vDn [g]+=vrho[2*g+1];
        out.vsUU[g]+=vsigma[3*g  ];
        out.vsUD[g]+=vsigma[3*g+1];
        out.vsDD[g]+=vsigma[3*g+2];
    }
}

std::ostream& Libxc_GGA::Write(std::ostream& os) const
{
    return os << itsPol.info->name << " (libxc thresholds: density " << itsPol.dens_threshold << ", sigma " << itsPol.sigma_threshold << ")";
}

} //namespace
