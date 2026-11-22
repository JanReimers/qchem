// File: Hamiltonian/Internal/Imp/GGAFunctional.C  The forgotten-gradient detector's body (moved out of the interface at GGA stage 2:
// a body in an interface rewrites its BMI on every edit, CLAUDE.md).
module;
#include <string>
#include <stdexcept>
module qchem.Hamiltonian.Internal.GGAFunctional;

namespace qchem::Hamiltonian
{

void RequireGradientSampled(const GGAInput& in, const std::string& route)
{
    const size_t n=in.size();
    if (n==0) return;
    double sMax=0.0, rMin=in.up[0]+in.dn[0], rMax=rMin;
    for (size_t g=0;g<n;++g)
    {
        const double s=in.sigUU[g]+in.sigDD[g], r=in.up[g]+in.dn[g];
        if (s>sMax) sMax=s;
        if (r<rMin) rMin=r;
        if (r>rMax) rMax=r;
    }
    if (sMax==0.0 && rMax-rMin>1e-12*rMax)
        throw std::logic_error("GGA on "+route+": sigma is zero at all "+std::to_string(n)+" points while rho varies ("+std::to_string(rMin)+
                               " .. "+std::to_string(rMax)+") -- the density gradient was never computed (or cached), and the run would be an LDA");
}

} //namespace
