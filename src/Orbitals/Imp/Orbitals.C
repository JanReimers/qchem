// File: Orbitals/Imp/Orbitals.C  Implementation unit for module qchem.Orbitals: the shared Fermi-fill numeric
// kernels (occupancy + μ-solver).  These are pure scalar math, factored out of TOrbitalsImp::TakeElectronsFermi
// so the composite cross-k global-μ fill (doc/GPWPlan1.md item 3) reuses the EXACT same bisection over the
// aggregated k-block levels -- keeping per-block and global μ bit-identical at a single k (the Γ invariant).
module;
#include <cassert>
#include <algorithm>   // std::min/std::max (not re-exported by qchem.Math)
module qchem.Orbitals;

import qchem.Math;    // exp, log, std::min/max, fabs
import qchem.Blaze;   // rvec_t

namespace qchem::Orbitals
{

double FermiOccupancy(double e, double mu, double kT)
{
    const double x=(e-mu)/kT;
    if (x> 40.0) return 0.0;   // overflow guard: exp(40)~2e17 already underflows f to ~0
    if (x<-40.0) return 1.0;
    return 1.0/(1.0+exp(x));
}

double FermiLevel(const rvec_t& e, const rvec_t& g, double target, double kT)
{
    assert(kT>0.0);
    const size_t n=e.size();
    assert(n>0 && g.size()==n);
    // Total capacity Σ g_i must cover the target (else too few levels -- same contract as aufbau's assert),
    // and the level range brackets μ once padded by ±50 kT.
    double cap=0.0, emin=e[0], emax=e[0];
    for (size_t i=0;i<n;++i){ cap+=g[i]; emin=std::min(emin,e[i]); emax=std::max(emax,e[i]); }
    assert(target<=cap+1e-9 && "Fermi fill: too few orbitals for the electron count (basis/mesh too small)");
    auto count=[&](double mu)->double
    {
        double N=0.0;
        for (size_t i=0;i<n;++i) N += g[i]*FermiOccupancy(e[i],mu,kT);
        return N;
    };
    // Bisect μ.  N(lo)≈0≤target and N(hi)≈cap≥target by the padding + the capacity assert above.
    double lo=emin-50.0*kT, hi=emax+50.0*kT;
    for (int it=0; it<100 && (hi-lo)>1e-14*(1.0+fabs(hi)); ++it)
    {
        const double mu=0.5*(lo+hi);
        if (count(mu)<target) lo=mu; else hi=mu;
    }
    return 0.5*(lo+hi);
}

} // namespace qchem::Orbitals
