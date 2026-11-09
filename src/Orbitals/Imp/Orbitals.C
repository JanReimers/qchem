// File: Orbitals/Imp/Orbitals.C  Implementation unit for module qchem.Orbitals: the shared Fermi-fill numeric
// kernels (occupancy + μ-solver).  These are pure scalar math, factored out of TOrbitalsImp::TakeElectronsFermi
// so the composite cross-k global-μ fill (doc/GPWPlan1.md item 3) reuses the EXACT same bisection over the
// aggregated k-block levels -- keeping per-block and global μ bit-identical at a single k (the Γ invariant).
module;
#include <cassert>
#include <algorithm>   // std::min/std::max (not re-exported by qchem.Math)
#include <numeric>
#include <stdexcept>
#include <vector>
module qchem.Orbitals;


namespace qchem::Orbitals
{

double FermiOccupancy(double e, double mu, double kT)
{
    const double x=(e-mu)/kT;
    if (x> 40.0) return 0.0;   // overflow guard: exp(40)~2e17 already underflows f to ~0
    if (x<-40.0) return 1.0;
    return 1.0/(1.0+exp(x));
}

rvec_t ShareDegenerateLevels(const rvec_t& e, const rvec_t& g, const rvec_t& rank, double ne, double levelTol)
{
    const size_t n=e.size();
    if (g.size()!=n || (rank.size()!=0 && rank.size()!=n))
        throw std::logic_error("ShareDegenerateLevels: energies, capacities and ranks must have one entry per orbital");
    if (!(ne>=0.0)) throw std::logic_error("ShareDegenerateLevels: a negative (or NaN) electron count");
    std::vector<std::vector<size_t>> levels;
    if (levelTol<=0.0 || rank.size())   // a STATED order (MOM, a held fill) wins exactly: nothing is shared
    {
        std::vector<size_t> order(n);
        std::iota(order.begin(),order.end(),size_t(0));
        if (rank.size()) std::stable_sort(order.begin(),order.end(),[&](size_t a,size_t b){return rank[a]>rank[b];});
        // else the STORED order: a diagonalisation stores ascending energies, and a direct minimiser (GDM) stores its rotated
        // orbitals with the occupied set as the LEADING block (rotated diag(F) need not ascend: re-sorting it moved the fill off
        // GDM's manifold, 2026-11-08 panel).  A caller whose orbitals arrive in another order (the molecular reservoir, irrep
        // by irrep) states the order it wants as a rank.
        for (size_t i : order) levels.push_back({i});
    }
    else
    {
        std::vector<size_t> byE(n);
        std::iota(byE.begin(),byE.end(),size_t(0));
        std::stable_sort(byE.begin(),byE.end(),[&](size_t a,size_t b){return e[a]<e[b];});
        for (size_t i : byE)
        {
            if (levels.empty() || e[i]-e[levels.back().front()]>levelTol) levels.push_back({});
            levels.back().push_back(i);
        }
    }
    rvec_t occ(n, 0.0);
    double left=ne;
    for (const auto& L : levels)
    {
        if (left<=0.0) break;
        double G=0.0; for (size_t i : L) G+=g[i];
        if (G<=0.0) continue;
        if (left>=G*(1.0-1e-14)) { for (size_t i : L) occ[i]=g[i]; left=std::max(0.0, left-G); }   // the whole level
        else                     { for (size_t i : L) occ[i]=g[i]*left/G; left=0.0; }            // the partial one: equal shares
    }
    return occ;
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
