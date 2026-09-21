// File: Hamiltonian/Internal/Imp/ACBN0.C  ACBN0 -- implementation.  See the module interface for the design.
module;
#include <cassert>
#include <complex>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
module qchem.Hamiltonian.Internal.ACBN0;
import qchem.Blaze;
import qchem.Math;

namespace qchem::Hamiltonian
{

ACBN0::ACBN0(const HubbardProjection& term) : itsTerm(term)
{
    // Both channels always exist here: an unpolarized run's folded doublet is split into them, and the
    // paper's eqs 12/13 are written per channel.
    for (Spin s : {Spin::Up, Spin::Down}) itsChannels[s];
}

template <class U> void ACBN0::EnsureIntegrals(const BasisSet::Orbital_DFT_IBS<U,dcmplx>& block)
{
    if (!itsERI.empty()) return;
    const auto* src=dynamic_cast<const BasisSet::BareCoulombSource*>(&block);
    if (!src) throw std::runtime_error("ACBN0: the orbital block cannot deliver bare two-electron integrals over a "
                                       "function subset (no BareCoulombSource face) -- no on-site U from it");
    const auto cols=itsTerm.ManifoldFunctions(block);
    for (const auto& c : cols) itsERI.push_back(src->BareCoulomb(c));
    for (auto& [s,ch] : itsChannels)
    {
        ch.P.clear(); ch.Pbare.clear();
        for (const auto& c : cols) { ch.P.emplace_back(c.size(), c.size(), 0.0); ch.Pbare.emplace_back(c.size(), c.size(), 0.0); }
    }
}

template <class U> void ACBN0::AccumulateT(const BasisSet::Orbital_DFT_IBS<U,dcmplx>& block, const Spin& s, double w,
                                           const mat_t<U>& C, const rvec_t& f)
{
    if (C.columns()!=f.size()) throw std::invalid_argument("ACBN0::Accumulate: one occupation per coefficient column");
    EnsureIntegrals(block);
    if (f.size()==0) return;
    const std::vector<mat_t<U>> ell=itsTerm.LowdinCoefficients(block, C);     // per manifold: m x nOrb
    const size_t nM=itsERI.size(), nOrb=f.size();
    // The renormalised occupation of each orbital on each manifold's EQUIVALENCE set (same species and l).
    std::vector<std::vector<size_t>> eq(nM);
    for (size_t M=0;M<nM;M++) eq[M]=itsTerm.EquivalentManifolds(M);
    auto norm2=[&](size_t M, size_t i){ double t=0; for (size_t a=0;a<ell[M].rows();a++) t+=std::norm(ell[M](a,i)); return t; };
    // Which channels this block feeds, and with what share of f.
    std::vector<std::pair<Spin,double>> feed;
    if (s==Spin::None) feed={{Spin::Up,0.5},{Spin::Down,0.5}}; else feed={{s,1.0}};
    for (size_t M=0;M<nM;M++)
    {
        const size_t m=ell[M].rows();
        for (size_t i=0;i<nOrb;i++)
        {
            double Nbar=0.0; for (size_t K : eq[M]) Nbar+=norm2(K,i);
            for (const auto& [sp,share] : feed)
            {
                Channel& ch=itsChannels[sp];
                const double wf=w*f[i]*share;
                for (size_t a=0;a<m;a++) for (size_t b=0;b<m;b++)
                {
                    const double x=std::real(ell[M](a,i)*blazem::conjs(ell[M](b,i)));   // Hermitian after the k-sum
                    ch.P    [M](a,b)+=wf*Nbar*x;
                    ch.Pbare[M](a,b)+=wf*x;
                }
            }
        }
    }
    itsFed=true;
}

void ACBN0::Accumulate(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& block, const Spin& s, double w, const mat_t<double>& C, const rvec_t& f)
{AccumulateT<double>(block,s,w,C,f);}
void ACBN0::Accumulate(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& block, const Spin& s, double w, const mat_t<dcmplx>& C, const rvec_t& f)
{AccumulateT<dcmplx>(block,s,w,C,f);}

void ACBN0::Averages(const rmat_t& Pa, const rmat_t& Pb, const BasisSet::ERI4Block& eri, double& Ubar, double& Jbar)
{
    const size_t m=eri.Size();
    assert(Pa.rows()==m && Pb.rows()==m);
    const rmat_t Pt=Pa+Pb;
    double numU=0.0, numJ=0.0;
    for (size_t a=0;a<m;a++) for (size_t b=0;b<m;b++)
    {
        if (Pt(a,b)==0.0 && Pa(a,b)==0.0 && Pb(a,b)==0.0) continue;
        for (size_t c=0;c<m;c++) for (size_t d=0;d<m;d++)
        {
            numU+=Pt(a,b)*Pt(c,d)*eri(a,b,c,d);                       // (12|34)
            numJ+=(Pa(a,b)*Pa(c,d)+Pb(a,b)*Pb(c,d))*eri(a,d,c,b);     // (14|32)
        }
    }
    double Na=0, Nb=0, Na2=0, Nb2=0;
    for (size_t a=0;a<m;a++) { Na+=Pa(a,a); Nb+=Pb(a,a); Na2+=Pa(a,a)*Pa(a,a); Nb2+=Pb(a,a)*Pb(a,a); }
    const double denU=(Na+Nb)*(Na+Nb)-Na2-Nb2;      // Sum_{m!=m'} NaNa' + Sum NaNb' + Sum NbNa' + Sum_{m!=m'} NbNb'
    const double denJ=Na*Na-Na2+Nb*Nb-Nb2;          // Sum_{m!=m'} (NaNa' + NbNb')
    Ubar = denU>1e-12 ? numU/denU : 0.0;
    Jbar = denJ>1e-12 ? numJ/denJ : 0.0;
}

std::vector<HubbardEstimate> ACBN0::Evaluate() const
{
    if (!itsFed) throw std::logic_error("ACBN0::Evaluate before any orbital block was fed");
    std::vector<HubbardEstimate> out;
    const Channel& up=itsChannels.at(Spin::Up);
    const Channel& dn=itsChannels.at(Spin::Down);
    for (size_t M=0;M<itsERI.size();M++)
    {
        HubbardEstimate e;
        e.site=itsTerm.Manifolds()[M].site; e.l=itsTerm.Manifolds()[M].l;
        Averages(up.P[M],     dn.P[M],     itsERI[M], e.Ubar,     e.Jbar);
        Averages(up.Pbare[M], dn.Pbare[M], itsERI[M], e.UbarBare, e.JbarBare);
        const size_t m=itsERI[M].Size();
        e.Nup=rvec_t(m); e.Ndn=rvec_t(m);
        for (size_t a=0;a<m;a++) { e.Nup[a]=up.P[M](a,a); e.Ndn[a]=dn.P[M](a,a); e.chargeUp+=e.Nup[a]; e.chargeDn+=e.Ndn[a]; }
        out.push_back(std::move(e));
    }
    return out;
}

std::ostream& ACBN0::Write(std::ostream& os) const
{
    const double eV=27.211386245988;
    for (const HubbardEstimate& e : Evaluate())
        os<<"[ACBN0] site "<<e.site<<" l="<<e.l<<": U="<<std::fixed<<std::setprecision(4)<<e.Ubar*eV<<" J="<<e.Jbar*eV
          <<" U_eff="<<e.Ueff()*eV<<" eV  (bare, unrenormalised: U="<<e.UbarBare*eV<<" J="<<e.JbarBare*eV<<")"
          <<"  N_ren up/dn="<<e.chargeUp<<"/"<<e.chargeDn<<std::endl;
    return os;
}

} // namespace
