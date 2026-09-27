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

ACBN0::ACBN0(HubbardProjection& term) : itsTerm(term)
{
    // Both channels always exist here: an unpolarized run's folded doublet is split into them, and the
    // paper's eqs 12/13 are written per channel.
    for (Spin s : {Spin::Up, Spin::Down}) itsChannels[s];
}

template <class U> void ACBN0::EnsureIntegrals(const BasisSet::Orbital_DFT_IBS<U,dcmplx>& block)
{
    if (!itsERI.empty()) return;
    itsERI=itsTerm.ManifoldIntegrals(block);                    // the term knows what its manifolds are built over
    for (auto& [s,ch] : itsChannels)
    {
        ch.P.clear(); ch.Pbare.clear(); ch.L.clear(); ch.Lbare.clear();
        for (const auto& e : itsERI)
        {
            const size_t m=e.Size();
            ch.P.emplace_back(m, m, 0.0); ch.Pbare.emplace_back(m, m, 0.0);
            ch.L.emplace_back(m, m, 0.0); ch.Lbare.emplace_back(m, m, 0.0);
        }
    }
}

template <class U> void ACBN0::AccumulateT(const BasisSet::Orbital_DFT_IBS<U,dcmplx>& block, const Spin& s, double w,
                                           const mat_t<U>& C, const rvec_t& f)
{
    if (C.columns()!=f.size()) throw std::invalid_argument("ACBN0::Accumulate: one occupation per coefficient column");
    EnsureIntegrals(block);
    if (f.size()==0) return;
    const std::vector<mat_t<U>> ell=itsTerm.ProjectorAmplitudes(block, C);     // per manifold: m x nOrb (Löwdin basis)
    const std::vector<mat_t<U>> coef=itsTerm.ManifoldCoefficients(block, C);  // per manifold: m x nOrb (on the manifold's functions)
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
                    const double x=std::real(coef[M](a,i)*blazem::conjs(coef[M](b,i)));   // on the manifold's functions; Hermitian after the k-sum
                    const double y=std::real(ell[M](a,i)*blazem::conjs(ell[M](b,i)));             // Löwdin basis
                    ch.P    [M](a,b)+=wf*Nbar*x;  ch.Pbare[M](a,b)+=wf*x;
                    ch.L    [M](a,b)+=wf*Nbar*y;  ch.Lbare[M](a,b)+=wf*y;
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

void ACBN0::Averages(const rmat_t& Pa, const rmat_t& Pb, const rvec_t& Nav, const rvec_t& Nbv,
                     const BasisSet::ERI4Block& eri, double& Ubar, double& Jbar)
{
    const size_t m=eri.Size();
    assert(Pa.rows()==m && Pb.rows()==m && Nav.size()==m && Nbv.size()==m);
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
    for (size_t a=0;a<m;a++) { Na+=Nav[a]; Nb+=Nbv[a]; Na2+=Nav[a]*Nav[a]; Nb2+=Nbv[a]*Nbv[a]; }
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
        const size_t m=itsERI[M].Size();
        rvec_t NupBare(m), NdnBare(m);
        e.Nup=rvec_t(m); e.Ndn=rvec_t(m);
        for (size_t a=0;a<m;a++)
        {
            e.Nup[a]=up.L[M](a,a); e.Ndn[a]=dn.L[M](a,a); e.chargeUp+=e.Nup[a]; e.chargeDn+=e.Ndn[a];
            NupBare[a]=up.Lbare[M](a,a); NdnBare[a]=dn.Lbare[M](a,a);
        }
        // THE SCREENING LIVES IN THIS ASYMMETRY (paper eqs 10a-10c, re-read 2026-09-21 after a symmetric first
        // version cancelled it): the NUMERATOR carries the renormalised P-bar twice (∝ N-bar^2), the pair-count
        // DENOMINATOR the UNRENORMALISED populations (eq 10c has no N-bar) -- so a manifold the KS states only
        // partly live in (N-bar < 1) is screened by ~N-bar^2, and a d^0 site's U goes to zero with its d weight.
        Averages(up.P[M],     dn.P[M],     NupBare, NdnBare, itsERI[M], e.Ubar,     e.Jbar);
        Averages(up.Pbare[M], dn.Pbare[M], NupBare, NdnBare, itsERI[M], e.UbarBare, e.JbarBare);
        e.chargeUpBare=0; e.chargeDnBare=0;
        for (size_t a=0;a<m;a++) { e.chargeUpBare+=NupBare[a]; e.chargeDnBare+=NdnBare[a]; }
        out.push_back(std::move(e));
    }
    return out;
}

void ACBN0::Apply(const std::vector<HubbardEstimate>& e)
{
    if (e.size()!=itsTerm.Manifolds().size()) throw std::invalid_argument("ACBN0::Apply: one estimate per manifold");
    for (size_t M=0;M<e.size();M++) itsTerm.SetU(M, e[M].Ueff());
    // Start the next feed clean: the orbitals of the run just estimated belong to the previous U.
    for (auto& [s,ch] : itsChannels) for (auto* v : {&ch.P,&ch.Pbare,&ch.L,&ch.Lbare}) for (rmat_t& m : *v) m=rmat_t(m.rows(), m.columns(), 0.0);
    itsFed=false;
}

std::ostream& ACBN0::Write(std::ostream& os) const
{
    const double eV=27.211386245988;
    for (const HubbardEstimate& e : Evaluate())
        os<<"[ACBN0] site "<<e.site<<" l="<<e.l<<": U="<<std::fixed<<std::setprecision(4)<<e.Ubar*eV<<" J="<<e.Jbar*eV
          <<" U_eff="<<e.Ueff()*eV<<" eV  (bare, unrenormalised: U="<<e.UbarBare*eV<<" J="<<e.JbarBare*eV<<")"
          <<"  N up/dn="<<e.chargeUpBare<<"/"<<e.chargeDnBare<<" (renormalised "<<e.chargeUp<<"/"<<e.chargeDn<<")"<<std::endl;
    return os;
}

} // namespace
