// File: Hamiltonian/Internal/Imp/Hubbard.C  DFT+U -- implementation.  See the module interface for the design.
module;
#include <algorithm>   // std::max_element
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
module qchem.Hamiltonian.Internal.Hubbard;
import qchem.Energy;
import qchem.ChargeDensity;                     // cDM_CD, ChannelOf, tDM_Sourced_CD (the DM-backed source)
import qchem.BasisSet.AoShellSource;            // the shell layout a manifold is selected from
import qchem.Symmetry.Molecule.OperationRep;    // AoShell (+ ShellRep::L / Monomials)
import qchem.Reporting;                         // the term's own [+U] line (pin 17: contemporaneous)
import qchem.Blaze;
import qchem.Math;
import qchem.Vector3D;                          // norm(rvec3_t) -- the site match
import qchem.Symmetry.Irrep;                    // Spin, SpinIrreps

namespace qchem::Hamiltonian
{

namespace
{
template <class U> double RealPart(const U& x)
{
    if constexpr (std::is_floating_point_v<U>) return x; else return std::real(x);
}
}

//============================================================================================= LowdinProjector

template <class U> LowdinProjector<U>::LowdinProjector(const hmat_t<U>& S, const std::vector<std::vector<size_t>>& columns)
    : itsN(0)
{
    // S^{1/2} = V diag(sqrt w) V^dagger from the block's overlap.  A near-null direction of S (a
    // rank-deficient diffuse span) gets sqrt(w)->0 and simply does not project, which is the honest
    // answer for a function the block cannot resolve.
    const size_t n=S.rows();
    rvec_t w; mat_t<U> V;
    blazem::eigen(S, w, V);
    mat_t<U> half(n,n);
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<n;j++)
        {
            U s=U(0);
            for (size_t k=0;k<n;k++) s+=V(i,k)*std::sqrt(std::max(w[k],0.0))*blazem::conjs(V(j,k));
            half(i,j)=s;
        }
    for (const auto& cols : columns)
    {
        mat_t<U> T(n, cols.size());
        for (size_t c=0;c<cols.size();c++)
        {
            if (cols[c]>=n) throw std::out_of_range("LowdinProjector: manifold column outside the block");
            for (size_t i=0;i<n;i++) T(i,c)=half(i,cols[c]);
        }
        itsN+=cols.size()*cols.size();
        itsT.push_back(std::move(T));
    }
}

template <class U> rvec_t LowdinProjector<U>::Forward(const hmat_t<U>& D) const
{
    rvec_t out(itsN, 0.0);
    size_t at=0;
    for (const mat_t<U>& T : itsT)
    {
        const size_t m=T.columns();
        mat_t<U> DT=D*T;                                          // n x m
        for (size_t a=0;a<m;a++)
            for (size_t b=0;b<m;b++)
            {
                U s=U(0);
                for (size_t i=0;i<T.rows();i++) s+=blazem::conjs(T(i,a))*DT(i,b);   // (T^dagger D T)_{ab}
                out[at+a*m+b]=RealPart(s);   // Hermitian after the k-sum; the imaginary part is checked by the term
            }
        at+=m*m;
    }
    return out;
}

template <class U> hmat_t<U> LowdinProjector<U>::Adjoint(const rvec_t& W) const
{
    assert(W.size()==itsN && "LowdinProjector::Adjoint: one entry per occupation-matrix element");
    const size_t n=itsT.empty() ? 0 : itsT[0].rows();
    mat_t<U> V(n,n,U(0));
    size_t at=0;
    for (const mat_t<U>& T : itsT)
    {
        const size_t m=T.columns();
        mat_t<U> Wm(m,m);
        for (size_t a=0;a<m;a++) for (size_t b=0;b<m;b++) Wm(a,b)=U(W[at+a*m+b]);
        V+=T*Wm*blazem::ctrans(T);                                // T W T^dagger
        at+=m*m;
    }
    hmat_t<U> H(n);
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++) H(i,j)=U(0.5)*(V(i,j)+blazem::conjs(V(j,i)));   // symmetrise the roundoff away
    return H;
}

template <class U> double LowdinProjector<U>::Integrate(const rvec_t& f) const
{
    double tr=0.0; size_t at=0;
    for (const mat_t<U>& T : itsT)
    {
        const size_t m=T.columns();
        for (size_t a=0;a<m;a++) tr+=f[at+a*m+a];
        at+=m*m;
    }
    return tr;
}

template class LowdinProjector<double>;
template class LowdinProjector<dcmplx>;

//================================================================================================= Hubbard_U

Hubbard_U::Hubbard_U(const std::shared_ptr<const Structure>& st, std::vector<HubbardManifold> manifolds, SpinGroup g)
    : itsSt(st), itsManifolds(std::move(manifolds)), itsGroup(g)
{
    if (itsManifolds.empty()) throw std::invalid_argument("Hubbard_U: no manifold -- a +U term with nothing to correct");
    itsSt->ForEachSite([this](int, const rvec3_t& R, bool){itsSites.push_back(R);});
    for (const HubbardManifold& M : itsManifolds)
    {
        if (M.site>=itsSites.size())
            throw std::out_of_range("Hubbard_U: manifold site "+std::to_string(M.site)+" -- the cell has "
                                    +std::to_string(itsSites.size())+" sites");
        if (M.l<0) throw std::invalid_argument("Hubbard_U: negative l");
    }
    for (Spin s : SpinIrreps(g)) itsChannels[s];   // one slot per spin irrep, from construction (R1.0h)
    // itsNCoeff is fixed by the FIRST block seen (every block of one run has the same shell layout).
}

Hubbard_U::~Hubbard_U() = default;

//------------------------------------------------------------------------------ the manifold's functions
namespace
{
template <class U> std::vector<std::vector<size_t>>
ColumnsOf(const BasisSet::Orbital_1E_IBS<U>& orb, const std::vector<HubbardManifold>& manifolds,
          const std::vector<rvec3_t>& sites)
{
    // "I am built from atom-centred shells" -- an abstract face, so this is abstract->abstract.
    const auto* src=dynamic_cast<const BasisSet::AoShellSource*>(&orb);
    if (!src) throw std::runtime_error("Hubbard_U: the orbital block is not built from atom-centred shells "
                                       "(no AoShellSource face) -- a Hubbard manifold cannot be selected on it");
    const std::vector<Symmetry::Molecule::AoShell> shells=src->GetAoShells();
    std::vector<std::vector<size_t>> cols;
    for (const HubbardManifold& M : manifolds)
    {
        std::vector<size_t> c;
        for (const auto& sh : shells)
        {
            if (norm(sh.center-sites[M.site])>1e-8 || sh.rep->L()!=M.l) continue;
            // CP2K PARITY (dft_plus_u.F): the manifold is EVERY shell of angular momentum l on the site --
            // all nsb contractions, an nsb(2l+1) x nsb(2l+1) occupation matrix -- not one "3d" picked out.
            if (M.l>=2 && !sh.rep->Monomials().empty())
                throw std::runtime_error("Hubbard_U: the l="+std::to_string(M.l)+" shell at site "
                    +std::to_string(M.site)+" is CARTESIAN ("+std::to_string(sh.nComponents())
                    +" components, the s-contaminant among them) -- a Hubbard manifold needs the 2l+1 "
                     "real harmonics; run the spherical view (GPW_SPHERICAL=1)");
            for (size_t k=0;k<sh.nComponents();k++) c.push_back(sh.offset+k);
        }
        if (c.empty())
            throw std::runtime_error("Hubbard_U: no l="+std::to_string(M.l)+" shell on site "+std::to_string(M.site));
        cols.push_back(std::move(c));
    }
    return cols;
}
}

std::vector<std::vector<size_t>> Hubbard_U::Columns(const BasisSet::Orbital_1E_IBS<double>& orb) const
{return ColumnsOf<double>(orb, itsManifolds, itsSites);}
std::vector<std::vector<size_t>> Hubbard_U::Columns(const BasisSet::Orbital_1E_IBS<dcmplx>& orb) const
{return ColumnsOf<dcmplx>(orb, itsManifolds, itsSites);}

template <> SymMap<LowdinProjector<double>>& Hubbard_U::Cache<double>() const {return itsProjR;}
template <> SymMap<LowdinProjector<dcmplx>>& Hubbard_U::Cache<dcmplx>() const {return itsProj;}

template <class U> const LowdinProjector<U>& Hubbard_U::Projector(const BasisSet::Orbital_DFT_IBS<U,dcmplx>& orb) const
{
    SymMap<LowdinProjector<U>>& cache=Cache<U>();
    const sym_t& id=orb.GetSymt();
    auto it=cache.find(id);
    if (it!=cache.end()) return it->second;
    LowdinProjector<U> p(orb.Overlap(), Columns(orb));
    if (itsNCoeff==0) itsNCoeff=p.NumCoefficients();
    else if (itsNCoeff!=p.NumCoefficients())
        throw std::logic_error("Hubbard_U: blocks disagree on the manifold size -- the shell layout is not uniform across k");
    return cache.emplace(id, std::move(p)).first->second;
}

const qcMesh::MatrixForward<double>& Hubbard_U::Forward(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& orb) const
{return Projector<double>(orb);}
const qcMesh::MatrixForward<dcmplx>& Hubbard_U::Forward(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& orb) const
{return Projector<dcmplx>(orb);}

//------------------------------------------------------------------------------------ the occupations
double Hubbard_U::Analyse(const rvec_t& n, std::vector<rvec_t>& occ, rvec_t& W) const
{
    occ.clear(); occ.resize(itsManifolds.size());
    W=rvec_t(n.size(), 0.0);
    double EU=0.0;
    size_t at=0;
    // The manifold sizes are the projectors' -- read them off the first cached block.
    std::vector<size_t> sizes;
    if (!itsProj.empty())  for (size_t M=0;M<itsManifolds.size();M++) sizes.push_back(itsProj .begin()->second.Size(M));
    else if (!itsProjR.empty()) for (size_t M=0;M<itsManifolds.size();M++) sizes.push_back(itsProjR.begin()->second.Size(M));
    else throw std::logic_error("Hubbard_U::Analyse before any block was seen");
    for (size_t M=0;M<itsManifolds.size();M++)
    {
        const size_t m=sizes[M];
        rsmat_t nM(m);
        for (size_t a=0;a<m;a++) for (size_t b=a;b<m;b++) nM(a,b)=0.5*(n[at+a*m+b]+n[at+b*m+a]);
        rvec_t lam; rmat_t v;
        blazem::eigen(nM, lam, v);                                // n = Sum_i lam_i v_i v_i^T
        occ[M]=lam;
        const double U=itsManifolds[M].U;
        // Dudarev in the eigenbasis (Macke eq 6 with U_i == U): E = Sum U/2 lam(1-lam),
        // W = Sum U (1/2 - lam) v v^T -- rotated back to the manifold's own basis here.
        for (size_t i=0;i<m;i++)
        {
            EU+=0.5*U*lam[i]*(1.0-lam[i]);
            const double w=U*(0.5-lam[i]);
            for (size_t a=0;a<m;a++) for (size_t b=0;b<m;b++) W[at+a*m+b]+=w*v(a,i)*v(b,i);
        }
        at+=m*m;
    }
    return EU;
}

void Hubbard_U::EnsureOccupations(const cChargeDensity* cd) const
{
    if (!cd || itsFrozen) return;
    if (cd->Version()==itsOccVersion) return;
    // The DM-backed density: cd itself, or the source a mixed density retains; a matrix-free seed has none.
    const cDM_CD* dm=dynamic_cast<const cDM_CD*>(cd);
    std::shared_ptr<const cDM_CD> held;
    if (!dm)
        if (auto* src=dynamic_cast<const ChargeDensity::cDM_Sourced_CD*>(cd)) { held=src->DMSource(); dm=held.get(); }
    itsEU=0.0;
    // A matrix-free density (the SAD seed) before any block was seen: there is no projector to size the
    // occupations by, and nothing to occupy -- E_U = 0, the potential stays zero (MakeMatrixT), and the
    // version is NOT stamped so the first DM-backed density does the real work.
    if (!dm && itsProj.empty() && itsProjR.empty()) return;
    for (auto& [s,ch] : itsChannels)
    {
        rvec_t n(NumCoefficients(), 0.0);
        if (dm)
        {
            const cChargeDensity* chan=ChannelOf<dcmplx>(dm, s);
            const cDM_CD* cdm=chan ? dynamic_cast<const cDM_CD*>(chan) : nullptr;
            if (cdm) n=cdm->ProjectOnto(*this);
            else if (s==Spin::None || !chan) n=dm->ProjectOnto(*this);
            if (s==Spin::None) n*=0.5;                   // the zeta=0 collapse: n_sigma = n_tot/2, both channels
        }
        ch.n=n;
        itsEU+=Analyse(ch.n, ch.occ, ch.W);
    }
    if (itsGroup==SpinGroup::UnPolarized) itsEU*=2.0;    // both (identical) channels
    itsOccVersion=cd->Version();
    // The term reports at its own activity (pin 17).
    {
        std::ostringstream os;
        os<<"[+U]";
        for (const auto& [s,ch] : itsChannels)
            for (size_t M=0;M<itsManifolds.size();M++)
            {
                double tr=0.0; for (double l : ch.occ[M]) tr+=l;
                double mx=0.0; for (double l : ch.occ[M]) mx=std::max(mx,l);
                os<<"  site "<<itsManifolds[M].site<<" l="<<itsManifolds[M].l
                  <<(s==Spin::Up ? " up" : s==Spin::Down ? " dn" : "")<<": N="<<std::fixed<<std::setprecision(4)<<tr
                  <<" max(lam)="<<mx;
            }
        os<<"  E_U="<<std::setprecision(8)<<itsEU;
        // A heartbeat on the attached console; QCHEM_U_TRACE=1 puts the same line on stdout (gtest runs
        // attach no console, and the MnO gate is read by eye against CP2K's per-step occupation print).
        static const bool trace = std::getenv("QCHEM_U_TRACE")!=nullptr;
        if (trace) std::cout<<os.str()<<std::endl; else report::Log(os.str());
    }
}

void Hubbard_U::PrepareSlots(const cbs_t* bs) const
{
    cDynamic_HT_Imp::PrepareSlots(bs);
    Dynamic_HT_RealBlock_Imp::PrepareSlots(bs);
    // The mixed-aware walk (RealComplexPlan step 3c-2): a real TRIM block inside the complex-faced set
    // answers GetRealIBS; every other block is a Bloch block.  The pairs are geometry-fixed, so build them all now.
    for (size_t i=0;i<bs->GetNumIBS();i++)
        if (const auto* r=bs->GetRealIBS(i)) Projector<double>(dynamic_cast<const BasisSet::Orbital_DFT_IBS<double,dcmplx>&>(*r));
        else                                 Projector<dcmplx>(dynamic_cast<const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>&>(*(*bs)[i]));
}

void Hubbard_U::RefreshForDensity(const cChargeDensity* cd) const
{
    EnsureOccupations(cd);
}

const rvec_t& Hubbard_U::Occupations(size_t M, const Spin& s) const
{
    auto it=itsChannels.find(s);
    if (it==itsChannels.end()) throw std::logic_error("Hubbard_U::Occupations: no such spin channel on this term");
    if (it->second.occ.size()<=M) throw std::out_of_range("Hubbard_U::Occupations: no such manifold, or no refresh yet");
    return it->second.occ[M];
}

//------------------------------------------------------------------------------------------ the matrix
template <class U> hmat_t<U> Hubbard_U::MakeMatrixT(const tobs_t<U>* bs, const Spin& s, const cChargeDensity* cd) const
{
    newCD(cd);
    EnsureOccupations(cd);                                       // ordinarily a lookup: the phase warmed it
    const auto& orb=dynamic_cast<const BasisSet::Orbital_DFT_IBS<U,dcmplx>&>(*bs);
    const Spin chan = (itsGroup==SpinGroup::UnPolarized) ? Spin::None : s;
    auto it=itsChannels.find(chan);
    if (it==itsChannels.end())
        throw std::logic_error("Hubbard_U: asked for a spin block this term was not built for (SpinIrreps(g))");
    const Channel& ch=it->second;
    if (ch.W.size()==0 || ch.W.size()!=NumCoefficients()) return blazem::zeroH<U>(bs->GetNumFunctions());   // before any density (both 0 on the seed)
    return Projector<U>(orb).Adjoint(ch.W);                       // V_k = T_k W T_k^dagger
}
chmat_t Hubbard_U::MakeMatrix (const cobs_t* bs, const Spin& s, const cChargeDensity* cd) const {return MakeMatrixT<dcmplx>(bs,s,cd);}
rsmat_t Hubbard_U::MakeMatrixR(const robs_t* bs, const Spin& s, const cChargeDensity* cd) const {return MakeMatrixT<double>(bs,s,cd);}

void Hubbard_U::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    EnsureOccupations(cd);                                        // the energy pass's density (rho_out)
    te.Add("E_U", itsEU, EnergyRole::Potential);
}

std::ostream& Hubbard_U::Write(std::ostream& os) const
{
    os<<"Hubbard +U (Lowdin, shell-averaged):";
    for (const HubbardManifold& M : itsManifolds)
        os<<" site "<<M.site<<" l="<<M.l<<" U="<<M.U*27.211386245988<<" eV;";
    return os<<std::endl;
}

} // namespace
