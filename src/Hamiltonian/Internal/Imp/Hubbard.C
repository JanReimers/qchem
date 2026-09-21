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
#include <typeinfo>
#include <vector>
module qchem.Hamiltonian.Internal.Hubbard;
import qchem.Energy;
import qchem.RunPolicy;                         // theRunPolicy().HubbardEigen() -- the form, read once here
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
    : itsSt(st), itsManifolds(std::move(manifolds)), itsGroup(g), itsEigenForm(theRunPolicy().HubbardEigen())
{
    if (itsManifolds.empty()) throw std::invalid_argument("Hubbard_U: no manifold -- a +U term with nothing to correct");
    itsSt->ForEachSite([this](int Z, const rvec3_t& R, bool){itsSites.push_back(R); itsSiteZ.push_back(Z);});
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
struct Selection
{
    std::vector<std::vector<size_t>> cols;
    std::vector<std::vector<Symmetry::Molecule::AoShell>> shells;
};
template <class U> Selection
ColumnsOf(const BasisSet::Orbital_1E_IBS<U>& orb, const std::vector<HubbardManifold>& manifolds,
          const std::vector<rvec3_t>& sites)
{
    // "I am built from atom-centred shells" -- an abstract face, so this is abstract->abstract.
    const auto* src=dynamic_cast<const BasisSet::AoShellSource*>(&orb);
    if (!src) throw std::runtime_error("Hubbard_U: the orbital block is not built from atom-centred shells "
                                       "(no AoShellSource face) -- a Hubbard manifold cannot be selected on it");
    const std::vector<Symmetry::Molecule::AoShell> shells=src->GetAoShells();
    Selection sel;
    for (const HubbardManifold& M : manifolds)
    {
        std::vector<size_t> c;
        std::vector<Symmetry::Molecule::AoShell> picked;
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
            picked.push_back(sh);
        }
        if (c.empty())
            throw std::runtime_error("Hubbard_U: no l="+std::to_string(M.l)+" shell on site "+std::to_string(M.site));
        sel.cols.push_back(std::move(c));
        sel.shells.push_back(std::move(picked));
    }
    return sel;
}

} // anonymous namespace

//================================================================================= ManifoldSymmetry
// A representation is CLUSTERED, not looked up: the molecular character tables are abelian-only, and a
// d shell on a cubic or trigonal site carries 3-D / 2-D irreps.  Symmetrise a generic matrix under the
// group, eigen-decompose, and a degenerate cluster IS an irrep copy; its character vector chi(g) = Tr(P D(g))
// names it -- two clusters with equal characters are the same irrep.
namespace
{
using Sig = ManifoldSymmetry::IrrepSig;

rsmat_t GenericSymmetric(size_t m)
{
    rsmat_t n(m);
    unsigned long long x=88172645463325252ULL;                       // a fixed-seed xorshift: reproducible
    auto rnd=[&]{ x^=x<<13; x^=x>>7; x^=x<<17; return double(x%1000003)/1000003.0-0.5; };
    for (size_t a=0;a<m;a++) for (size_t b=a;b<m;b++) n(a,b)=rnd();
    return n;
}
rsmat_t Symmetrise(const rsmat_t& n, const std::vector<rmat_t>& D)
{
    if (D.size()<=1) return n;
    const size_t m=n.rows();
    rmat_t acc(m,m,0.0);
    for (const rmat_t& d : D) acc += d*rmat_t(n)*blazem::trans(d);
    acc /= double(D.size());
    rsmat_t out(m);
    for (size_t a=0;a<m;a++) for (size_t b=a;b<m;b++) out(a,b)=0.5*(acc(a,b)+acc(b,a));
    return out;
}
//! Cluster the (ascending) eigenvalues by degeneracy; each cluster's members are eigen indices.
std::vector<std::vector<size_t>> ClusterDegenerate(const rvec_t& lam, double tol)
{
    std::vector<std::vector<size_t>> cl;
    for (size_t i=0;i<lam.size();i++)
    {
        if (!cl.empty() && std::abs(lam[i]-lam[cl.back().back()])<=tol*std::max(1.0,std::abs(lam[i]))) cl.back().push_back(i);
        else cl.push_back({i});
    }
    return cl;
}
rvec_t Characters(const std::vector<size_t>& members, const rmat_t& v, const std::vector<rmat_t>& D)
{
    rvec_t chi(D.size(), 0.0);
    for (size_t g=0; g<D.size(); g++)
    {
        double t=0.0;
        for (size_t i : members)
        {
            const auto vi=blazem::column(v,i);
            t += blazem::dot(vi, D[g]*vi);
        }
        chi[g]=t;
    }
    return chi;
}
bool SameSig(const rvec_t& a, const rvec_t& b, double tol=1e-6)
{
    if (a.size()!=b.size()) return false;
    for (size_t i=0;i<a.size();i++) if (std::abs(a[i]-b[i])>tol) return false;
    return true;
}
//! The irreps a representation contains, in a FIXED order (dimension, then characters lexicographically),
//! read off a generic symmetric matrix symmetrised under it -- so the table does not depend on any density.
std::vector<Sig> IrrepsOf(const std::vector<rmat_t>& D, size_t m)
{
    rsmat_t nb=Symmetrise(GenericSymmetric(m), D);
    rvec_t lam; rmat_t v; blazem::eigen(nb, lam, v);
    std::vector<Sig> out;
    for (const auto& c : ClusterDegenerate(lam, 1e-8))
    {
        Sig sg{Characters(c, v, D), c.size()};
        bool seen=false; for (const Sig& o : out) if (SameSig(o.chi, sg.chi)) seen=true;
        if (!seen) out.push_back(std::move(sg));
    }
    std::sort(out.begin(), out.end(), [](const Sig& a, const Sig& b)
    {
        if (a.dim!=b.dim) return a.dim<b.dim;
        for (size_t i=0;i<a.chi.size();i++) if (std::abs(a.chi[i]-b.chi[i])>1e-6) return a.chi[i]>b.chi[i];
        return false;
    });
    return out;
}
//! The isotypic projectors, one per irrep: A_k = Sum_g chi_k(g) D(g) is PROPORTIONAL to the projector onto
//! irrep k's isotypic component (exactly (|G|/d_k) P_k for an absolutely irreducible irrep; a complex-type
//! pair seen as one real 2-D cluster carries a different factor), so the scale is read off A itself:
//! Tr(A^2)/Tr(A) = c for A = c P.  Symmetric because the characters are real and D orthogonal.
std::vector<rmat_t> Projectors(const std::vector<Sig>& sigs, const std::vector<rmat_t>& D, size_t m)
{
    std::vector<rmat_t> P;
    for (const Sig& g : sigs)
    {
        rmat_t A(m,m,0.0);
        for (size_t o=0;o<D.size();o++) A += g.chi[o]*D[o];
        double trA=0.0, trA2=0.0;
        for (size_t a=0;a<m;a++) { trA+=A(a,a); for (size_t b=0;b<m;b++) trA2+=A(a,b)*A(b,a); }
        if (trA<=1e-12) throw std::logic_error("ManifoldSymmetry: an irrep with an empty isotypic component");
        A *= trA/trA2;
        P.push_back(std::move(A));
    }
    return P;
}
} // anonymous namespace

std::vector<rmat_t> ManifoldSymmetry::Rep(const std::vector<Symmetry::Molecule::AoShell>& shells, const std::vector<rmat3d_t>& ops)
{
    size_t m=0; for (const auto& sh : shells) m+=sh.nComponents();
    std::vector<rmat3d_t> R=ops; if (R.empty()) R.push_back(rmat3d_t(1.0,0,0, 0,1.0,0, 0,0,1.0));
    std::vector<rmat_t> D;
    for (const rmat3d_t& r : R)
    {
        rmat_t d(m,m,0.0); size_t at=0;
        for (const auto& sh : shells)
        {
            const rmat_t b=sh.rep->Rep(r); const size_t k=sh.nComponents();
            assert(sh.norm.size()==k && "ManifoldSymmetry::Rep: a shell without its per-component norms");
            // normalised rep (BuildOperationRep's convention): phi_a(R^-1 r) = Sum_b (N_a/N_b) Rep(b,a) phi_b(r)
            for (size_t bb=0;bb<k;bb++) for (size_t a=0;a<k;a++) d(at+bb,at+a)=(sh.norm[a]/sh.norm[bb])*b(bb,a);
            at+=k;
        }
        D.push_back(std::move(d));
    }
    return D;
}

ManifoldSymmetry::ManifoldSymmetry(std::vector<rmat_t> D, std::vector<rmat_t> Dgrey)
    : itsM(D.empty() ? 0 : D[0].rows()), itsD(std::move(D)), itsDgrey(std::move(Dgrey))
{
    if (itsD.empty() || itsDgrey.empty()) throw std::invalid_argument("ManifoldSymmetry: an empty group (the identity is an op too)");
    itsIrreps=IrrepsOf(itsD,     itsM);
    itsGrey  =IrrepsOf(itsDgrey, itsM);
    itsPsite =Projectors(itsIrreps, itsD,     itsM);
    itsPgrey =Projectors(itsGrey,   itsDgrey, itsM);
    // THE SLOTS, from group theory alone: dim = Tr(P_site,k P_grey,p), an integer because the two commute
    // (the site group is a subgroup of the grey one) and their product is the projector onto the intersection.
    for (size_t k=0;k<itsPsite.size();k++)
        for (size_t p=0;p<itsPgrey.size();p++)
        {
            double t=0.0;
            for (size_t a=0;a<itsM;a++) for (size_t b=0;b<itsM;b++) t+=itsPsite[k](a,b)*itsPgrey[p](b,a);
            const long d=std::lround(t);
            if (std::abs(t-double(d))>1e-6)
                throw std::logic_error("ManifoldSymmetry: Tr(P_site P_grey) is not an integer -- the site group is not a subgroup of the grey one");
            if (d>0) itsSlots.push_back({k, p, size_t(d)});
        }
    std::sort(itsSlots.begin(), itsSlots.end(), [](const Slot& a, const Slot& b)
    { return a.dim!=b.dim ? a.dim<b.dim : a.parent!=b.parent ? a.parent<b.parent : a.irrep<b.irrep; });
}

ManifoldSymmetry::Labelling ManifoldSymmetry::Label(const rsmat_t& n, rvec_t& lam, rmat_t& v) const
{
    blazem::eigen(n, lam, v);                                         // the density's OWN occupations
    const size_t m=lam.size();
    // A degenerate cluster's eigenbasis is arbitrary: rotate it to diagonalise the projectors, so that a
    // vector inside the cluster is EITHER in one slot or another (n = Sum lam v v^T is unchanged by any
    // rotation inside the cluster).  Q's spectrum separates every (site irrep, grey parent) pair.
    for (const auto& c : ClusterDegenerate(lam, 1e-7))
    {
        if (c.size()<2) continue;
        rmat_t Q(c.size(), c.size(), 0.0);
        auto add=[&](const std::vector<rmat_t>& P, double scale)
        {
            for (size_t k=0;k<P.size();k++)
            {
                const double w=scale*double(k+1);
                for (size_t i=0;i<c.size();i++)
                    for (size_t j=0;j<c.size();j++)
                        Q(i,j)+=w*blazem::dot(blazem::column(v,c[i]), P[k]*blazem::column(v,c[j]));
            }
        };
        add(itsPsite, 1.0);
        add(itsPgrey, 1.0/double(itsPgrey.size()+1));   // the grey part stays below one site step
        rsmat_t Qs(c.size());
        for (size_t i=0;i<c.size();i++) for (size_t j=i;j<c.size();j++) Qs(i,j)=0.5*(Q(i,j)+Q(j,i));
        rvec_t mu; rmat_t Y; blazem::eigen(Qs, mu, Y);
        rmat_t Vc(m, c.size());
        for (size_t i=0;i<c.size();i++) for (size_t a=0;a<m;a++) Vc(a,i)=v(a,c[i]);
        rmat_t Vr=Vc*Y;
        for (size_t i=0;i<c.size();i++) for (size_t a=0;a<m;a++) v(a,c[i])=Vr(a,i);
    }
    Labelling L;
    L.slot.assign(m, size_t(-1));
    L.parentage=rvec_t(itsSlots.size(), 1.0);
    auto dominant=[&](const std::vector<rmat_t>& P, size_t i, double& weight)
    {
        const auto vi=blazem::column(v,i);
        size_t best=0; weight=-1.0;
        for (size_t k=0;k<P.size();k++) { const double w=blazem::dot(vi, P[k]*vi); if (w>weight) { weight=w; best=k; } }
        return best;
    };
    for (size_t i=0;i<m;i++)
    {
        double wSite=0.0, wGrey=0.0;
        const size_t irr   =dominant(itsPsite, i, wSite);
        const size_t parent=dominant(itsPgrey, i, wGrey);
        L.purity=std::min(L.purity, wSite);
        // The slot (irr, parent); an impure vector whose dominant pair is not an intersection takes the
        // slot of its site irrep with the heaviest allowed parent.
        size_t slot=size_t(-1); double bestW=-1.0;
        for (size_t k=0;k<itsSlots.size();k++)
        {
            if (itsSlots[k].irrep!=irr) continue;
            if (itsSlots[k].parent==parent) { slot=k; bestW=wGrey; break; }
            const double w=blazem::dot(blazem::column(v,i), itsPgrey[itsSlots[k].parent]*blazem::column(v,i));
            if (w>bestW) { bestW=w; slot=k; }
        }
        assert(slot!=size_t(-1) && "ManifoldSymmetry::Label: every site irrep has at least one slot");
        L.slot[i]=slot;
        L.parentage[slot]=std::min(L.parentage[slot], bestW);
    }
    return L;
}

std::ostream& ManifoldSymmetry::WriteSlots(std::ostream& os, double U, const std::vector<double>& Uirrep) const
{
    for (size_t k=0;k<itsSlots.size();k++)
    {
        const Slot& S=itsSlots[k];
        os<<"  ["<<k<<"] dim "<<S.dim<<" (site irrep "<<S.irrep<<" dim "<<itsIrreps[S.irrep].dim<<")";
        if (Graded()) os<<" < grey irrep "<<S.parent<<" dim "<<itsGrey[S.parent].dim;
        os<<" U="<<std::setprecision(3)<<(Uirrep.empty() ? U : Uirrep[k])*27.211386245988<<" eV";
    }
    return os;
}

double DudarevInEigenbasis(const rvec_t& lam, const rmat_t& v, const std::vector<size_t>& slot,
                           const std::vector<double>& Uirrep, double U, rmat_t& W)
{
    const size_t m=lam.size();
    W=rmat_t(m,m,0.0);
    double EU=0.0;
    for (size_t i=0;i<m;i++)
    {
        const double Ui = (!Uirrep.empty() && i<slot.size() && slot[i]!=size_t(-1)) ? Uirrep[slot[i]] : U;
        EU+=0.5*Ui*lam[i]*(1.0-lam[i]);
        const double w=Ui*(0.5-lam[i]);
        for (size_t a=0;a<m;a++) for (size_t b=0;b<m;b++) W(a,b)+=w*v(a,i)*v(b,i);
    }
    return EU;
}

//------------------------------------------------------------------------------ the manifold's functions (cont.)
std::vector<std::vector<size_t>> Hubbard_U::Columns(const BasisSet::Orbital_1E_IBS<double>& orb) const
{return ColumnsOf<double>(orb, itsManifolds, itsSites).cols;}
std::vector<std::vector<size_t>> Hubbard_U::Columns(const BasisSet::Orbital_1E_IBS<dcmplx>& orb) const
{return ColumnsOf<dcmplx>(orb, itsManifolds, itsSites).cols;}
std::vector<std::vector<Symmetry::Molecule::AoShell>> Hubbard_U::Shells(const BasisSet::Orbital_1E_IBS<double>& orb) const
{return ColumnsOf<double>(orb, itsManifolds, itsSites).shells;}
std::vector<std::vector<Symmetry::Molecule::AoShell>> Hubbard_U::Shells(const BasisSet::Orbital_1E_IBS<dcmplx>& orb) const
{return ColumnsOf<dcmplx>(orb, itsManifolds, itsSites).shells;}

void Hubbard_U::BuildSymmetry(size_t M, const std::vector<Symmetry::Molecule::AoShell>& shells) const
{
    if (itsSym.size()<itsManifolds.size()) itsSym.resize(itsManifolds.size());
    if (itsSym[M]) return;                                            // built once
    const HubbardManifold& man=itsManifolds[M];
    auto sym=std::make_shared<const ManifoldSymmetry>(ManifoldSymmetry::Rep(shells, man.siteOps),
                                                      ManifoldSymmetry::Rep(shells, man.greyOps.empty() ? man.siteOps : man.greyOps));
    if (!man.Uirrep.empty() && man.Uirrep.size()!=sym->Slots().size())
        throw std::invalid_argument("Hubbard_U: manifold (site "+std::to_string(man.site)+", l="+std::to_string(man.l)
            +") has "+std::to_string(sym->Slots().size())+" U slots (site irreps x grey parents) but Uirrep carries "
            +std::to_string(man.Uirrep.size())+" -- see the [+U] label table");
    itsSym[M]=sym;
    // Say it once (pin 17): the table a user reads Uirrep against.  A DECLARATION, like the run banner,
    // so it goes to stdout as well as the console heartbeat.
    std::ostringstream os;
    os<<"[+U] site "<<man.site<<" l="<<man.l<<": "<<sym->Size()<<" functions, site group of "<<sym->NumOps()<<" ops (grey "
      <<sym->NumGreyOps()<<");  U slots:";
    sym->WriteSlots(os, man.U, man.Uirrep);
    std::cout<<os.str()<<std::endl;
    report::Log(os.str());
}

const std::vector<ManifoldSymmetry::Slot>& Hubbard_U::Slots(size_t M) const
{
    static const std::vector<ManifoldSymmetry::Slot> none;
    return (M<itsSym.size() && itsSym[M]) ? itsSym[M]->Slots() : none;
}

const ManifoldSymmetry::Labelling& Hubbard_U::LastLabelling(size_t M, const Spin& s) const
{
    auto it=itsChannels.find(s);
    if (it==itsChannels.end()) throw std::logic_error("Hubbard_U::LastLabelling: no such spin channel on this term");
    if (it->second.lab.size()<=M) throw std::out_of_range("Hubbard_U::LastLabelling: no such manifold, or no refresh yet");
    return it->second.lab[M];
}

const rvec_t& Hubbard_U::OccupationByLabel(size_t M, const Spin& s) const
{
    auto it=itsChannels.find(s);
    if (it==itsChannels.end()) throw std::logic_error("Hubbard_U::OccupationByLabel: no such spin channel on this term");
    if (it->second.byLabel.size()<=M) throw std::out_of_range("Hubbard_U::OccupationByLabel: no such manifold, or no refresh yet");
    return it->second.byLabel[M];
}

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
double Hubbard_U::Analyse(Channel& ch) const
{
    const rvec_t& n=ch.n;
    ch.occ.clear();     ch.occ.resize(itsManifolds.size());
    ch.byLabel.clear(); ch.byLabel.resize(itsManifolds.size());
    ch.lab.clear();     ch.lab.resize(itsManifolds.size());
    ch.W=rvec_t(n.size(), 0.0);
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
        std::vector<size_t> slot;                                  // per eigenvalue: its U slot (increment 2)
        const HubbardManifold& man=itsManifolds[M];
        if (!itsEigenForm)
        {   // CP2K parity (RunPolicy::HubbardEigen): the POPULATIONS are the "eigenvalues", the basis is the identity.
            lam=rvec_t(m); v=rmat_t(m,m,0.0);
            for (size_t a=0;a<m;a++) {lam[a]=nM(a,a); v(a,a)=1.0;}
        }
        else if (M<itsSym.size() && itsSym[M])
        {   // the density's own n, its eigenvectors NAMED by the site group (never symmetrised)
            ch.lab[M]=itsSym[M]->Label(nM, lam, v);
            slot=ch.lab[M].slot;
            rvec_t sums(itsSym[M]->Slots().size(), 0.0);
            for (size_t i=0;i<m;i++) sums[slot[i]]+=lam[i];
            ch.byLabel[M]=sums;
        }
        else blazem::eigen(nM, lam, v);                           // n = Sum_i lam_i v_i v_i^T
        ch.occ[M]=lam;
        // Dudarev in the eigenbasis (Macke eq 6), one U per slot (Uirrep) or the manifold's U.
        rmat_t W;
        EU+=DudarevInEigenbasis(lam, v, slot, man.Uirrep, man.U, W);
        for (size_t a=0;a<m;a++) for (size_t b=0;b<m;b++) ch.W[at+a*m+b]=W(a,b);
        at+=m*m;
    }
    return EU;
}

namespace
{
// THE DM BEHIND A DENSITY: the density itself when it carries a D, else the DM-backed SOURCE a mixed field
// retains (cDM_Sourced_CD, one cast away -- the XC cusp-deficit route's idiom); the returned shared_ptr keeps
// a retained source alive across the call.  Null = matrix-free with no source (the seed).
const cDM_CD* DMBehind(const cChargeDensity* cd, std::vector<std::shared_ptr<const cDM_CD>>& keep)
{
    if (!cd) return nullptr;
    if (const auto* dm=dynamic_cast<const cDM_CD*>(cd)) return dm;
    if (const auto* src=dynamic_cast<const ChargeDensity::cDM_Sourced_CD*>(cd))
        if (auto held=src->DMSource()) { keep.push_back(held); return held.get(); }
    return nullptr;
}
}

void Hubbard_U::EnsureOccupations(const cChargeDensity* cd) const
{
    if (!cd || itsFrozen) return;
    if (cd->Version()==itsOccVersion) return;
    // PER CHANNEL, because that is where the D lives on a polarized run: the polarized MIXED density
    // (PolarizedMixCD) carries no D of its own and answers no source face -- each of its channel views does
    // (the mixer seats the split D on them).  The first MnO run zeroed the occupations on every Fock build
    // because it asked the total (2026-09-20); the trace read n=0 on alternate refreshes.
    // A channel with NO source (the matrix-free seed) KEEPS whatever occupations it has: a density that
    // cannot say its D says nothing about n -- only before any occupations exist is n=0 the answer.
    std::vector<std::shared_ptr<const cDM_CD>> keep;        // retained sources stay alive across the call
    std::map<Spin,const cDM_CD*> dms;
    bool any=false;
    for (const auto& [s,ch] : itsChannels)
    {
        const cChargeDensity* chan = (s==Spin::None) ? cd : ChannelOf<dcmplx>(cd, s);
        const cDM_CD* dm=DMBehind(chan, keep);
        if (!dm && s!=Spin::None) dm=DMBehind(cd, keep);    // a spin-agnostic total under a polarized term
        if (dm) any=true;
        dms[s]=dm;
    }
    if (!any)
    {
        if (itsProj.empty() && itsProjR.empty()) return;    // before any block: nothing to size n by (seed)
        for (const auto& [s,ch] : itsChannels) if (ch.W.size()==NumCoefficients()) return;   // keep the last n
        // nothing yet: n = 0 in every channel (V = U/2 P, E_U = 0); the version is NOT stamped so the first
        // DM-backed density does the real work.
        for (auto& [s,ch] : itsChannels) { ch.n=rvec_t(NumCoefficients(),0.0); Analyse(ch); }
        return;
    }
    itsEU=0.0;
    for (auto& [s,ch] : itsChannels)
    {
        const cDM_CD* dm=dms[s];
        rvec_t n = dm ? dm->ProjectOnto(*this) : ch.n;         // a source-less channel keeps its n
        if (n.size()!=NumCoefficients()) n=rvec_t(NumCoefficients(),0.0);
        if (dm && s==Spin::None) n*=0.5;                       // the zeta=0 collapse: n_sigma = n_tot/2
        ch.n=n;
        itsEU+=Analyse(ch);
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
                if (ch.byLabel.size()>M && ch.byLabel[M].size()>0)
                {   // the orbital resolution: Sum lambda per U slot, in the [+U] label table's order
                    os<<" by slot:";
                    for (size_t k=0;k<ch.byLabel[M].size();k++) os<<" ["<<k<<"]"<<std::setprecision(3)<<ch.byLabel[M][k];
                    os<<" purity "<<std::setprecision(3)<<ch.lab[M].purity;
                    if (itsSym[M]->Graded())
                    {   // parentage: how much of each slot's eigenvectors sits in the grey parent it is named by
                        os<<" parentage";
                        for (double w : ch.lab[M].parentage) os<<" "<<std::setprecision(2)<<100*w<<"%";
                    }
                }
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
    // The site groups on the manifolds (increment 2): the shell layout is uniform across blocks, so the first
    // block's shells give every manifold's representation.
    if (bs->GetNumIBS()>0 && itsEigenForm)
    {
        auto shells = bs->GetRealIBS(0) ? Shells(*bs->GetRealIBS(0)) : Shells(*(*bs)[0]);
        for (size_t M=0;M<itsManifolds.size();M++) BuildSymmetry(M, shells[M]);
    }
}

void Hubbard_U::RefreshForDensity(const cChargeDensity* cd) const
{
    EnsureOccupations(cd);
}

//------------------------------------------------------------------------------ HubbardProjection (ACBN0)
std::vector<size_t> Hubbard_U::EquivalentManifolds(size_t M) const
{
    std::vector<size_t> eq;
    for (size_t K=0;K<itsManifolds.size();K++)
        if (itsManifolds[K].l==itsManifolds[M].l && itsSiteZ[itsManifolds[K].site]==itsSiteZ[itsManifolds[M].site]) eq.push_back(K);
    return eq;
}
template <class U> static std::vector<mat_t<U>> LowdinOf(const LowdinProjector<U>& P, const mat_t<U>& C)
{
    std::vector<mat_t<U>> out;
    for (size_t M=0;M<P.NumManifolds();M++) out.push_back(mat_t<U>(blazem::ctrans(P.T(M))*C));   // T^dagger C
    return out;
}
std::vector<mat_t<double>> Hubbard_U::LowdinCoefficients(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& orb, const mat_t<double>& C) const
{return LowdinOf<double>(Projector<double>(orb), C);}
std::vector<mat_t<dcmplx>> Hubbard_U::LowdinCoefficients(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& orb, const mat_t<dcmplx>& C) const
{return LowdinOf<dcmplx>(Projector<dcmplx>(orb), C);}

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
    os<<"Hubbard +U (Lowdin):";
    for (size_t M=0;M<itsManifolds.size();M++)
    {
        const HubbardManifold& man=itsManifolds[M];
        os<<" site "<<man.site<<" l="<<man.l;
        if (man.Uirrep.empty()) os<<" U="<<man.U*27.211386245988<<" eV (shell-averaged)";
        else { os<<" Uirrep="; for (size_t k=0;k<man.Uirrep.size();k++) os<<(k?",":"")<<man.Uirrep[k]*27.211386245988; os<<" eV"; }
        if (M<itsSym.size() && itsSym[M]) { os<<" slots:"; itsSym[M]->WriteSlots(os, man.U, man.Uirrep); }
        os<<";";
    }
    os<<(itsEigenForm ? "  [eigenvalue form]" : "  [DIAGONAL populations, CP2K's form]");
    return os<<std::endl;
}

} // namespace
