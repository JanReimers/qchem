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
#include <mutex>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <vector>
module qchem.Hamiltonian.Internal.Hubbard;
import qchem.Diagnostics;                     // the ONE diagnostics registry (D-ENV)
import qchem.Energy;
import qchem.RunPolicy;                         // theRunPolicy().HubbardEigen() -- the form, read once here
import qchem.ChargeDensity;                     // cDM_CD, ChannelOf, tDM_Sourced_CD (the DM-backed source)
import qchem.BasisSet.AoShellSource;
import qchem.BasisSet.BareCoulombSource;        // the bare integrals a manifold is built over (ACBN0)            // the shell layout a manifold is selected from
import qchem.Symmetry.Molecule.OperationRep;    // AoShell (+ ShellRep::L / Monomials)
import qchem.Symmetry.Molecule.SphericalRep;    // SphericalShellRep -- the site group on a pseudo-atomic manifold's harmonics
import qchem.Math.Angular;                      // Math::SphericalShell (the set's m convention)
import qchem.BasisSet.AtomCentredSet;           // the pseudo-atomic projector functions + the block's overlaps with them
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
    for (size_t M=0;M<columns.size();M++)
    {
        const std::vector<size_t>& cols=columns[M];
        for (size_t c : cols) if (c>=n) throw std::out_of_range("LowdinProjector: manifold column outside the block");
        const size_t m=cols.size();
        // CP2K's LOWDIN: T = S^{1/2}[:,M] (the whole block symmetrically orthogonalised, then the manifold's columns read
        // off), so T^dagger T = I.  Q is the SELECTOR -- an orbital's raw coefficients on those columns (ACBN0 eq 9).
        mat_t<U> T(n, m), Q(m, n, U(0));
        for (size_t a=0;a<m;a++) for (size_t i=0;i<n;i++) T(i,a)=half(i,cols[a]);
        for (size_t a=0;a<m;a++) Q(a,cols[a])=U(1);
        itsN+=m*m;
        itsT.push_back(std::move(T));
        itsQ.push_back(std::move(Q));
        itsCols.push_back(cols);
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

// THE ATOMIC-SET PAIR (HubbardProjectorPlan D4).  H^{-1/2} of a small Hermitian matrix by eigen-decomposition;
// a (near-)singular ortho group is two manifolds spanning the same functions -- a broken input, thrown.
template <class U> static mat_t<U> InverseSqrt(const hmat_t<U>& O, const char* what)
{
    const size_t m=O.rows();
    rvec_t w; mat_t<U> V; blazem::eigen(O, w, V);
    if (w[0]<=1e-8*w[m-1]) throw std::runtime_error(what);
    mat_t<U> R(m,m);
    for (size_t a=0;a<m;a++) for (size_t b=0;b<m;b++)
    { U t=U(0); for (size_t k=0;k<m;k++) t+=V(a,k)*(1.0/std::sqrt(w[k]))*blazem::conjs(V(b,k)); R(a,b)=t; }
    return R;
}

template <class U> LowdinProjector<U> LowdinProjector<U>::Atomic(const mat_t<U>& B, const hmat_t<U>& O,
                                                              const std::vector<size_t>& sizes, const std::vector<bool>& ortho)
{
    const size_t n=B.rows();
    size_t mtot=0; for (size_t m : sizes) mtot+=m;
    if (mtot!=B.columns() || O.rows()!=B.columns())
        throw std::invalid_argument("LowdinProjector::Atomic: the manifold sizes must partition B's columns, and O must be square over them");
    if (!ortho.empty() && ortho.size()!=sizes.size())
        throw std::invalid_argument("LowdinProjector::Atomic: one ortho flag per manifold (or none)");
    std::vector<size_t> off(sizes.size(),0);
    for (size_t M=1;M<sizes.size();M++) off[M]=off[M-1]+sizes[M-1];

    LowdinProjector P;
    for (size_t M=0;M<sizes.size();M++)                         // atomic: T_M = B_M
    {
        mat_t<U> T(n, sizes[M]);
        for (size_t i=0;i<n;i++) for (size_t a=0;a<sizes[M];a++) T(i,a)=B(i,off[M]+a);
        P.itsT.push_back(std::move(T));
        P.itsCols.emplace_back();
        P.itsN+=sizes[M]*sizes[M];
    }
    std::vector<size_t> members;                                // ortho-atomic: T = B_o O~^{-1/2}
    for (size_t M=0;M<sizes.size();M++) if (M<ortho.size() && ortho[M]) members.push_back(M);
    if (!members.empty())
    {
        std::vector<size_t> c;                                  // the members' columns of B / O
        for (size_t M : members) for (size_t a=0;a<sizes[M];a++) c.push_back(off[M]+a);
        const size_t mo=c.size();
        hmat_t<U> Ot(mo);
        for (size_t a=0;a<mo;a++) for (size_t b=a;b<mo;b++) Ot(a,b)=O(c[a],c[b]);
        const mat_t<U> Oih=InverseSqrt<U>(Ot, "LowdinProjector::Atomic: the ortho-atomic set is linearly dependent (two manifolds span the same functions?)");
        mat_t<U> Bo(n,mo);
        for (size_t i=0;i<n;i++) for (size_t a=0;a<mo;a++) Bo(i,a)=B(i,c[a]);
        const mat_t<U> To=Bo*Oih;
        size_t at=0;
        for (size_t M : members)
        {
            for (size_t i=0;i<n;i++) for (size_t a=0;a<sizes[M];a++) P.itsT[M](i,a)=To(i,at+a);
            at+=sizes[M];
        }
    }
    for (const mat_t<U>& T : P.itsT) P.itsQ.push_back(blazem::ctrans(T));
    return P;
}

template <class U> hmat_t<U> LowdinProjector<U>::Represented(const mat_t<U>& B, const hmat_t<U>& S)
{
    const size_t n=S.rows(), m=B.columns();
    rvec_t w; mat_t<U> V; blazem::eigen(S, w, V);
    const double floor=1e-12*w[n-1];
    mat_t<U> X=blazem::ctrans(V)*B;                             // B in S's eigenbasis
    hmat_t<U> R(m);
    for (size_t a=0;a<m;a++) for (size_t b=a;b<m;b++)
    {
        U t=U(0);
        for (size_t k=0;k<n;k++) if (w[k]>floor) t+=blazem::conjs(X(k,a))*X(k,b)/w[k];
        R(a,b)= a==b ? U(RealPart(t)) : t;   // a Hermitian diagonal is real
    }
    return R;
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
    // THE PSEUDO-ATOMIC ROUTE (HubbardProjectorPlan D6): one projector kind per run -- the two pairs are built
    // from different inputs (the block's own columns vs its overlaps with a set), so a mixed run would need a
    // joint ortho-atomic metric across both, which no physics asks for.
    size_t nAtomic=0;
    for (const HubbardManifold& M : itsManifolds) if (M.projector==HubbardProjectorKind::PseudoAtomic) nAtomic++;
    if (nAtomic>0)
    {
        if (nAtomic!=itsManifolds.size())
            throw std::invalid_argument("Hubbard_U: every manifold of a run must use the same projector kind (PseudoAtomic or SiteShells)");
        auto set=std::make_shared<BasisSet::GaussianAtomCentredSet>();
        for (const HubbardManifold& M : itsManifolds)
        {
            if (M.atomicFunction.empty())
                throw std::invalid_argument("Hubbard_U: pseudo-atomic manifold (site "+std::to_string(M.site)+", l="+std::to_string(M.l)
                                            +") has no atomicFunction -- the facade fills it from the pseudo-atom");
            set->Add(itsSites[M.site], M.l, M.atomicFunction);
        }
        itsSet=set;
    }
    else
        for (const HubbardManifold& M : itsManifolds)
            if (M.orthoAtomic)
                throw std::invalid_argument("Hubbard_U: orthoAtomic on a site-shell (LOWDIN) manifold (site "+std::to_string(M.site)
                                            +") -- LOWDIN is orthogonal over the whole block already; ortho-atomic is a pseudo-atomic option");
    for (Spin s : SpinIrreps(g)) itsChannels[s];   // one slot per spin irrep, from construction (R1.0h)
    // itsNCoeff is fixed by the FIRST block seen (every block of one run has the same shell layout).
}

Hubbard_U::~Hubbard_U() = default;

//------------------------------------------------------------------------------ the manifold's functions
namespace
{
template <class U> Hubbard_U::Selection
SelectOn(const BasisSet::Orbital_1E_IBS<U>& orb, const std::vector<HubbardManifold>& manifolds,
         const std::vector<rvec3_t>& sites)
{
    // "I am built from atom-centred shells" -- an abstract face, so this is abstract->abstract.
    const auto* src=dynamic_cast<const BasisSet::AoShellSource*>(&orb);
    if (!src) throw std::runtime_error("Hubbard_U: the orbital block is not built from atom-centred shells "
                                       "(no AoShellSource face) -- a Hubbard manifold cannot be selected on it");
    const std::vector<Symmetry::Molecule::AoShell> shells=src->GetAoShells();
    Hubbard_U::Selection sel;
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

//------------------------------------------------------------------------------ the pseudo-atomic route (D6)
namespace
{
// int over the unit sphere of x^a y^b z^c: 0 unless all even, else 2 G((a+1)/2) G((b+1)/2) G((c+1)/2) / G((a+b+c+3)/2).
double SphereMoment(int a, int b, int c)
{
    if (a%2 || b%2 || c%2) return 0.0;
    return 2.0*std::tgamma(0.5*(a+1))*std::tgamma(0.5*(b+1))*std::tgamma(0.5*(c+1))/std::tgamma(0.5*(a+b+c+3));
}
}

// ONE synthetic shell per manifold: the set's functions are SphericalShell(l)[l+m] x (positive unit-sphere scale)
// x R(r) -- an AoShell over that raw harmonic with norm[m] = 1/|raw_m| describes exactly them (the radial is
// common to all m and cancels in ManifoldSymmetry::Rep's N_a/N_b), so the site group and the slots are built by
// the SAME code a basis shell uses.
std::vector<Symmetry::Molecule::AoShell> Hubbard_U::SetShells(size_t M) const
{
    const int l=itsSet->L(M);
    Symmetry::Molecule::AoShell sh;
    sh.shellType=l; sh.center=itsSet->Centre(M); sh.offset=0;
    const auto harm=Math::SphericalShell(l);
    sh.rep=std::make_shared<Symmetry::Molecule::SphericalShellRep>(harm);
    sh.norm=rvec_t(size_t(2*l+1));
    for (size_t c=0;c<sh.norm.size();c++)
    {
        double n2=0.0;
        for (const auto& a : harm[c]) for (const auto& b : harm[c]) n2+=a.c*b.c*SphereMoment(a.p.n+b.p.n, a.p.l+b.p.l, a.p.m+b.p.m);
        sh.norm[c]=1.0/std::sqrt(n2);
    }
    return {sh};
}

// D4 + D6: B and O from the block (AtomCentredOverlapSource, abstract->abstract), the pair from them; the
// captured norm B^dag S^-1 B / O is said once -- how much of each projector the orbital basis spans, the number
// that made the site-shell route basis-dependent.
template <class U> LowdinProjector<U> Hubbard_U::AtomicProjector(const BasisSet::Orbital_DFT_IBS<U,dcmplx>& orb) const
{
    const auto* src=dynamic_cast<const BasisSet::AtomCentredOverlapSource<U>*>(&orb);
    if (!src) throw std::runtime_error("Hubbard_U: a pseudo-atomic +U projector needs the block's overlap with an atom-centred set "
                                       "(no AtomCentredOverlapSource face on this orbital block)");
    const mat_t<U>  B=src->MakeAtomCentredOverlap(*itsSet);
    const hmat_t<U> O=src->MakeSetOverlap(*itsSet);
    std::vector<size_t> sizes; std::vector<bool> ortho;
    for (size_t M=0;M<itsManifolds.size();M++) { sizes.push_back(itsSet->Size(M)); ortho.push_back(itsManifolds[M].orthoAtomic); }
    if (!itsCapturedSaid)
    {
        itsCapturedSaid=true;
        const hmat_t<U> R=LowdinProjector<U>::Represented(B, orb.Overlap());
        std::ostringstream os;
        os<<"[+U projector] pseudo-atomic (independent of the orbital basis); captured by the block, min over m of B'S^-1B/O:";
        for (size_t M=0;M<itsManifolds.size();M++)
        {
            double mn=1e300;
            for (size_t a=itsSet->Offset(M); a<itsSet->Offset(M)+itsSet->Size(M); a++) mn=std::min(mn, RealPart(R(a,a))/RealPart(O(a,a)));
            os<<"  site "<<itsManifolds[M].site<<" l="<<itsManifolds[M].l<<(itsManifolds[M].orthoAtomic ? " (ortho)" : "")
              <<" "<<std::fixed<<std::setprecision(4)<<mn;
        }
        std::cout<<os.str()<<std::endl;
        report::Log(os.str());
    }
    return LowdinProjector<U>::Atomic(B, O, sizes, ortho);
}

//------------------------------------------------------------------------------ the manifold's functions (cont.)
Hubbard_U::Selection Hubbard_U::Select(const BasisSet::Orbital_1E_IBS<double>& orb) const {return SelectOn<double>(orb, itsManifolds, itsSites);}
Hubbard_U::Selection Hubbard_U::Select(const BasisSet::Orbital_1E_IBS<dcmplx>& orb) const {return SelectOn<dcmplx>(orb, itsManifolds, itsSites);}

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
    // The functions' names, for a report that must say what a slot IS: Cartesian monomials as written, real solid
    // harmonics in SphericalShell's m = -l..l order.
    if (itsComponents.size()<itsManifolds.size()) itsComponents.resize(itsManifolds.size());
    static const std::vector<std::vector<std::string>> harmonic={{"s"},{"y","z","x"},{"xy","yz","z2","xz","x2-y2"},
        {"y(3x2-y2)","xyz","yz2","z3","xz2","z(x2-y2)","x(x2-3y2)"}};
    for (const auto& sh : shells)
    {
        const auto monos=sh.rep->Monomials();
        const int l=sh.rep->L();
        for (size_t c=0;c<sh.nComponents();c++)
        {
            std::string name;
            if (!monos.empty()) { for (int i=0;i<monos[c].n;i++) name+="x"; for (int i=0;i<monos[c].l;i++) name+="y"; for (int i=0;i<monos[c].m;i++) name+="z"; if (name.empty()) name="s"; }
            else if (l<int(harmonic.size())) name=harmonic[l][c];
            else name="m"+std::to_string(int(c)-l);
            itsComponents[M].push_back(name);
        }
    }
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
    // LOCKED: normally the refresh's occupation pass builds every block's projector first, but the SEED iteration has no D,
    // so its (concurrent) Fock blocks are the first to ask.  std::map nodes are address-stable: the reference outlives the lock.
    std::lock_guard<std::mutex> lock(itsProjLock.m);
    SymMap<LowdinProjector<U>>& cache=Cache<U>();
    const sym_t& id=orb.GetSymt();
    auto it=cache.find(id);
    if (it!=cache.end()) return it->second;
    if (PseudoAtomic())
    {
        LowdinProjector<U> p=AtomicProjector<U>(orb);
        if (itsNCoeff==0) itsNCoeff=p.NumCoefficients();
        return cache.emplace(id, std::move(p)).first->second;
    }
    const Selection sel=Select(orb);
    LowdinProjector<U> p(orb.Overlap(), sel.cols);
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
        if (man.alpha!=0.0)
        {   // the STATIC shift alpha*P (QE's Hubbard_alpha): W += alpha*1, E += alpha*Tr n -- the LR-cDFT perturbation
            for (size_t a=0;a<m;a++) {W(a,a)+=man.alpha; EU+=man.alpha*nM(a,a);}
        }
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

// PER CHANNEL, because that is where the D lives on a polarized run: the polarized MIXED density
// (PolarizedMixCD) carries no D of its own and answers no source face -- each of its channel views does
// (the mixer seats the split D on them).  The first MnO run zeroed the occupations on every Fock build
// because it asked the total (2026-09-20); the trace read n=0 on alternate refreshes.
// An entry is EMPTY for a channel with NO source (the matrix-free seed): a density that cannot say its D says
// nothing about n.
std::map<Spin,rvec_t> Hubbard_U::DensityOccupations(const cChargeDensity* cd) const
{
    std::vector<std::shared_ptr<const cDM_CD>> keep;        // retained sources stay alive across the call
    std::map<Spin,rvec_t> out;
    for (const auto& [s,ch] : itsChannels)
    {
        const cChargeDensity* chan = (s==Spin::None) ? cd : ChannelOf<dcmplx>(cd, s);
        const cDM_CD* dm=DMBehind(chan, keep);
        if (!dm && s!=Spin::None) dm=DMBehind(cd, keep);    // a spin-agnostic total under a polarized term
        rvec_t n;
        if (dm)
        {
            n=dm->ProjectOnto(*this);
            if (n.size()!=NumCoefficients()) n=rvec_t(NumCoefficients(),0.0);
            if (s==Spin::None) n*=0.5;                         // the zeta=0 collapse: n_sigma = n_tot/2
        }
        out[s]=std::move(n);
    }
    return out;
}

void Hubbard_U::EnsureOccupations(const cChargeDensity* cd) const
{
    if (!cd || itsFrozen) return;
    if (cd->Version()==itsOccVersion) return;
    // A channel with NO source (the matrix-free seed) KEEPS whatever occupations it has -- only before any
    // occupations exist is n=0 the answer.
    const std::map<Spin,rvec_t> dn=DensityOccupations(cd);
    bool any=false;
    for (const auto& [s,n] : dn) any = any || n.size()>0;
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
        const rvec_t& n=dn.at(s);
        if (n.size()>0) ch.n=n;                                // a source-less channel keeps its n
        if (ch.n.size()!=NumCoefficients()) ch.n=rvec_t(NumCoefficients(),0.0);
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
        static const bool trace = qchem::Diagnostics::Enabled("u_trace");
        if (trace) std::cout<<os.str()<<std::endl; else report::Log(os.str());
    }
}

void Hubbard_U::PrepareSlots(const cbs_t* bs, const cChargeDensity* cd) const
{
    cDynamic_HT_Imp::PrepareSlots(bs, cd);
    Dynamic_HT_RealBlock_Imp::PrepareSlots(bs, cd);
    if (!bs) return;   // a refresh with no block list (energy evaluation): the stamp only (KAxisPlan P2)
    // The mixed-aware walk (RealComplexPlan step 3c-2): a real TRIM block inside the complex-faced set
    // answers GetRealIBS; every other block is a Bloch block.  The pairs are geometry-fixed, so build them all now.
    for (size_t i=0;i<bs->GetNumIBS();i++)
        if (const auto* r=bs->GetRealIBS(i)) Projector<double>(dynamic_cast<const BasisSet::Orbital_DFT_IBS<double,dcmplx>&>(*r));
        else                                 Projector<dcmplx>(dynamic_cast<const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>&>(*(*bs)[i]));
    // The site groups on the manifolds (increment 2): the shell layout is uniform across blocks, so the first
    // block's shells give every manifold's representation.
    if (bs->GetNumIBS()>0 && itsEigenForm)
    {
        if (PseudoAtomic())
            for (size_t M=0;M<itsManifolds.size();M++) BuildSymmetry(M, SetShells(M));
        else
        {
            const Selection sel = bs->GetRealIBS(0) ? Select(*bs->GetRealIBS(0)) : Select(*(*bs)[0]);
            for (size_t M=0;M<itsManifolds.size();M++) BuildSymmetry(M, sel.shells[M]);
        }
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
// The bare integrals pair with ManifoldCoefficients: over the functions Q_M c is expressed on.  A SITE-SHELL manifold's are
// the block's columns carried through the on-site contraction; a PSEUDO-ATOMIC manifold's are the set's own functions (one
// centre, orthonormal within the manifold, central cell -- HubbardProjectorPlan D6), which the set answers itself.
template <class U> static std::vector<BasisSet::ERI4Block> IntegralsOf(const LowdinProjector<U>& P, const BasisSet::Orbital_DFT_IBS<U,dcmplx>& orb,
                                                                      const BasisSet::GaussianAtomCentredSet* set)
{
    if (set)
    {
        std::vector<BasisSet::ERI4Block> out;
        for (size_t M=0;M<P.NumManifolds();M++)
        {
            std::vector<size_t> cols(set->Size(M));
            for (size_t a=0;a<cols.size();a++) cols[a]=set->Offset(M)+a;
            out.push_back(set->BareCoulomb(cols));
        }
        return out;
    }
    const auto* src=dynamic_cast<const BasisSet::BareCoulombSource*>(&orb);
    if (!src) throw std::runtime_error("Hubbard_U::ManifoldIntegrals: the orbital block cannot deliver bare two-electron "
                                       "integrals over a function subset (no BareCoulombSource face)");
    std::vector<BasisSet::ERI4Block> out;
    for (size_t M=0;M<P.NumManifolds();M++)
    {
        out.push_back(src->BareCoulomb(P.Columns(M)));   // a LOWDIN manifold's functions ARE its columns
    }
    return out;
}
std::vector<BasisSet::ERI4Block> Hubbard_U::ManifoldIntegrals(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& orb) const
{return IntegralsOf<double>(Projector<double>(orb), orb, itsSet.get());}
std::vector<BasisSet::ERI4Block> Hubbard_U::ManifoldIntegrals(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& orb) const
{return IntegralsOf<dcmplx>(Projector<dcmplx>(orb), orb, itsSet.get());}
template <class U> static std::vector<mat_t<U>> CoefficientsOf(const LowdinProjector<U>& P, const mat_t<U>& C)
{
    std::vector<mat_t<U>> out;
    for (size_t M=0;M<P.NumManifolds();M++) out.push_back(mat_t<U>(P.Coefficients(M)*C));   // Q_M C
    return out;
}
std::vector<mat_t<double>> Hubbard_U::ManifoldCoefficients(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& orb, const mat_t<double>& C) const
{return CoefficientsOf<double>(Projector<double>(orb), C);}
std::vector<mat_t<dcmplx>> Hubbard_U::ManifoldCoefficients(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& orb, const mat_t<dcmplx>& C) const
{return CoefficientsOf<dcmplx>(Projector<dcmplx>(orb), C);}
void Hubbard_U::SetU(size_t M, double U)
{
    if (M>=itsManifolds.size()) throw std::out_of_range("Hubbard_U::SetU: no such manifold");
    itsManifolds[M].U=U;
    for (double& u : itsManifolds[M].Uirrep) u=U;          // shell-averaged: every slot takes the new U
    ApplyWrite();
}
void Hubbard_U::SetPerturbation(size_t M, double alpha)
{
    if (M>=itsManifolds.size()) throw std::out_of_range("Hubbard_U::SetPerturbation: no such manifold");
    itsManifolds[M].alpha=alpha;
    ApplyWrite();
}
// UNFREEZING re-arms the refresh: the version stamp is the FREEZE-time density's, and a later density must not be
// mistaken for it (nor the cached V_k blocks, built from the held n, for the density's own).
void Hubbard_U::FreezeOccupations(bool frozen)
{
    const bool was=itsFrozen;
    itsFrozen=frozen;
    if (was && !frozen)
    {
        itsOccVersion=size_t(-1);
        cDynamic_HT_Imp::InvalidateCache();
        Dynamic_HT_RealBlock_Imp::InvalidateRealCache();
    }
}
// ⚠ FROZEN WAS A TRAP BEFORE 2026-09-29: EnsureOccupations returns at once when frozen, so a SetU on a frozen
// term invalidated the version stamp and the caches -- and the next refresh rebuilt the SAME old W.  The write
// reached nothing.  Frozen means "hold n", never "hold U and alpha", so the held n is re-analysed here.
void Hubbard_U::ApplyWrite()
{
    itsOccVersion=size_t(-1);                               // W depends on U/alpha: the next refresh rebuilds it from n
    if (itsFrozen)
    {
        bool held=!itsChannels.empty();
        for (const auto& [s,ch] : itsChannels) held = held && ch.n.size()==NumCoefficients() && NumCoefficients()>0;
        if (held)
        {
            itsEU=0.0;
            for (auto& [s,ch] : itsChannels) itsEU+=Analyse(ch);
            if (itsGroup==SpinGroup::UnPolarized) itsEU*=2.0;
        }
    }
    cDynamic_HT_Imp::InvalidateCache();                     // the cached V_k blocks are the OLD W's (same density serial!)
    Dynamic_HT_RealBlock_Imp::InvalidateRealCache();
}
template <class U> static std::vector<mat_t<U>> LowdinOf(const LowdinProjector<U>& P, const mat_t<U>& C)
{
    std::vector<mat_t<U>> out;
    for (size_t M=0;M<P.NumManifolds();M++) out.push_back(mat_t<U>(blazem::ctrans(P.T(M))*C));   // T^dagger C
    return out;
}
std::vector<mat_t<double>> Hubbard_U::ProjectorAmplitudes(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& orb, const mat_t<double>& C) const
{return LowdinOf<double>(Projector<double>(orb), C);}
std::vector<mat_t<dcmplx>> Hubbard_U::ProjectorAmplitudes(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& orb, const mat_t<dcmplx>& C) const
{return LowdinOf<dcmplx>(Projector<dcmplx>(orb), C);}
std::vector<HubbardChannel> Hubbard_U::Channels() const
{
    std::vector<HubbardChannel> out;
    for (const auto& M : itsManifolds)
    {
        bool carries = (M.U!=0.0);
        for (double u : M.Uirrep) carries = carries || (u!=0.0);
        out.push_back({M.site, M.l, carries});
    }
    return out;
}

HubbardSlots Hubbard_U::SlotProjectors(size_t M) const
{
    if (M>=itsManifolds.size()) throw std::out_of_range("Hubbard_U::SlotProjectors: no such manifold");
    HubbardSlots out;
    if (!itsEigenForm) return out;                                   // the population form: no symmetry, one subspace
    if (M>=itsSym.size() || !itsSym[M])
        throw std::logic_error("Hubbard_U::SlotProjectors: the term has not seen a basis yet (slots are built in PrepareSlots)");
    const ManifoldSymmetry& S=*itsSym[M];
    out.components=itsComponents[M];
    for (const auto& s : S.Slots())
    {
        std::ostringstream lab; lab<<"Γ"<<s.irrep; if (S.Graded()) lab<<"<"<<s.parent;
        out.slots.push_back({lab.str(), S.Irreps()[s.irrep].dim, rmat_t(S.SiteProjector(s.irrep)*S.GreyProjector(s.parent))});
    }
    return out;
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
    EnsureOccupations(cd);   // a LOOKUP: the refresh (the contract, KAxisPlan P2) filled it for this density
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

// FROZEN, the energy is the +U functional LINEARISED about the held occupations n0:
//     E = E_U(n0) + Tr[W0 (n - n0)]      (W0 includes the alpha*1 shift)
// -- whose derivative in D is exactly the frozen Fock term W0, so the SCF's energy and its Fock matrix are the
// SAME functional again.  Reporting E_U(n0) alone (a constant) left the energy blind to the alpha*P and V_Hub
// terms the Fock carries: the energy-guarded re-damp then fought the Fock and a frozen +-alpha SCF on Si
// oscillated for 80 iterations (2026-09-29).  Unfrozen, n0 = n and the correction is identically zero.
void Hubbard_U::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    EnsureOccupations(cd);                                        // the energy pass's density (rho_out)
    double E=itsEU;
    if (itsFrozen && cd)
    {
        const std::map<Spin,rvec_t> dn=DensityOccupations(cd);
        double lin=0.0;
        for (const auto& [s,ch] : itsChannels)
        {
            const rvec_t& n=dn.at(s);
            if (n.size()!=ch.W.size() || ch.n.size()!=ch.W.size()) continue;   // no D: nothing moved
            for (size_t i=0;i<n.size();i++) lin+=ch.W[i]*(n[i]-ch.n[i]);
        }
        if (itsGroup==SpinGroup::UnPolarized) lin*=2.0;           // both (identical) channels, as itsEU
        E+=lin;
    }
    te.Add("E_U", E, EnergyRole::Potential);
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
        if (man.alpha!=0.0) os<<" alpha="<<man.alpha<<" Ha (static shift)";
        if (M<itsSym.size() && itsSym[M]) { os<<" slots:"; itsSym[M]->WriteSlots(os, man.U, man.Uirrep); }
        os<<";";
    }
    os<<(itsEigenForm ? "  [eigenvalue form]" : "  [DIAGONAL populations, CP2K's form]");
    return os<<std::endl;
}

// THE RESPONSE (R2 / H4): zero when frozen or U = 0 -- the only two cases a linear response asks for today.
void Hubbard_U::RefreshForDensity(const cbs_t*, const cChargeDensity*, const TransitionDensity<dcmplx>&) const
{
    if (itsFrozen) return;
    for (const auto& M : itsManifolds)
    {
        bool zero = (M.U==0.0);
        for (double u : M.Uirrep) zero = zero && (u==0.0);
        if (!zero)
            throw std::logic_error("Hubbard_U: the UNFROZEN +U response kernel is not built -- freeze the occupations "
                                   "(Timrov eq 20: V_Hub held at its ground-state value) or run at U = 0 (the U_0 case)");
    }
}

mat_t<dcmplx> Hubbard_U::GetMatrix(const cobs_t* bra, const cobs_t* ket, const Spin&, const TransitionDensity<dcmplx>& delta) const
{
    RefreshForDensity(nullptr, nullptr, delta);   // the same gate, so a direct caller cannot skip it
    return mat_t<dcmplx>(bra->GetNumFunctions(), ket->GetNumFunctions(), dcmplx(0.0));
}

} // namespace
