// File: BasisSet/Molecule/Evaluators/PG_LibCint/Imp/Evaluator.C
//
// libcint packing + matrix assembly for the Cartesian PG basis.  See Evaluator.C for the why; this file
// is the how:
//   1. Regroup PGData's flattened (radial x polarization) components back into SHELLS (a run of components
//      sharing one radial and one total-L).  Each shell becomes one libcint contracted Cartesian basis
//      function (NCTR=1) with the radial's exponents + (normalization-folded) contraction coeffs.
//   2. Per shell, build the libcint-cart -> PG-order component permutation (lc2pg) by matching (lx,ly,lz).
//   3. Renormalize every component by 1/sqrt(native self-overlap) so the delivered matrices use PG's
//      unit-self-overlap convention (== PGData::ns, but computed from libcint to stay convention-free).
// Matrices are then assembled by looping libcint shell tuples (the library's natural grain), scattering
// each dense block into the PG-ordered result with the permutation + scales.
module;
#include <vector>
#include <cmath>
#include <cassert>
#include <ostream>
#include <string>
extern "C" {
#include "cint_funcs.h"   // pulls cint.h: slot constants, FINT/CACHE_SIZE_T, the int*_cart functions
// cint_funcs.h omits prototypes for the 3-centre integrals (the symbols are in libcint); declare them.
CINTIntegralFunction int3c1e_cart;   // 3-centre overlap  <ab|c>
CINTIntegralFunction int3c2e_cart;   // 3-centre Coulomb  (ab|c)
}
module qchem.BasisSet.Molecule.Evaluators.PG_LibCint;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.PGData;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Polarization;
import qchem.Structure;
import qchem.Types;

namespace BasisSet::Molecule::Evaluators::PG_LibCint
{
using PG_Cart_MnD::PGData;
using PG_Cart_MnD::Polarization;
using PG_Cart_MnD::GaussianRF;
using IntFn = CINTIntegralFunction*;   // (out,dims,shls,atm,natm,bas,nbas,env,opt,cache)

static int ncart(int l) {return (l+1)*(l+2)/2;}

// --- internal layout ------------------------------------------------------------------------------
struct AtomRec { int Z; double x,y,z; };
struct ShellDesc
{
    int                 atom;   // index into the (shared) atom table
    int                 l;
    std::vector<double> exp;    // primitive exponents
    std::vector<double> coef;   // PG normalization-folded contraction coefficients
    std::vector<size_t> lc2pg;  // libcint component -> global output index (Cartesian: PG order; spherical:
                                // sequential = libcint's native real-spherical order)
};

struct NR_Evaluator::Imp
{
    bool                   spherical=false;
    std::vector<AtomRec>   atoms;
    std::vector<ShellDesc> shells;
    rvec_t                 scale;   // per global component: 1/sqrt(native self-overlap)
    size_t                 N=0;     // number of delivered components ((l+1)(l+2)/2 cart, 2l+1 sph, per shell)

    int dim(int l) const {return spherical ? 2*l+1 : (l+1)*(l+2)/2;}   // components per shell

    // Build atm/bas/env for the shared atom table + a concatenation of shell lists (slots).  Shell ids in
    // the returned bas run in the order the ShellDesc pointers are given.
    static void BuildMol(const std::vector<AtomRec>& atoms, const std::vector<const ShellDesc*>& sh,
                         std::vector<int>& atm, std::vector<int>& bas, std::vector<double>& env)
    {
        atm.assign(atoms.size()*ATM_SLOTS, 0);
        env.assign(PTR_ENV_START, 0.0);
        for (size_t i=0;i<atoms.size();++i)
        {
            atm[CHARGE_OF+ATM_SLOTS*i]=atoms[i].Z;
            atm[PTR_COORD+ATM_SLOTS*i]=(int)env.size();
            env.push_back(atoms[i].x); env.push_back(atoms[i].y); env.push_back(atoms[i].z);
        }
        bas.assign(sh.size()*BAS_SLOTS, 0);
        for (size_t s=0;s<sh.size();++s)
        {
            bas[ATOM_OF +BAS_SLOTS*s]=sh[s]->atom;
            bas[ANG_OF  +BAS_SLOTS*s]=sh[s]->l;
            bas[NPRIM_OF+BAS_SLOTS*s]=(int)sh[s]->exp.size();
            bas[NCTR_OF +BAS_SLOTS*s]=1;
            bas[PTR_EXP +BAS_SLOTS*s]=(int)env.size();
            for (double e:sh[s]->exp ) env.push_back(e);
            bas[PTR_COEFF+BAS_SLOTS*s]=(int)env.size();
            for (double c:sh[s]->coef) env.push_back(c);
        }
    }
};

// libcint Cartesian component order for momentum l (CINTcart_comp: lx desc, ly desc), as Polarizations.
static std::vector<Polarization> LibcintCartOrder(int l)
{
    std::vector<Polarization> v;
    for (int lx=l; lx>=0; --lx)
        for (int ly=l-lx; ly>=0; --ly)
            v.push_back(Polarization(lx, ly, l-lx-ly));   // (n,l,m) = (x,y,z) powers
    return v;
}

NR_Evaluator::NR_Evaluator() = default;

NR_Evaluator::NR_Evaluator(const PGData& data, const Structure* cl)
{
    static_cast<PGData&>(*this) = data;             // own the component layout + order (size/Norm/ns)
    Init(cl);
}

// Build the libcint shell + atom tables from this evaluator's (already-populated) PGData and the structure.
void NR_Evaluator::Init(const Structure* cl, bool spherical)
{
    itsImp.reset(new Imp);
    Imp& m = *itsImp;
    m.spherical = spherical;
    size_t Ncart = PGData::size();                  // # Cartesian components flattened in the PGData

    // Atom table (shell centres + nuclear charges) from the structure.
    for (auto a:*cl) m.atoms.push_back(AtomRec{a->itsZ, a->itsR.x, a->itsR.y, a->itsR.z});
    auto findAtom = [&](const rvec3_t& R)->int
    {
        for (size_t i=0;i<m.atoms.size();++i)
            if (std::fabs(m.atoms[i].x-R.x)+std::fabs(m.atoms[i].y-R.y)+std::fabs(m.atoms[i].z-R.z) < 1e-8)
                return (int)i;
        assert(false && "shell centre does not match a nucleus"); return 0;
    };

    // Regroup the flattened Cartesian components into shells: a new shell begins when the radial pointer OR
    // the total L changes (PGData lays out block-by-block, MakePolarizations groups by L, so runs are
    // contiguous).  Packing into libcint is the SAME (Cartesian shell defs) for both modes; only the
    // delivered component count + order differ (Cartesian -> PG permutation; spherical -> 2l+1 sequential).
    size_t gout=0;
    for (size_t i=0;i<Ncart;)
    {
        const GaussianRF* r = radials[i];
        int L = pols[i].GetTotalL();
        size_t j=i;
        std::vector<size_t> gidx;                   // global Cartesian component indices (PG order)
        while (j<Ncart && radials[j]==r && pols[j].GetTotalL()==L) gidx.push_back(j++);
        assert((int)gidx.size()==ncart(L) && "shell is not a complete Cartesian set");

        ShellDesc sh;
        sh.atom = findAtom(r->GetCenter());
        sh.l    = L;
        { rvec_t es=r->GetExponents(); const rvec_t& cs=r->GetCoeffs();
          for (size_t k=0;k<es.size();++k){ sh.exp.push_back(es[k]); sh.coef.push_back(cs[k]); } }
        if (spherical)
            for (int p=0;p<m.dim(L);++p) sh.lc2pg.push_back(gout++);    // libcint-native sph order, sequential
        else
            for (const Polarization& lp:LibcintCartOrder(L))           // libcint-cart order -> PG global index
            {
                size_t hit=gidx.size();
                for (size_t k=0;k<gidx.size();++k) if (pols[gidx[k]]==lp){ hit=k; break; }
                assert(hit<gidx.size() && "missing Cartesian component"); sh.lc2pg.push_back(gidx[hit]);
            }
        m.shells.push_back(std::move(sh));
        i=j;
    }
    m.N = spherical ? gout : Ncart;

    // Per-component scale = 1/sqrt(native self-overlap) (libcint does not unit-normalize Cartesians; the
    // _sph harmonics are already normalized, so this is ~1 there, but compute it anyway -- convention-free).
    m.scale.resize(m.N);
    std::vector<const ShellDesc*> self; for (auto& s:m.shells) self.push_back(&s);
    std::vector<int> atm,bas; std::vector<double> env; Imp::BuildMol(m.atoms, self, atm,bas,env);
    int natm=(int)m.atoms.size(), nbas=(int)m.shells.size();
    IntFn ovlp = spherical ? int1e_ovlp_sph : int1e_ovlp_cart;
    std::vector<double> buf;
    for (size_t s=0;s<m.shells.size();++s)
    {
        int d=m.dim(m.shells[s].l); buf.assign((size_t)d*d, 0.0);
        int shls[2]={(int)s,(int)s};
        ovlp(buf.data(), nullptr, shls, atm.data(), natm, bas.data(), nbas, env.data(), nullptr, nullptr);
        for (int p=0;p<d;++p) m.scale[m.shells[s].lc2pg[p]] = 1.0/std::sqrt(buf[p+d*p]);
    }
}

NR_Evaluator::~NR_Evaluator() = default;

size_t NR_Evaluator::size() const {return itsImp ? itsImp->N     : PGData::size();}
rvec_t NR_Evaluator::Norm() const {return itsImp ? itsImp->scale : ns;}

std::ostream& NR_Evaluator::Write(std::ostream& os) const
{
    return os << "PG_LibCint::NR_Evaluator[" << size() << "]";
}

// --- 1E: loop shell pairs, scatter each block into the symmetric PG-ordered matrix ----------------
rsmat_t NR_Evaluator::Build1(int which) const
{
    const Imp& m=*itsImp;
    rsmat_t M(m.N);
    std::vector<const ShellDesc*> self; for (auto& s:m.shells) self.push_back(&s);
    std::vector<int> atm,bas; std::vector<double> env; Imp::BuildMol(m.atoms, self, atm,bas,env);
    int natm=(int)m.atoms.size(), nbas=(int)m.shells.size();
    IntFn fn   = m.spherical
        ? (which==0 ? int1e_ovlp_sph  : which==1 ? int1e_kin_sph  : int1e_nuc_sph)
        : (which==0 ? int1e_ovlp_cart : which==1 ? int1e_kin_cart : int1e_nuc_cart);
    double pre = which==1 ? 2.0 : 1.0;   // KineticMatrix is <-nabla^2> = 2 x libcint's T = -1/2 nabla^2
    std::vector<double> buf;
    for (size_t si=0;si<m.shells.size();++si)
        for (size_t sj=si;sj<m.shells.size();++sj)
        {
            int di=m.dim(m.shells[si].l), dj=m.dim(m.shells[sj].l);
            buf.assign((size_t)di*dj, 0.0);
            int shls[2]={(int)si,(int)sj};
            fn(buf.data(), nullptr, shls, atm.data(), natm, bas.data(), nbas, env.data(), nullptr, nullptr);
            for (int p=0;p<di;++p){ size_t gi=m.shells[si].lc2pg[p];
                for (int q=0;q<dj;++q){ size_t gj=m.shells[sj].lc2pg[q];
                    M(gi,gj) = pre * buf[p+di*q] * m.scale[gi]*m.scale[gj]; } }
        }
    return M;
}
rsmat_t NR_Evaluator::OverlapMatrix()             const {return Build1(0);}
rsmat_t NR_Evaluator::KineticMatrix()             const {return Build1(1);}
rsmat_t NR_Evaluator::NuclearMatrix(const Structure*) const {return Build1(2);} // atm already carries the nuclei

// --- 3-centre <ab|c>: one symmetric (ia,ib) block per fit component ic ----------------------------
ERI3<double> NR_Evaluator::Build3(const NR_Evaluator& fit, int which) const
{
    const Imp& a=*itsImp; const Imp& c=*fit.itsImp;
    assert(!a.spherical && "3-centre (DFT) is Cartesian-only for the libcint evaluator");
    std::vector<const ShellDesc*> sh; for (auto& s:a.shells) sh.push_back(&s);
    size_t cOff=sh.size();           for (auto& s:c.shells) sh.push_back(&s);
    std::vector<int> atm,bas; std::vector<double> env; Imp::BuildMol(a.atoms, sh, atm,bas,env);
    int natm=(int)a.atoms.size(), nbas=(int)sh.size();
    IntFn fn = which==0 ? int3c1e_cart : int3c2e_cart;   // overlap : Coulomb

    ERI3<double> s3; s3.reserve(c.N); for (size_t i=0;i<c.N;++i) s3.push_back(rsmat_t(a.N));
    std::vector<double> buf;
    for (size_t sk=0;sk<c.shells.size();++sk)
        for (size_t si=0;si<a.shells.size();++si)
            for (size_t sj=si;sj<a.shells.size();++sj)
            {
                int di=a.dim(a.shells[si].l), dj=a.dim(a.shells[sj].l), dk=c.dim(c.shells[sk].l);
                buf.assign((size_t)di*dj*dk, 0.0);
                int shls[3]={(int)si,(int)sj,(int)(cOff+sk)};
                fn(buf.data(), nullptr, shls, atm.data(), natm, bas.data(), nbas, env.data(), nullptr, nullptr);
                for (int p=0;p<di;++p){ size_t gia=a.shells[si].lc2pg[p];
                    for (int q=0;q<dj;++q){ size_t gib=a.shells[sj].lc2pg[q];
                        for (int rr=0;rr<dk;++rr){ size_t gic=c.shells[sk].lc2pg[rr];
                            s3[gic](gia,gib) = buf[p+di*(q+dj*rr)] * a.scale[gia]*a.scale[gib]*c.scale[gic]; } } }
            }
    return s3;
}
ERI3<double> NR_Evaluator::OverlapThreeC_Matrix  (const NR_Evaluator& fit) const {return Build3(fit,0);}
ERI3<double> NR_Evaluator::RepulsionThreeC_Matrix(const NR_Evaluator& fit) const {return Build3(fit,1);}

// --- 4-centre (ab|cd): dense reordered+renormalized tensor over four (this/partner) slots ----------
std::vector<double> NR_Evaluator::Compute4(const NR_Evaluator& B, const NR_Evaluator& C,
                                           const NR_Evaluator& D) const
{
    const Imp& A=*itsImp; const Imp& mb=*B.itsImp; const Imp& mc=*C.itsImp; const Imp& md=*D.itsImp;
    size_t Na=A.N,Nb=mb.N,Nc=mc.N,Nd=md.N;
    std::vector<const ShellDesc*> sh;
    size_t oA=0;            for (auto& s:A.shells ) sh.push_back(&s);
    size_t oB=sh.size();    for (auto& s:mb.shells) sh.push_back(&s);
    size_t oC=sh.size();    for (auto& s:mc.shells) sh.push_back(&s);
    size_t oD=sh.size();    for (auto& s:md.shells) sh.push_back(&s);
    std::vector<int> atm,bas; std::vector<double> env; Imp::BuildMol(A.atoms, sh, atm,bas,env);
    int natm=(int)A.atoms.size(), nbas=(int)sh.size();

    std::vector<double> G(Na*Nb*Nc*Nd, 0.0);
    auto idx=[&](size_t ia,size_t ib,size_t ic,size_t id){return ((ia*Nb+ib)*Nc+ic)*Nd+id;};
    std::vector<double> buf;
    for (size_t sa=0;sa<A.shells.size();++sa)
     for (size_t sb=0;sb<mb.shells.size();++sb)
      for (size_t sc=0;sc<mc.shells.size();++sc)
       for (size_t sd=0;sd<md.shells.size();++sd)
       {
           int da=A.dim(A.shells[sa].l), db=mb.dim(mb.shells[sb].l),
               dc=mc.dim(mc.shells[sc].l), dd=md.dim(md.shells[sd].l);
           buf.assign((size_t)da*db*dc*dd, 0.0);
           int shls[4]={(int)(oA+sa),(int)(oB+sb),(int)(oC+sc),(int)(oD+sd)};
           IntFn fn = A.spherical ? int2e_sph : int2e_cart;
           fn(buf.data(), nullptr, shls, atm.data(), natm, bas.data(), nbas, env.data(), nullptr, nullptr);
           for (int p=0;p<da;++p){ size_t gia=A.shells[sa].lc2pg[p];
            for (int q=0;q<db;++q){ size_t gib=mb.shells[sb].lc2pg[q];
             for (int rr=0;rr<dc;++rr){ size_t gic=mc.shells[sc].lc2pg[rr];
              for (int s4=0;s4<dd;++s4){ size_t gid=md.shells[sd].lc2pg[s4];
                G[idx(gia,gib,gic,gid)] = buf[p+da*(q+db*(rr+dc*s4))]
                    * A.scale[gia]*mb.scale[gib]*mc.scale[gic]*md.scale[gid]; } } } }
       }
    return G;
}

ERI4 NR_Evaluator::DirectMatrix(const NR_Evaluator& p) const
{
    size_t Na=itsImp->N, Nc=p.itsImp->N;
    std::vector<double> G=Compute4(*this, p, p);   // (a a | c c): A=this,B=this,C=p,D=p
    auto g=[&](size_t ia,size_t ib,size_t ic,size_t id){return G[((ia*Na+ib)*Nc+ic)*Nc+id];};
    ERI4 J(Na,Nc);
    for (size_t ia=0;ia<Na;++ia) for (size_t ib=ia;ib<Na;++ib)
    {
        rsmat_t& Jab=J(ia,ib);
        for (size_t ic=0;ic<Nc;++ic) for (size_t id=ic;id<Nc;++id) Jab(ic,id)=g(ia,ib,ic,id);
    }
    return J;
}

ERI4 NR_Evaluator::ExchangeMatrix(const NR_Evaluator& pb) const
{
    size_t Na=itsImp->N, Nb=pb.itsImp->N;
    std::vector<double> G=Compute4(pb, *this, pb);   // slots (a b | a b): A=this,B=pb,C=this,D=pb
    auto g=[&](size_t ia,size_t ib,size_t ic,size_t id){return G[((ia*Nb+ib)*Na+ic)*Nb+id];};
    ERI4 K(Na,Nb);
    for (size_t ia=0;ia<Na;++ia) for (size_t ib=0;ib<Nb;++ib) for (size_t ic=ia;ic<Na;++ic)
    {
        rsmat_t& Kac=K(ia,ic);
        for (size_t id=0;id<Nb;++id)
        {
            double v=g(ia,ib,ic,id);
            if (ib==id)     Kac(ib,id) = v;
            else if (ib<id) Kac(ib,id)+= 0.5*v;
            else            Kac(id,ib)+= 0.5*v;
        }
    }
    return K;
}

} //namespace
