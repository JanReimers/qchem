// File: GaussianRF.C  Implementation of the single radial type (contracted Gaussian) + the internal
// primitive helper PrimGaussian that carries the M&D integral kernels.
module;
#include <iostream>
#include <iomanip>
#include <cassert>
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <utility>

module qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Internal.Hermite;   // Hermite1/2/3
import qchem.BasisSet.Molecule.Evaluators.Internal.MnD;   // generic MnD core (RNLM, ...)
import qchem.BasisSet.Internal.DB_Cache;   // theCache<double>(), Register/GetCache{2,3}
import qchem.BasisSet.Internal.Cache2;     // Cache2, Cacheable2, Cache2_Client
import qchem.BasisSet.Internal.Cache3;     // Cache3, Cacheable3, Cache3_Client

import qchem.Blaze;     // rvec_t (itsCoeff)
import qchem.Structure;
import qchem.Math;

namespace BasisSet::Molecule::Evaluators::PG_Cart_MnD
{
using Evaluators::Internal::MnD::RNLM;
using Evaluators::Internal::MnD::Index3;   // plain (n,l,m) Hermite index for the RNLM lookup

//#######################################################################
//
//   Content-based identity registry -- replaces the old per-object UniqueID.  Each distinct primitive
//   (exponent, centre, max-L) and each distinct charge-distribution pair gets a stable index from ONE
//   shared counter, so primitive ids and Ω ids never collide in the RNLM cache (which keys on a mix of
//   the two).  Content keys mean two equal primitives in different IBSs share a Cache2/3 entry (the
//   reuse win, as the atom ExponentGrouper does for Cache4).  The registry is a serialisable symbol
//   table: a future per-job disk cache would persist it next to Cache{2,3} -- nothing here blocks that.
//
namespace
{
    struct PrimKey
    {
        double e, x, y, z; int L;
        bool operator<(const PrimKey& o) const
        {
            if (e!=o.e) return e<o.e;
            if (x!=o.x) return x<o.x;
            if (y!=o.y) return y<o.y;
            if (z!=o.z) return z<o.z;
            return L<o.L;
        }
    };
    class IdRegistry
    {
    public:
        size_t Primitive  (double e, const rvec3_t& R, int L) {return intern(itsPrims, PrimKey{e,R.x,R.y,R.z,L});}
        size_t ChargeDist (size_t a, size_t b)                {return intern(itsPairs, std::make_pair(a,b));}
    private:
        template <class Map, class Key> size_t intern(Map& m, const Key& k)
        {
            auto [it, fresh] = m.try_emplace(k, itsNext);
            if (fresh) ++itsNext;
            return it->second;
        }
        size_t                                    itsNext = 0;
        std::map<PrimKey,size_t>                   itsPrims;
        std::map<std::pair<size_t,size_t>,size_t>  itsPairs;
    };
    IdRegistry& theIds() {static IdRegistry r; return r;}
}

//#######################################################################
//
//   Charge distribution Ω + the global Cache2/Cache3 access points.
//   All PG charge-distribution caching lives in process-global Cache2s, keyed by registry indices so entries
//   are shared wherever the underlying objects are shared (e.g. SALC irreps over one raw basis):
//     - Ω keyed by the primitive pair (primA.ID, primB.ID)
//     - 3/4-centre RNLM keyed by the two charge-distribution ids (ab.ID, c.ID)
//   (The 2-centre self-RNLM is a lazy member of each Ω, not cached here.)
//
namespace
{
    struct OmegaClient : Cache2_Client
    {
        std::string RadialType() const override {return "PG.Omega";}
        Cache2*     MakeCache2() const override {return new Cache2;}
    };
    // Memoize the Cache2* per theGlobalCache instance.  GetCache2 is a std::map<std::string,...>::find
    // (string compare), and findΩ runs it on EVERY call from the innermost contraction loops -- profiling
    // (callgrind, water/dzvp) put that per-call string lookup at ~13% of MnD integral time.  Caching the
    // pointer turns it into a pointer compare; the owner guard re-fetches if theGlobalCache is swapped
    // (e.g. a test teardown), so it stays correct.
    const Cache2* OmegaCache()
    {
        static decltype(&BasisSet::theCache<double>()) owner = nullptr;  // typed cache-instance guard (no void*)
        static const Cache2* cache = nullptr;
        if (owner != &BasisSet::theCache<double>())
        {
            static OmegaClient c;
            BasisSet::theCache<double>().Register(&c);
            cache = BasisSet::theCache<double>().GetCache2("PG.Omega");
            owner = &BasisSet::theCache<double>();
        }
        return cache;
    }

    // Cacheable2 wrapper so an RNLM can live in a Cache2 (RNLM itself is a low-level MnD type and
    // must not depend on the cache framework).
    struct RNLM_C2 : public Cacheable2
    {
        RNLM rnlm;
        RNLM_C2(int L, double alpha, const rvec3_t& R) : rnlm(L, alpha, R) {}
        bool   isSupported(const Cache2_Client*) const override {return false;}
        size_t RAMsize() const override {return sizeof(RNLM_C2);}
    };
    struct RNLMClient : Cache2_Client
    {
        std::string RadialType() const override {return "PG.RNLM";}
        Cache2*     MakeCache2() const override {return new Cache2;}
    };
    const Cache2* RNLMCache()   // memoized per theGlobalCache instance (see OmegaCache)
    {
        static decltype(&BasisSet::theCache<double>()) owner = nullptr;  // typed cache-instance guard (no void*)
        static const Cache2* cache = nullptr;
        if (owner != &BasisSet::theCache<double>())
        {
            static RNLMClient c;
            BasisSet::theCache<double>().Register(&c);
            cache = BasisSet::theCache<double>().GetCache2("PG.RNLM");
            owner = &BasisSet::theCache<double>();
        }
        return cache;
    }

    // Cacheable3 wrapper holding the 3-centre Hermite block.
    struct H3_C3 : public Cacheable3
    {
        std::unique_ptr<Hermite3> h3;
        explicit H3_C3(Hermite3* p) : h3(p) {}
        bool   isSupported(const Cache3_Client*) const override {return false;}
        size_t RAMsize() const override {return sizeof(H3_C3);}
    };
    struct H3Client : Cache3_Client
    {
        std::string RadialType() const override {return "PG.H3";}
        Cache3*     MakeCache3() const override {return new Cache3;}
    };
    const Cache3* H3Cache()   // memoized per theGlobalCache instance (see OmegaCache)
    {
        static decltype(&BasisSet::theCache<double>()) owner = nullptr;  // typed cache-instance guard (no void*)
        static const Cache3* cache = nullptr;
        if (owner != &BasisSet::theCache<double>())
        {
            static H3Client c;
            BasisSet::theCache<double>().Register(&c);
            cache = BasisSet::theCache<double>().GetCache3("PG.H3");
            owner = &BasisSet::theCache<double>();
        }
        return cache;
    }
}

const Ω& findΩ(const GData& a,const GData& b)
{
    const Cacheable2& cd = OmegaCache()->get(a.ID, b.ID,
        [&a,&b]() -> const Cacheable2* { return new Ω(a,b); });
    return static_cast<const Ω&>(cd);
}

const RNLM& findRNLM(const GData& ab,const GData& c)
{
    const Cacheable2& w = RNLMCache()->get(ab.ID, c.ID,
        [&ab,&c]() -> const Cacheable2*
        {
            double alpha = ab.α*c.α/(ab.α+c.α);
            return new RNLM_C2(ab.L+c.L, alpha, ab.R-c.R);
        });
    return static_cast<const RNLM_C2&>(w).rnlm;
}

// 3-centre Hermite block, cached by primitive triple.  Now that Ω lives in this module there is no
// Ω<->radial module cycle, so the block is built inline here (gc is centre C) -- no make-lambda.
const Hermite3& findH3(const PrimGaussian* ga, const PrimGaussian* gb, const PrimGaussian* gc)
{
    const Cacheable3& w = H3Cache()->get(ga->GetID(), gb->GetID(), gc->GetID(),
        [&]() -> const Cacheable3* { return new H3_C3(gc->GetH3(*ga,*gb)); });
    return *static_cast<const H3_C3&>(w).h3;
}

//--------------------------------------------------------------------------------------------------
//  Ω (charge distribution for a primitive pair) methods.
//--------------------------------------------------------------------------------------------------
Ω::~Ω() { delete itsSelfRNLM; }

// 2-centre self-auxiliary: RNLM(Ltotal, ab/αₚ, AB) -- built once, owned by this Ω.
const RNLM& Ω::SelfRNLM() const
{
    if (!itsSelfRNLM) itsSelfRNLM = new RNLM(Ltotal, ab/αₚ, AB);
    return *itsSelfRNLM;
}

std::vector<std::vector<Polarization>> Ω::theNMLs;

static std::vector<Polarization> MakeAllPolarizations(int Lmax)
{
    std::vector<Polarization> list;
    for (int n=0; n<=Lmax; n++)
        for (int l=0; l<=Lmax-n; l++)
            for (int m=0; m<=Lmax-n-l; m++)
                list.push_back(Polarization(n,l,m));
    return list;
}

void Ω::MakeNMLs()
{
    for (int L=0; L<=10; L++) theNMLs.push_back(MakeAllPolarizations(L));
}

Ω::Ω(const GData& g1,const GData& g2)
    : itsIndex(theIds().ChargeDist(g1.ID, g2.ID))   // content-based id (disjoint from primitive ids)
    , Ltotal(g1.L + g2.L)
    , a     (g1.α)
    , b     (g2.α)
    , ab    (a * b)
    , αₚ(a + b)
    , AB    (g1.R - g2.R)
    , P     ( (a*g1.R + b*g2.R) / αₚ)
    , Eij   ( exp(-ab / αₚ * (AB*AB)) )
    , H2    (αₚ, P - g1.R, P - g2.R, g1.L+1, g2.L+1)
{
    if (theNMLs.size()==0) MakeNMLs();
}

size_t Ω::RAMsize() const {return sizeof(Ω);} // includes the by-value Hermite2 block

//#######################################################################
//
//   PrimGaussian: primitive Gaussian + M&D kernels
//

PrimGaussian::PrimGaussian(double theExponent, const rvec3_t& theCenter, int theL)
    : itsIndex   (theIds().Primitive(theExponent, theCenter, theL))
    , itsExponent(theExponent)
    , itsCenter  (theCenter)
    , itsL       (theL)
    , itsH1      (0)
{
    if (itsExponent < 0)
    {
        std::cerr << "PrimGaussian exponent < 0" << std::endl;
        exit(-1);
    }
}

PrimGaussian::~PrimGaussian()
{
    delete itsH1;
}

const Hermite1& PrimGaussian::GetH1() const
{
    if (!itsH1) itsH1 = new Hermite1(itsExponent, itsL);
    return *itsH1;
}

//
//  Symmetric "gradient-dot-gradient" form of the kinetic-energy building block for a pair of
//  Cartesian Gaussians \f$\phi_a,\phi_b\f$.  The (non-relativistic) kinetic-energy operator is
//  \f$\hat T = -\tfrac12\nabla^2\f$, so the matrix element ultimately wanted is
//  \f[ T_{ab} = -\tfrac12\,\langle\phi_a|\nabla^2|\phi_b\rangle . \f]
//  Integrating by parts (the surface term vanishes for localized Gaussians) gives the symmetric,
//  manifestly positive-definite form returned here:
//  \f[ \texttt{Grad2}_{ab} \;=\; \langle\phi_a|-\nabla^2|\phi_b\rangle
//        \;=\; \int \nabla\phi_a\cdot\nabla\phi_b \, d^3r . \f]
//  So \b Grad2 is the matrix of \f$-\nabla^2\f$ (NOT \f$+\nabla^2\f$) -- the name reads as
//  \f$\nabla\!\cdot\!\nabla\f$, which is why it is positive.  The factor \f$\tfrac12\f$ is applied
//  later, at the boundary where the word "Kinetic" is introduced: \f$T = \tfrac12\,\texttt{Grad2}\f$.
//
//  Per Cartesian axis (x shown; y,z analogous) with \f$l_a=\f$ p1.n, \f$l_b=\f$ p2.n,
//  \f$\alpha=\f$ ab.a, \f$\beta=\f$ ab.b, and \f$S(\cdot,\cdot)=\f$ ab.H2(0,\f$\cdot,\cdot\f$) the
//  overlap building block:
//  \f[ t_{xx} = l_a l_b\,S(l_a\!-\!1,l_b\!-\!1) - 2\beta l_a\,S(l_a\!-\!1,l_b\!+\!1)
//             - 2\alpha l_b\,S(l_a\!+\!1,l_b\!-\!1) + 4\alpha\beta\,S(l_a\!+\!1,l_b\!+\!1) , \f]
//  and \f$\texttt{Grad2}=t_{xx}+t_{yy}+t_{zz}\f$.
//
static double GetGrad2(const Polarization& p1, const Polarization& p2, const Ω& ab)
{
    static Polarization p0(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
    const Polarization& P1 = p1;
    const Polarization& P2 = p2;

    double txx =   P1.n * P2.n * ab.H2(p0,P1-x,P2-x)
                 - 2*P1.n * ab.b * ab.H2(p0,P1-x,P2+x)
                 - 2*P2.n * ab.a * ab.H2(p0,P1+x,P2-x)
                 + 4*ab.ab        * ab.H2(p0,P1+x,P2+x);

    double tyy =   P1.l * P2.l * ab.H2(p0,P1-y,P2-y)
                 - 2*P1.l * ab.b * ab.H2(p0,P1-y,P2+y)
                 - 2*P2.l * ab.a * ab.H2(p0,P1+y,P2-y)
                 + 4*ab.ab        * ab.H2(p0,P1+y,P2+y);

    double tzz =   P1.m * P2.m * ab.H2(p0,P1-z,P2-z)
                 - 2*P1.m * ab.b * ab.H2(p0,P1-z,P2+z)
                 - 2*P2.m * ab.a * ab.H2(p0,P1+z,P2-z)
                 + 4*ab.ab        * ab.H2(p0,P1+z,P2+z);

    return txx+tyy+tzz;
}

double PrimGaussian::Overlap2C(const PrimGaussian* a, const PrimGaussian* b,
                               const Polarization& pa, const Polarization& pb)
{
    const Ω& ab = findΩ(a->GetGData(), b->GetGData());
    Polarization zero(0,0,0);
    return pow(Pi/ab.αₚ,1.5)*ab.Eij*ab.H2(zero,pa,pb);
}

double PrimGaussian::Repulsion2C(const PrimGaussian* a, const PrimGaussian* b,
                                 const Polarization& pa, const Polarization& pb)
{
    const Ω& ab = findΩ(a->GetGData(), b->GetGData());
    auto NLMs = Ω::GetNMLs(a->GetL());
    const Hermite1& H1a = a->GetH1();
    const Hermite1& H1b = b->GetH1();
    const RNLM& R = ab.SelfRNLM();

    double factor = 1.0/(ab.ab*sqrt(ab.αₚ));
    factor = (pb.GetTotalL()%2) ? -factor : factor;

    double s = 0.0;
    for (auto bNLM:NLMs)
    {
        if (bNLM > pa) continue;
        double ha = H1a(bNLM,pa);
        if (ha==0.0) continue;
        double RR = 0.0;
        for (int n=0; n<=pb.n; n++)
            for (int l=0; l<=pb.l; l++)
                for (int m=0; m<=pb.m; m++)
                {
                    Polarization NLMp(n,l,m);
                    double hb = H1b(NLMp,pb);
                    if (hb!=0.0)
                        RR += hb*R(bNLM+NLMp);
                }
        s += ha*RR;
    }
    return s * 2*Pi52*factor;
}

// <p^2>=<-nabla^2> building block (no 1/2; the Hamiltonian's).
double PrimGaussian::Grad2(const PrimGaussian* a, const PrimGaussian* b,
                           const Polarization& pa, const Polarization& pb)
{
    const Ω& ab = findΩ(a->GetGData(), b->GetGData());
    double factor = pow(Pi/ab.αₚ,1.5)*ab.Eij;
    double h = GetGrad2(pa,pb,ab);
    return h!=0 ? factor*h : 0.0;
}

double PrimGaussian::Nuclear(const PrimGaussian* a, const PrimGaussian* b,
                             const Polarization& pa, const Polarization& pb, const Structure* cl)
{
    assert(cl);
    const Ω& ab = findΩ(a->GetGData(), b->GetGData());
    RNLM R;
    for (auto atom:*cl)
        R.Add(RNLM(ab.Ltotal,ab.αₚ,ab.P-atom->itsR), -1.0*(atom->itsZ));

    auto NLMs = Ω::GetNMLs(ab.Ltotal);
    const Polarization Pab = pa + pb;
    double s = 0.0;
    for (auto bNLM:NLMs)
    {
        if (bNLM > Pab) continue;
        if (double h = ab.H2(bNLM,pa,pb); h!=0)
            s += h*R(bNLM);
    }
    return s * 2*Pi/ab.αₚ*ab.Eij;
}

// 3-centre overlap <ab|c>: the Hermite block, cached by primitive triple in the global Cache3 (the build
// logic stays here; findH3 just caches the result, no new/delete per call).
double PrimGaussian::Overlap3C(const PrimGaussian* ga, const PrimGaussian* gb, const PrimGaussian* gc,
                               const Polarization& pa, const Polarization& pb, const Polarization& pc)
{
    const Hermite3& H3 = findH3(ga, gb, gc);
    return H3(pa,pb,pc);
}

double PrimGaussian::Repulsion3C(const PrimGaussian* ga, const PrimGaussian* gb, const PrimGaussian* gc,
                                 const Polarization& pa, const Polarization& pb, const Polarization& pc)
{
    const Ω& ab(findΩ(ga->GetGData(), gb->GetGData()));
    const RNLM&        R(findRNLM(ab.GetGData(), gc->GetGData()));

    auto NLMs = Ω::GetNMLs(ab.Ltotal);
    const Hermite1& Hc = gc->GetH1();
    const Polarization Pab = pa+pb;
    double s = 0.0;
    for (auto nlm:NLMs)
    {
        if (nlm > Pab) continue;
        double hab = ab.H2(nlm,pa,pb);
        if (hab==0.0) continue;
        double Rs = 0.0;
        for (int n=0; n<=pc.n; n++)
            for (int l=0; l<=pc.l; l++)
                for (int m=0; m<=pc.m; m++)
                {
                    Polarization NLMp(n,l,m);
                    if (double h = Hc(NLMp,pc); h!=0.0)
                        Rs += h*R(nlm+NLMp);
                }
        if (Rs!=0) s += hab*Rs;
    }
    double factor = 1.0/(ab.αₚ*gc->GetExponent()*sqrt(ab.αₚ+gc->GetExponent()));
    factor = (pc.GetTotalL()%2) ? -factor : factor;
    return s * 2*Pi52 * ab.Eij*factor;
}

double PrimGaussian::Repulsion4C(const Ω& ab, const Ω& cd,
                                 const Polarization& pa, const Polarization& pb,
                                 const Polarization& pc, const Polarization& pd)
{
    double lambda = 2*Pi52/(ab.αₚ*cd.αₚ*sqrt(ab.αₚ+cd.αₚ)); //M&D 3.31
    lambda *= ab.Eij*cd.Eij; //M&D 2.25
    const RNLM& rnlm(findRNLM(ab.GetGData(), cd.GetGData())); //M&D section 4A

    // Iterate the CONTRIBUTING Hermite (N,L,M) terms directly -- those with each component <= the bra/ket
    // total (Pab/Pcd) -- instead of scanning the full GetNMLs() list and discarding most via operator>
    // (which also copied a Polarization per term and built a temporary for the RNLM index).  Same set,
    // order-independent sum.  This mirrors Repulsion3C's inner loop; it removes the hot Polarization
    // operator>/operator+/copy ops the profiler flagged.
    const Polarization Pab = pa + pb;
    const Polarization Pcd = pc + pd;
    double s = 0.0;
    for (int an=0; an<=Pab.n; ++an)
     for (int al=0; al<=Pab.l; ++al)
      for (int am=0; am<=Pab.m; ++am)
      {
          const double hab = ab.H2(Polarization(an,al,am), pa, pb);
          if (hab==0.0) continue;
          for (int cn=0; cn<=Pcd.n; ++cn)
           for (int cl=0; cl<=Pcd.l; ++cl)
            for (int cm=0; cm<=Pcd.m; ++cm)
            {
                const double hcd = cd.H2(Polarization(cn,cl,cm), pc, pd);
                if (hcd==0.0) continue;
                const double r = rnlm(Index3{an+cn, al+cl, am+cm});
                if (r!=0.0)
                {
                    const double sign = ((cn+cl+cm)&1) ? -1.0 : 1.0;   // (-1)^(total of the cd term)
                    s += hab*hcd*r*sign;
                }
            }
      }
    return s*lambda;
}

//
//  Calculate 3 center hermite functions.  Here *this is treated as the third argument, center C.
//
Hermite3* PrimGaussian::GetH3(const PrimGaussian& g1, const PrimGaussian& g2) const
{
    const double  a = g1.itsExponent, b = g2.itsExponent, c = itsExponent;
    const rvec3_t A = g1.itsCenter,   B = g2.itsCenter,   C = itsCenter;
    const int    La = g1.itsL,       Lb = g2.itsL,       Lc = itsL;

    double alphaQ = a+b+c;

    rvec3_t AB = A-B;
    rvec3_t AC = A-C;
    rvec3_t BC = B-C;
    rvec3_t Q  = (a*A+b*B+c*C)/alphaQ;

    double Eabc = pow(Pi/alphaQ,1.5)*exp( -(a*b*AB*AB + a*c*AC*AC + b*c*BC*BC) / alphaQ );

    return new Hermite3(alphaQ,Q-A,Q-B,Q-C,La,Lb,Lc,Eabc);
}

//
// double factorial table, starts at -1, so you have to add 1 to the index.
//
static constexpr double DoubleFactData[14] = {1,1,1,2,3,8,15,48,105,384,945,3840,10395,46080};
static constexpr double DoubleFact(int i) {return DoubleFactData[i+1];}

double PrimGaussian::GetNormalization(const Polarization& p) const
{
    assert(2*p.n-1 <= 12);
    assert(2*p.l-1 <= 12);
    assert(2*p.m-1 <= 12);
    double s = pow(Pi/(2*itsExponent),1.5);
    double t = intpow(4*itsExponent,p.GetTotalL());
    double f = DoubleFact(2*p.n-1) * DoubleFact(2*p.l-1) * DoubleFact(2*p.m-1);
    return sqrt(t/(s*f));
}

double PrimGaussian::GetCharge(const Polarization& p) const
{
    assert(2*p.n-1 <= 12);
    assert(2*p.l-1 <= 12);
    assert(2*p.m-1 <= 12);
    if ((p.n%2) || (p.l%2) || (p.m%2)) return 0;
    double s = pow(Pi/itsExponent,1.5);
    double t = intpow(2*itsExponent,p.GetTotalL()/2);
    double f = DoubleFact(p.n-1) * DoubleFact(p.l-1) * DoubleFact(p.m-1);
    return s*f/t;
}

double PrimGaussian::operator()(const rvec3_t& r) const
{
    rvec3_t dr = itsCenter-r;
    return exp(-itsExponent*dr*dr);
}

rvec3_t PrimGaussian::Gradient(const rvec3_t& r) const
{
    rvec3_t dr = itsCenter-r;
    return -2*itsExponent* (*this)(r) * dr;
}

//#######################################################################
//
//   GaussianRF: the single radial function (a contraction of >=1 primitives)
//

GaussianRF::GaussianRF()
    : itsCenter(0,0,0)
    , itsL     (0)
{}

GaussianRF::GaussianRF(double theExponent, const rvec3_t& theCenter, int theL)
    : itsCenter(theCenter)
    , itsL     (theL)
{
    itsPrims.push_back(std::make_unique<PrimGaussian>(theExponent,theCenter,theL));
    itsCoeff.resize(1);
    itsCoeff[0]=1.0;
}

GaussianRF::GaussianRF(const rvec_t& coeffs, const rvec_t& exponents, const rvec3_t& theCenter, int theL)
    : itsCenter(theCenter)
    , itsL     (theL)
    , itsCoeff (coeffs)
{
    assert(coeffs.size()==exponents.size());
    for (size_t i=0;i<exponents.size();++i)
        itsPrims.push_back(std::make_unique<PrimGaussian>(exponents[i],theCenter,theL));
    // Absorb each primitive's normalization into its contraction coefficient.
    for (size_t i=0;i<itsPrims.size();++i)
        itsCoeff[i] *= itsPrims[i]->GetNormalization(Polarization(theL,0,0));
}

GaussianRF::~GaussianRF() {}

rvec_t GaussianRF::GetExponents() const
{
    rvec_t es(itsPrims.size());
    for (size_t i=0;i<itsPrims.size();++i) es[i]=itsPrims[i]->GetExponent();
    return es;
}

// Deep-copy the primitives (rebuilt from exponent/centre/L, so an equal copy interns to the SAME
// registry index -- and therefore shares cache entries, unlike the old per-object UniqueID).
GaussianRF::GaussianRF(const GaussianRF& o)
    : itsCenter(o.itsCenter)
    , itsL     (o.itsL)
    , itsCoeff (o.itsCoeff)
{
    for (auto& g:o.itsPrims)
        itsPrims.push_back(std::make_unique<PrimGaussian>(g->GetExponent(),o.itsCenter,o.itsL));
}

GaussianRF& GaussianRF::operator=(const GaussianRF& o)
{
    if (this!=&o)
    {
        itsCenter = o.itsCenter;
        itsL      = o.itsL;
        itsCoeff  = o.itsCoeff;
        itsPrims.clear();
        for (auto& g:o.itsPrims)
            itsPrims.push_back(std::make_unique<PrimGaussian>(g->GetExponent(),o.itsCenter,o.itsL));
    }
    return *this;
}

// The same radial (exponents + contraction) placed at a new centre.
GaussianRF GaussianRF::AtCenter(const rvec3_t& newCenter) const
{
    GaussianRF ret;
    ret.itsCenter = newCenter;
    ret.itsL      = itsL;
    ret.itsCoeff  = itsCoeff;
    for (auto& g:itsPrims)
        ret.itsPrims.push_back(std::make_unique<PrimGaussian>(g->GetExponent(),newCenter,itsL));
    return ret;
}

bool GaussianRF::operator==(const GaussianRF& g) const
{
    if (norm(itsCenter-g.itsCenter) >= 0.01) return false;          // same centre (0.01 a.u.)
    if (itsPrims.size()!=g.itsPrims.size()) return false;
    for (size_t i=0;i<itsPrims.size();++i)
    {
        double ea=itsPrims[i]->GetExponent(), eb=g.itsPrims[i]->GetExponent();
        if (fabs((ea-eb)/(ea+eb)*2.0) >= 0.001) return false;       // exponents within 0.1%
        if (itsCoeff[i]!=g.itsCoeff[i]) return false;
    }
    return true;
}

double GaussianRF::GetNormalization(const Polarization& p) const
{
    // Only meaningful for an uncontracted function; contracted normalization comes from the
    // self-overlap (see PGData).
    assert(itsPrims.size()==1);
    return itsPrims[0]->GetNormalization(p);
}

double GaussianRF::GetCharge(const Polarization& p) const
{
    double ret = 0.0;
    for (size_t i=0;i<itsPrims.size();++i) ret += itsCoeff[i]*itsPrims[i]->GetCharge(p);
    return ret;
}

// Centre-independent identity (L + per-primitive exponent,coeff): symmetry-equivalent shells on
// different atoms share it, so the symmetry code can match them without exposing raw exponents.
std::string GaussianRF::TypeID() const
{
    std::string key = std::to_string(itsL);
    for (size_t i=0;i<itsPrims.size();++i)
        key += ":" + std::to_string(itsPrims[i]->GetExponent()) + "," + std::to_string(itsCoeff[i]);
    return key;
}

// --- 2-centre: contract the primitive pairs over the matching named PrimGaussian kernel -------------
double GaussianRF::Overlap2C(rf_t& rb, po_t& pa, po_t& pb) const
{
    double s = 0.0;
    for (size_t i=0;i<itsPrims.size();++i)
        for (size_t j=0;j<rb.itsPrims.size();++j)
            s += itsCoeff[i]*rb.itsCoeff[j]
                 * PrimGaussian::Overlap2C(itsPrims[i].get(), rb.itsPrims[j].get(), pa, pb);
    return s;
}
double GaussianRF::Repulsion2C(rf_t& rb, po_t& pa, po_t& pb) const
{
    double s = 0.0;
    for (size_t i=0;i<itsPrims.size();++i)
        for (size_t j=0;j<rb.itsPrims.size();++j)
            s += itsCoeff[i]*rb.itsCoeff[j]
                 * PrimGaussian::Repulsion2C(itsPrims[i].get(), rb.itsPrims[j].get(), pa, pb);
    return s;
}
double GaussianRF::Grad2(rf_t& rb, po_t& pa, po_t& pb) const
{
    double s = 0.0;
    for (size_t i=0;i<itsPrims.size();++i)
        for (size_t j=0;j<rb.itsPrims.size();++j)
            s += itsCoeff[i]*rb.itsCoeff[j]
                 * PrimGaussian::Grad2(itsPrims[i].get(), rb.itsPrims[j].get(), pa, pb);
    return s;
}
double GaussianRF::Nuclear(rf_t& rb, po_t& pa, po_t& pb, const Structure* cl) const
{
    double s = 0.0;
    for (size_t i=0;i<itsPrims.size();++i)
        for (size_t j=0;j<rb.itsPrims.size();++j)
            s += itsCoeff[i]*rb.itsCoeff[j]
                 * PrimGaussian::Nuclear(itsPrims[i].get(), rb.itsPrims[j].get(), pa, pb, cl);
    return s;
}

// --- 3-centre <ab|c> (this is centre C) -------------------------------------------------------------
double GaussianRF::Overlap3C(rf_t& ra, rf_t& rb, po_t& pa, po_t& pb, po_t& pc) const
{
    double s = 0.0;
    for (size_t i=0;i<ra.itsPrims.size();++i)
        for (size_t j=0;j<rb.itsPrims.size();++j)
            for (size_t k=0;k<itsPrims.size();++k)
                s += ra.itsCoeff[i]*rb.itsCoeff[j]*itsCoeff[k]
                     * PrimGaussian::Overlap3C(ra.itsPrims[i].get(), rb.itsPrims[j].get(), itsPrims[k].get(),
                                               pa, pb, pc);
    return s;
}
double GaussianRF::Repulsion3C(rf_t& ra, rf_t& rb, po_t& pa, po_t& pb, po_t& pc) const
{
    double s = 0.0;
    for (size_t i=0;i<ra.itsPrims.size();++i)
        for (size_t j=0;j<rb.itsPrims.size();++j)
            for (size_t k=0;k<itsPrims.size();++k)
                s += ra.itsCoeff[i]*rb.itsCoeff[j]*itsCoeff[k]
                     * PrimGaussian::Repulsion3C(ra.itsPrims[i].get(), rb.itsPrims[j].get(), itsPrims[k].get(),
                                                 pa, pb, pc);
    return s;
}

// --- 4-centre (ab|cd) (this is centre D) ------------------------------------------------------------
double GaussianRF::Repulsion4C(rf_t& ra, rf_t& rb, rf_t& rc, po_t& pa, po_t& pb, po_t& pc, po_t& pd) const
{
    const size_t ni=ra.itsPrims.size(), nj=rb.itsPrims.size(),
                 nk=rc.itsPrims.size(), nl=itsPrims.size();
    // Hoist the Ω (charge-distribution) lookups out of the inner contraction: ab is constant across the
    // (k,l) loop and the cd table is reused across every (i,j), so findΩ runs ni*nj+nk*nl times instead of
    // ni*nj*nk*nl (each is a cached red-black-tree lookup).  Cache entries are never evicted, so the Ω
    // addresses stay valid for the duration of this call.
    std::vector<const Ω*> cdΩ(nk*nl);
    for (size_t k=0;k<nk;++k)
        for (size_t l=0;l<nl;++l)
            cdΩ[k*nl+l] = &findΩ(rc.itsPrims[k]->GetGData(), itsPrims[l]->GetGData());

    double s = 0.0;
    for (size_t i=0;i<ni;++i)
        for (size_t j=0;j<nj;++j)
        {
            const Ω& ab = findΩ(ra.itsPrims[i]->GetGData(), rb.itsPrims[j]->GetGData());
            const double cij = ra.itsCoeff[i]*rb.itsCoeff[j];
            for (size_t k=0;k<nk;++k)
                for (size_t l=0;l<nl;++l)
                    s += cij*rc.itsCoeff[k]*itsCoeff[l]
                         * PrimGaussian::Repulsion4C(ab, *cdΩ[k*nl+l], pa, pb, pc, pd);
        }
    return s;
}

std::ostream& GaussianRF::Write(std::ostream& os) const
{
    if (itsPrims.size()==1)
        os << "Primative  " << std::setw(8) << itsPrims[0]->GetExponent();
    else
    {
        os << "Contracted {";
        for (auto& g:itsPrims) os << g->GetExponent() << " ";
        os << "}";
    }
    return os;
}

double GaussianRF::operator()(const rvec3_t& r) const
{
    double ret = 0.0;
    for (size_t i=0;i<itsPrims.size();++i) ret += itsCoeff[i]*(*itsPrims[i])(r);
    return ret;
}

rvec3_t GaussianRF::Gradient(const rvec3_t& r) const
{
    rvec3_t ret(0,0,0);
    for (size_t i=0;i<itsPrims.size();++i) ret += itsCoeff[i]*itsPrims[i]->Gradient(r);
    return ret;
}

} //namespace BasisSet::Molecule::Evaluators::PG_Cart_MnD
