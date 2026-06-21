// File: Omega.C  Charge distribution Ω for a primitive pair + the global Ω/RNLM/H3 Cache access points.
module;
#include <string>
#include <vector>
#include <memory>
#include <functional>
module qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Internal.Omega;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Polarization;
import qchem.BasisSet.Molecule.Evaluators.Internal.MnD;   // generic MnD core (RNLM, ...)
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Internal.Hermite3;
import qchem.BasisSet.Internal.DB_Cache;   // theGlobalCache, Register/GetCache{2,3}
import qchem.BasisSet.Internal.Cache2;     // Cache2, Cacheable2, Cache2_Client
import qchem.BasisSet.Internal.Cache3;     // Cache3, Cacheable3, Cache3_Client
import qchem.Math;                          // exp (Ω ctor)
namespace BasisSet::Molecule::Evaluators::PG_Cart_MnD
{

// All PG charge-distribution caching lives in process-global Cache2s, keyed by UniqueIDs so entries
// are shared wherever the underlying objects are shared (e.g. SALC irreps over one raw basis):
//   - Ω (charge distribution) keyed by the primitive pair (primA.ID, primB.ID)
//   - 3/4-centre RNLM keyed by the two charge-distribution ids (ab.ID, c.ID)
// (The 2-centre self-RNLM is a lazy member of each Ω, not cached here.)
namespace
{
    struct OmegaClient : Cache2_Client
    {
        std::string RadialType() const override {return "PG.Omega";}
        Cache2*     MakeCache2() const override {return new Cache2;}
    };
    const Cache2* OmegaCache()
    {
        static bool reg = []{ static OmegaClient c; BasisSet::theGlobalCache->Register(&c); return true; }();
        (void)reg;
        return BasisSet::theGlobalCache->GetCache2("PG.Omega");
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
    const Cache2* RNLMCache()
    {
        static bool reg = []{ static RNLMClient c; BasisSet::theGlobalCache->Register(&c); return true; }();
        (void)reg;
        return BasisSet::theGlobalCache->GetCache2("PG.RNLM");
    }

    // Cacheable3 wrapper holding the 3-centre Hermite block (GaussianH3, via the Hermite3 base).
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
    const Cache3* H3Cache()
    {
        static bool reg = []{ static H3Client c; BasisSet::theGlobalCache->Register(&c); return true; }();
        (void)reg;
        return BasisSet::theGlobalCache->GetCache3("PG.H3");
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
            double alpha = ab.Alpha*c.Alpha/(ab.Alpha+c.Alpha);
            return new RNLM_C2(ab.L+c.L, alpha, ab.R-c.R);
        });
    return static_cast<const RNLM_C2&>(w).rnlm;
}

const Hermite3& findH3(UniqueID::IDtype a, UniqueID::IDtype b, UniqueID::IDtype c,
                       std::function<Hermite3*()> make)
{
    const Cacheable3& w = H3Cache()->get(a, b, c,
        [&make]() -> const Cacheable3* { return new H3_C3(make()); });
    return *static_cast<const H3_C3&>(w).h3;
}

//------------------------------------------------------------------------------------------------
//  Ω (charge distribution for a primitive pair) methods.
//------------------------------------------------------------------------------------------------
Ω::~Ω()
{
    delete itsSelfRNLM;
}

// 2-centre self-auxiliary: RNLM(Ltotal, ab/AlphaP, AB) -- built once, owned by this Ω.
const RNLM& Ω::SelfRNLM() const
{
    if (!itsSelfRNLM) itsSelfRNLM = new RNLM(Ltotal, ab/AlphaP, AB);
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
    : Ltotal(g1.L + g2.L)
    , a     (g1.Alpha)
    , b     (g2.Alpha)
    , ab    (a * b)
    , AlphaP(a + b)
    , AB    (g1.R - g2.R)
    , P     ( (a*g1.R + b*g2.R) / AlphaP)
    , Eij   ( exp(-ab / AlphaP * (AB*AB)) )
    , H2    (AlphaP, P - g1.R, P - g2.R, g1.L+1, g2.L+1)
{
    if (theNMLs.size()==0) MakeNMLs();
}

size_t Ω::RAMsize() const
{
    return sizeof(Ω); // includes the by-value Hermite2 block
}

} //namespace BasisSet::Molecule::Evaluators::PG_Cart_MnD
