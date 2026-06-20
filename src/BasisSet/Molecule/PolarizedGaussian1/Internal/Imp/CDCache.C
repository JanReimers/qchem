// File: CDCache.C
module;
#include <string>
module qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.CDCache;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.MnD.RNLM;
import qchem.BasisSet.Internal.DB_Cache;   // theGlobalCache, Register/GetCache2
import qchem.BasisSet.Internal.Cache2;     // Cache2, Cacheable2, Cache2_Client
namespace BasisSet::Molecule::PolarizedGaussian1
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
        std::string RadialType() const override {return "PG1.Omega";}
        Cache2*     MakeCache2() const override {return new Cache2;}
    };
    const Cache2* OmegaCache()
    {
        static bool reg = []{ static OmegaClient c; BasisSet::theGlobalCache->Register(&c); return true; }();
        (void)reg;
        return BasisSet::theGlobalCache->GetCache2("PG1.Omega");
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
        std::string RadialType() const override {return "PG1.RNLM";}
        Cache2*     MakeCache2() const override {return new Cache2;}
    };
    const Cache2* RNLMCache()
    {
        static bool reg = []{ static RNLMClient c; BasisSet::theGlobalCache->Register(&c); return true; }();
        (void)reg;
        return BasisSet::theGlobalCache->GetCache2("PG1.RNLM");
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

} //namespace BasisSet::Molecule::PolarizedGaussian1
