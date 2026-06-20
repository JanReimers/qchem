// File: CDCache.C
module;
#include <iomanip>
module qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.CDCache;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.MnD.RNLM;
import qchem.BasisSet.Internal.DB_Cache;   // theGlobalCache, Register/GetCache2
import qchem.BasisSet.Internal.Cache2;     // Cache2, Cacheable2, Cache2_Client
import Common.UniqueID;
namespace BasisSet::Molecule::PolarizedGaussian1
{

// Omega_ab (Ω) lives in ONE process-global Cache2, keyed by the primitive UniqueID pair
// (primA.ID, primB.ID).  Because the key is the primitive identity, Omega is shared automatically
// wherever primitive objects are shared (e.g. SALC irreps over one raw basis).  Registered once.
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
}

CDCache::CDCache() : CDlookups(0), CDinserts(0), RNLMlookups(0), RNLMinserts(0) {};

CDCache::~CDCache()
{
    // Omega/Ω are owned by the global Cache2 now; only the RNLM caches are local.
    for (auto r:RNLMcache1) delete r.second;
    for (auto r:RNLMcache) delete r.second;
}

using std::setw;
void CDCache::Report(std::ostream& os) const
{
    os.precision(4);
    {
        double eff=Efficiency(CDinserts,CDlookups);
        os << "    Charge Distributions cache N=" << setw(10) << CDinserts << " lookups=" << setw(10) << CDlookups << " efficiencty=" << eff << "%" << std::endl;
    }
    {
        double eff=Efficiency(RNLMinserts,RNLMlookups);
        os << "    RNLM                 cache N=" << setw(10) << RNLMinserts << " lookups=" << setw(10) << RNLMlookups << " efficiencty=" << eff << "%" << std::endl;
    }
}

CDCache::ids_t CDCache::Sort(UniqueID::IDtype i1,UniqueID::IDtype i2)
{
    return i1<=i2 ? std::make_pair(i1,i2) : std::make_pair(i2,i1);
}

const Ω& CDCache::findCD(const GData& a,const GData& b)
{
    CDlookups++;
    // Delegate to the process-global Cache2 (keyed by the primitive UniqueID pair).  On a miss it
    // builds the Ω via the lambda; the Cache2 owns the stored object.
    const Cacheable2& cd = OmegaCache()->get(a.ID, b.ID,
        [&a,&b]() -> const Cacheable2* { return new Ω(a,b); });
    return static_cast<const Ω&>(cd);
}

const RNLM& CDCache::find(const GData& ab,const GData& c)
{
    RNLMlookups++;
    ids_t key=std::make_pair(ab.ID,c.ID);
    if (auto i=RNLMcache.find(key);i==RNLMcache.end())
    {
        RNLMinserts++;
        double alpha =ab.Alpha*c.Alpha/(ab.Alpha+c.Alpha);
        return *(RNLMcache[key]=new RNLM(ab.L+c.L,alpha,ab.R-c.R));
    }
    else
        return *(i->second);
}

const RNLM& CDCache::find(const Ω& ab)
{
    RNLMlookups++;
    id_t key=ab.GetID();
    if (auto i=RNLMcache1.find(key);i==RNLMcache1.end())
    {
        RNLMinserts++;
        return *(RNLMcache1[key]=new RNLM(ab.Ltotal,ab.ab/ab.AlphaP,ab.AB));
    }
    else
        return *(i->second);
}

} //namespace BasisSet::Molecule::PolarizedGaussian1
