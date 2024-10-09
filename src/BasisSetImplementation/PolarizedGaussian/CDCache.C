// File: CDCache.H

#include "CDCache.H"
#include "Gaussian/GaussianCD.H"
#include "Gaussian/GaussianRF.H"
#include "Auxillary/RNLM.H"


CDCache::CDCache() : CDlookups(0), RNLMlookups(0) {};

CDCache::~CDCache()
{
    for (auto c:cache) delete c.second;
}

void CDCache::Report(std::ostream& os) const
{
    for (auto c:cache) c.second->Report(os);
    double eff=(100.0*CDlookups)/(size()+CDlookups-1);
    os << "Charge Distributions cache N=" << size() << " lookups=" << CDlookups << " efficiencty=" << eff << "%" << std::endl;
}

CDCache::ids_t CDCache::Sort(UniqueID::IDtype i1,UniqueID::IDtype i2)
{
    return i1<=i2 ? std::make_pair(i1,i2) : std::make_pair(i2,i1);
}

const GaussianCD& CDCache::find(const GaussianRF* a,const GaussianRF* b)
{
    assert(a);
    assert(b);
    ids_t key=std::make_pair(a->GetID(),b->GetID());
//    ids_t key=Sort(a->GetID(),b->GetID());
    auto i=cache.find(key);
    if (i==cache.end())
    {
        cache[key]=new GaussianCD(*a,*b);
        i=cache.find(key);
    }
    CDlookups++;
    return *(i->second);
}

const RNLM& CDCache::find(const GaussianCD& ab,const GaussianRF* c)
{
    ids_t key=std::make_pair(ab.GetID(),c->GetID());
//    ids_t key=Sort(a->GetID(),b->GetID());
    auto i=RNLMcache.find(key);
    if (i==RNLMcache.end())
    {
        double alpha =ab.AlphaP*c->itsExponent/(ab.AlphaP+c->itsExponent);
        RNLMcache[key]=new RNLM(ab.Ltotal+c->GetL(),alpha,ab.P-c->GetCenter());
        i=RNLMcache.find(key);
    }
    RNLMlookups++;
    return *(i->second);
 
}

const RNLM& CDCache::find(const GaussianCD& ab,const GaussianCD& cd)
{
    ids_t key=std::make_pair(ab.GetID(),cd.GetID());
//    ids_t key=Sort(a->GetID(),b->GetID());
    auto i=RNLMcache.find(key);
    if (i==RNLMcache.end())
    {
        double alpha=ab.AlphaP*cd.AlphaP/(ab.AlphaP+cd.AlphaP); //M&D 3.32
        RVec3 PQ = ab.P-cd.P; //M&D 3.32
        RNLMcache[key]=new RNLM(ab.Ltotal+cd.Ltotal,alpha,PQ);
        i=RNLMcache.find(key);
    }
    RNLMlookups++;
    return *(i->second);
 
}
