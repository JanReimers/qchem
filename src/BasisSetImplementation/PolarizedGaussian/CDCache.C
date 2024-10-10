// File: CDCache.H

#include "CDCache.H"
#include "Gaussian/GaussianCD.H"
#include "Gaussian/GaussianRF.H"
#include "Auxillary/RNLM.H"


CDCache::CDCache() : CDlookups(0), RNLMlookups(0) {};

CDCache::~CDCache()
{
    for (auto c:cache) delete c.second;
    for (auto r:RNLMcache1) delete r.second;
    for (auto r:RNLMcache) delete r.second;
}

void CDCache::Report(std::ostream& os) const
{
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
    CDlookups++;
    ids_t key=std::make_pair(a->GetID(),b->GetID());
    if (auto i=cache.find(key);i==cache.end())
        return *(cache[key]=new GaussianCD(*a,*b));
    else
        return *(i->second);
}

const RNLM& CDCache::find(const GaussianCD& ab,const GaussianRF* c)
{
    RNLMlookups++;
    ids_t key=std::make_pair(ab.GetID(),c->GetID());
    if (auto i=RNLMcache.find(key);i==RNLMcache.end())
    {
        double alpha =ab.AlphaP*c->itsExponent/(ab.AlphaP+c->itsExponent);
        return *(RNLMcache[key]=new RNLM(ab.Ltotal+c->GetL(),alpha,ab.P-c->GetCenter()));
    }
    else
        return *(i->second);
}

const RNLM& CDCache::find(const GaussianCD& ab,const GaussianCD& cd)
{
    RNLMlookups++;
    ids_t key=std::make_pair(ab.GetID(),cd.GetID());
    if (auto i=RNLMcache.find(key);i==RNLMcache.end())
    {
        double alpha=ab.AlphaP*cd.AlphaP/(ab.AlphaP+cd.AlphaP); //M&D 3.32
        RVec3 PQ = ab.P-cd.P; //M&D 3.32
        return *(RNLMcache[key]=new RNLM(ab.Ltotal+cd.Ltotal,alpha,PQ)) ;
    }
    else
        return *(i->second);
}

const RNLM& CDCache::find(const GaussianCD& ab)
{
    RNLMlookups++;
    id_t key=ab.GetID();
    if (auto i=RNLMcache1.find(key);i==RNLMcache1.end())
        return *(RNLMcache1[key]=new RNLM(ab.Ltotal,ab.ab/ab.AlphaP,ab.AB));
    else
        return *(i->second);
}
