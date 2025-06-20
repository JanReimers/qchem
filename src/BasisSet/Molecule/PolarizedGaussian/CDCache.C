// File: CDCache.H

#include "PolarizedGaussian/CDCache.H"
#include "PolarizedGaussian/Radial/GaussianCD.H"
#include "PolarizedGaussian/Radial/GaussianRF.H"
#include "PolarizedGaussian/MnD/RNLM.H"
#include <iomanip>

namespace PolarizedGaussian
{

CDCache::CDCache() : CDlookups(0), CDinserts(0), RNLMlookups(0), RNLMinserts(0) {};

CDCache::~CDCache()
{
    for (auto c:GCDcache) delete c.second;
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

const GaussianCD& CDCache::find(const GaussianRF* a,const GaussianRF* b)
{
    assert(a);
    assert(b);
    CDlookups++;
    ids_t key=std::make_pair(a->GetID(),b->GetID());
    if (auto i=GCDcache.find(key);i==GCDcache.end())
    {
        CDinserts++;
        return *(GCDcache[key]=new GaussianCD(*a,*b));
    }
    else
        return *(i->second);
}

const RNLM& CDCache::find(const GaussianCD& ab,const GaussianRF* c)
{
    RNLMlookups++;
    ids_t key=std::make_pair(ab.GetID(),c->GetID());
    if (auto i=RNLMcache.find(key);i==RNLMcache.end())
    {
        RNLMinserts++;
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
        RNLMinserts++;
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
    {
        RNLMinserts++;
        return *(RNLMcache1[key]=new RNLM(ab.Ltotal,ab.ab/ab.AlphaP,ab.AB));
    }
    else
        return *(i->second);
}

} //namespace PolarizedGaussian
