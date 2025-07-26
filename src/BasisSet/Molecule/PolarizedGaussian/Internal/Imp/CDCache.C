// File: CDCache.C
module;
#include <iomanip>
#include "PolarizedGaussian/MnD/RNLM.H"
module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.CDCache;
import Common.UniqueID; 
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

const GaussianCD& CDCache::findCD(const GData& a,const GData& b)
{
    CDlookups++;
    ids_t key=std::make_pair(a.ID,b.ID);
    if (auto i=GCDcache.find(key);i==GCDcache.end())
    {
        CDinserts++;
        return *(GCDcache[key]=new GaussianCD(a,b));
    }
    else
        return *(i->second);
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
