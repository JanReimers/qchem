// File: CDCache.H
#ifndef _CDCache_H_
#define _CDCache_H_

#include <Common/UniqueID.H> 
#include <tuple>
#include <map>
#include <iosfwd>

namespace PolarizedGaussian
{
    
class GaussianCD;
class GaussianRF;
class RNLM;

class CDCache
{
public:
    CDCache();
    ~CDCache();
    const GaussianCD& find(const GaussianRF*,const GaussianRF*);
    const RNLM&       find(const GaussianCD&,const GaussianRF*);
    const RNLM&       find(const GaussianCD&,const GaussianCD&);
    const RNLM&       find(const GaussianCD&);
    
    void Report(std::ostream&) const;
private:
    static inline double Efficiency(size_t N, size_t Nl)
    {
        return 100*(1.0-N/(double)Nl);
    } 
    typedef UniqueID::IDtype id_t;
    typedef std::pair<id_t,id_t> ids_t;
    static ids_t Sort(id_t,id_t);
    
    size_t CDlookups,CDinserts;
    std::map<ids_t,const GaussianCD*> GCDcache;
    size_t RNLMlookups,RNLMinserts;
    std::map<id_t ,const RNLM*> RNLMcache1; // For 2 centres, or CD. 
    std::map<ids_t,const RNLM*> RNLMcache; // For 3 centres.
    
};

} //namespace PolarizedGaussian

#endif // _CDCache_H_
