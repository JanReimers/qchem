// File: CDCache.H
#ifndef _CDCache_H_
#define _CDCache_H_

#include "types.H"
#include "Misc/UniqueID/UniqueID.H"
#include <tuple>
#include <map>

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
    size_t size() const {return cache.size();}
    void Report(std::ostream&) const;
private:
    typedef std::pair<UniqueID::IDtype,UniqueID::IDtype> ids_t;
    static ids_t Sort(UniqueID::IDtype,UniqueID::IDtype);
    
    size_t CDlookups;
    std::map<ids_t,const GaussianCD*> cache;
    size_t RNLMlookups;
    std::map<ids_t,const RNLM*> RNLMcache;
    
};
#endif // _CDCache_H_
