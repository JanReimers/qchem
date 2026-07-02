// File: Cache4Imp.C Cache object based on four unsigned integer indices.
module;
#include <map>
#include <memory>
#include <iostream>
#include <iomanip>
#include <string>
#include <cassert>
module qchem.BasisSet.Internal.Cache4;

namespace qchem {


//---------------------------------------------------------------------------------------

Cache4::~Cache4()
{
    //  for (auto a:cache) 
    //     for (auto b:a.second) 
    //         for (auto c:b.second) 
    //             for (auto d:c.second) delete d.second;
}

//  When a client registers it may raise the maximum angular momentum of some shared exponent (grouper
//  bumps maxl in the sibling Register call that runs just before this one).  Any already-cached Rk whose
//  stored LMax is now too small for this client is STALE and must be dropped so loop_4 recreates it with
//  the (now larger) grouper LMax.  We evict exactly the UNsupported entries -- isSupported(client) is
//  (client.l <= Rk.LMax), so we keep those and erase the rest.
//
//  NB: erasing on the WRONG side here (evicting the supported entries) silently wipes the whole cache on
//  every new basis -- the registering basis's l=0 shell is <= every LMax -- destroying all cross-element
//  Rk reuse (a heavy atom can no longer warm the pool for its lighter, subset-exponent siblings).
void Cache4::Register(Cache4_Client* client)
{
    for (auto& a:cache)
        for (auto& b:a.second)
            for (auto& c:b.second)
                std::erase_if(c.second, [client](const auto& pair) //C++-20 magic!!
                {
                    return !pair.second->isSupported(client);
                });
}

size_t Cache4::RAMsize() const
{
    size_t ram=0;
    for (auto& a:cache) 
        for (auto& b:a.second) 
            for (auto& c:b.second) 
                for (auto& d:c.second)
                    ram+=d.second->RAMsize();
    return ram;
}

void Cache4::loop_1(size_t _i1) const
{
    i1=_i1;
    if (auto i=cache.find(i1);i==cache.end())
        i1_cache = &(cache[i1]=cache_2());
    else
        i1_cache = &i->second;
}

void Cache4::loop_2(size_t _i2) const
{
    i2=_i2;
    cache_2& c=*i1_cache; //De-reference for readability.
    if (auto i=c.find(i2);i==c.end())
        i2_cache = &(c[i2]=cache_3());
    else
        i2_cache = &i->second;
}
void Cache4::loop_3(size_t _i3) const
{
    i3=_i3;
    cache_3& c=*i2_cache; //De-reference for readability.
    if (auto i=c.find(i3);i==c.end())
        i3_cache = &(c[i3]=cache_4());
    else
        i3_cache = &i->second;    
}

const Cacheable4* Cache4::loop_4(size_t _i4)  const
{
    i4=_i4;
    ++itsLookups;
    cache_4& c=*i3_cache; //De-reference for readability.
    auto i=c.find(i4);
    if (i==c.end())
    {
        ++itsInserts;
        const auto [iterator,success]=c.insert({i4,std::unique_ptr<const Cacheable4>(Create(i1,i2,i3,i4))});
        assert(success);
        i=iterator;
    }
    return i->second.get();
}

void Cache4::Report(std::ostream& os, const std::string& name) const
{
    double reuse = itsLookups ? 100.0*(1.0 - double(itsInserts)/double(itsLookups)) : 0.0;
    os << "    " << std::left << std::setw(12) << name << std::right
       << " entries=" << std::setw(9) << itsInserts
       << " lookups=" << std::setw(10) << itsLookups
       << " reuse="   << std::setw(6) << std::setprecision(4) << reuse << "%" << std::endl;
}

} // namespace qchem