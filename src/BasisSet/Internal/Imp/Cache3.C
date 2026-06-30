// File: Cache3Imp.C  Cache object based on three unsigned integer indices.
module;
#include <map>
#include <memory>
#include <cassert>
#include <functional>
#include <iostream>
#include <iomanip>
#include <string>
module qchem.BasisSet.Internal.Cache3;

namespace qchem {


Cache3::~Cache3()
{
}

void Cache3::Register(Cache3_Client* client)
{
    for (auto& a:cache)
        for (auto& b:a.second)
            std::erase_if(b.second, [client](const auto& pair)
            {
                return pair.second->isSupported(client);
            });
}

size_t Cache3::RAMsize() const
{
    size_t ram=0;
    for (auto& a:cache)
        for (auto& b:a.second)
            for (auto& c:b.second)
                ram+=c.second->RAMsize();
    return ram;
}

void Cache3::Report(std::ostream& os, const std::string& name) const
{
    double reuse = itsLookups ? 100.0*(1.0 - double(itsInserts)/double(itsLookups)) : 0.0;
    os << "    " << std::left << std::setw(12) << name << std::right
       << " entries=" << std::setw(9) << itsInserts
       << " lookups=" << std::setw(10) << itsLookups
       << " reuse="   << std::setw(6) << std::setprecision(4) << reuse << "%" << std::endl;
}

const Cacheable3* Cache3::Create(size_t,size_t,size_t) const
{
    return nullptr; // facade (get) form does not use Create; override for the loop_3 descent form
}

void Cache3::loop_1(size_t _i1) const
{
    i1=_i1;
    if (auto i=cache.find(i1);i==cache.end())
        i1_cache = &(cache[i1]=cache_2());
    else
        i1_cache = &i->second;
}

void Cache3::loop_2(size_t _i2) const
{
    i2=_i2;
    cache_2& c=*i1_cache;
    if (auto i=c.find(i2);i==c.end())
        i2_cache = &(c[i2]=cache_3());
    else
        i2_cache = &i->second;
}

const Cacheable3* Cache3::loop_3(size_t _i3) const
{
    i3=_i3;
    ++itsLookups;
    cache_3& c=*i2_cache;
    auto i=c.find(i3);
    if (i==c.end())
    {
        ++itsInserts;
        const auto [iterator,success]=c.insert({i3,std::unique_ptr<const Cacheable3>(Create(i1,i2,i3))});
        assert(success);
        i=iterator;
    }
    return i->second.get();
}

} // namespace qchem