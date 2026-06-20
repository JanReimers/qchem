// File: Cache2Imp.C  Cache object based on two unsigned integer indices.
module;
#include <map>
#include <memory>
#include <cassert>
#include <functional>
#include <iostream>
#include <iomanip>
#include <string>
module qchem.BasisSet.Internal.Cache2;


//---------------------------------------------------------------------------------------

Cache2::~Cache2()
{
}

//  All unsupported entries are removed.  They are automatically recreated next time loop_2 is called.
//  Subclasses override Register to assign their integer indices, then delegate to this base body.
void Cache2::Register(Cache2_Client* client)
{
    for (auto& a:cache)
        std::erase_if(a.second, [client](const auto& pair) //C++-20 magic!!
        {
            return pair.second->isSupported(client);
        });
}

size_t Cache2::RAMsize() const
{
    size_t ram=0;
    for (auto& a:cache)
        for (auto& b:a.second)
            ram+=b.second->RAMsize();
    return ram;
}

const Cacheable2& Cache2::get(size_t i1, size_t i2, std::function<const Cacheable2*()> make) const
{
    ++itsLookups;
    cache_2& sub = cache[i1];                 // creates an empty sub-map on first use
    auto it = sub.find(i2);
    if (it==sub.end())
    {
        ++itsInserts;
        const auto [iterator,success]=sub.insert({i2,std::unique_ptr<const Cacheable2>(make())});
        assert(success);
        it=iterator;
    }
    return *it->second;
}

void Cache2::Report(std::ostream& os, const std::string& name) const
{
    double reuse = itsLookups ? 100.0*(1.0 - double(itsInserts)/double(itsLookups)) : 0.0;
    os << "    " << std::left << std::setw(12) << name << std::right
       << " entries=" << std::setw(9) << itsInserts
       << " lookups=" << std::setw(10) << itsLookups
       << " reuse="   << std::setw(6) << std::setprecision(4) << reuse << "%" << std::endl;
}

const Cacheable2* Cache2::Create(size_t,size_t) const
{
    return nullptr; // facade (get) form does not use Create; override for the loop_2 descent form
}

void Cache2::loop_1(size_t _i1) const
{
    i1=_i1;
    if (auto i=cache.find(i1);i==cache.end())
        i1_cache = &(cache[i1]=cache_2());
    else
        i1_cache = &i->second;
}

const Cacheable2* Cache2::loop_2(size_t _i2) const
{
    i2=_i2;
    ++itsLookups;
    cache_2& c=*i1_cache; //De-reference for readability.
    auto i=c.find(i2);
    if (i==c.end())
    {
        ++itsInserts;
        const auto [iterator,success]=c.insert({i2,std::unique_ptr<const Cacheable2>(Create(i1,i2))});
        assert(success);
        i=iterator;
    }
    return i->second.get();
}
