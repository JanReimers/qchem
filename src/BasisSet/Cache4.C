// File: Cache4.H Cache object based on four unsigned integer indices.

#include "Cache4.H"

Cache4::~Cache4()
{
     for (auto a:cache) 
        for (auto b:a.second) 
            for (auto c:b.second) 
                for (auto d:c.second) delete d.second;
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

const Cacheable* Cache4::loop_4(size_t _i4)  const
{
    i4=_i4;
    cache_4& c=*i3_cache; //De-reference for readability.
    if (auto i=c.find(i4);i==c.end())
        return c[i4]=Create(i1,i2,i3,i4);
    else
        return i->second;    
}
