// File: Cache4.H Cache object based on four unsigned integer indices.
#ifndef _Cache4_H_
#define _Cache4_H_

#include <map>
//
//  Abstract base for the type being cached.
//
class Cacheable
{
public:
    virtual ~Cacheable() {};
};
//
//  Cache object based on four unsigned integer indices.  Client code using the caching
//  should derive rom this class and overload all the protected functions with function forwarding.
//  Use covariant return return types for the loop_4 overload. 
//  Derived class also needs to supply a Create function.
//
class Cache4
{
protected:
    ~Cache4();
    void       loop_1(size_t i1) const;
    void       loop_2(size_t i2) const;
    void       loop_3(size_t i3) const;
    virtual const Cacheable* loop_4(size_t i4) const;
    
private:
    virtual const Cacheable* Create(size_t i1,size_t i2,size_t i3,size_t i4) const=0;

    typedef std::map<size_t,const Cacheable*> cache_4; 
    typedef std::map<size_t,cache_4> cache_3; 
    typedef std::map<size_t,cache_3> cache_2; 
    typedef std::map<size_t,cache_2> cache_t; 
    
    mutable cache_t cache;
    mutable cache_2* i1_cache;
    mutable cache_3* i2_cache;
    mutable cache_4* i3_cache;
    mutable size_t i1,i2,i3,i4; //Current indexes
};

#endif // _Cache4_H_
