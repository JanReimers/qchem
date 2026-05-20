// File: Cache4.C Cache object based on four unsigned integer indices.
module;
#include <map>
#include <string>
export module qchem.BasisSet.Internal.Cache4;


export class Cache4; //fwd

//
//  Derive from this class if you know how to create a Cache4 derived object.
//
export class Cache4_Client
{
public:
    virtual ~Cache4_Client() {};
    virtual std::string RadialType() const=0; //Used as the map key in the database.
    virtual Cache4*    MakeCache4() const=0;
};
//
//  Abstract base for the type being cached.  Usually these are the Slater integrals Rk(abcd)
//
export class Cacheable
{
public:
    virtual ~Cacheable() {};
    virtual bool isSupported(const Cache4_Client*) const=0;
};

//
//  Derive from this class if you need to run four index loops in order to calculation HF Direct and Exchange
//  tables.  Cached object is based on four unsigned integer indices.  
//  Use covariant return types for the loop_4 overload. 
//
export class Cache4
{
public: 
    virtual ~Cache4();
    virtual void Register(Cache4_Client*)=0;

    void       loop_1(size_t i1) const;
    void       loop_2(size_t i2) const;
    void       loop_3(size_t i3) const;
    virtual const Cacheable* loop_4(size_t i4) const;
//
//     template <class T> const T* Tloop_4(size_t i4) const
//     {
//         const Cacheable1* c=loop_4(i4);
//         assert(c);
//         const T* Tc=dynamic_cast<const T*>(c);
//         assert(Tc);
//         return Tc;
//     } 
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