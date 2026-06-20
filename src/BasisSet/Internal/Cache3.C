// File: Cache3.C  Cache object based on three unsigned integer indices.
//
// The 3-index sibling of Cache2/Cache4.  For molecules it caches the 3-centre Hermite block
// (GaussianH3) used by the DFT Overlap3C path, keyed by the primitive triple -- the analogue of
// Cache2 (2-centre Ω) and Cache4 (4-index atomic Slater Rk).  As with those, the string RadialType
// key is paid once (in the global cache) and the three indices are integers descended by
// loop_1/loop_2/loop_3 so the inner loop is a single-index lookup.
module;
#include <map>
#include <memory>
#include <string>
#include <iosfwd>
#include <functional>
export module qchem.BasisSet.Internal.Cache3;


export class Cache3; //fwd

export class Cache3_Client
{
public:
    virtual ~Cache3_Client() {};
    virtual std::string RadialType() const=0;
    virtual Cache3*     MakeCache3() const=0;
};

export class Cacheable3
{
public:
    virtual ~Cacheable3() {};
    virtual bool   isSupported(const Cache3_Client*) const=0;
    virtual size_t RAMsize() const=0;
};

export class Cache3
{
public:
    virtual ~Cache3();
    virtual void Register(Cache3_Client*);     //base body erases unsupported; subclass may extend

    // Facade form: re-entrant find-or-make on (i1,i2,i3).  const because the storage is mutable.
    const Cacheable3& get(size_t i1, size_t i2, size_t i3, std::function<const Cacheable3*()> make) const;

    // Descent form: loop_1/loop_2 descend, loop_3 returns/creates via Create.
    void                      loop_1(size_t i1) const;
    void                      loop_2(size_t i2) const;
    virtual const Cacheable3* loop_3(size_t i3) const;

    virtual size_t            RAMsize() const;
    virtual const Cacheable3* Create(size_t i1,size_t i2,size_t i3) const; //default: unused by facade

    // Hit/miss stats (see Cache2).  Report takes the cache's name (its RadialType key).
    size_t Lookups() const {return itsLookups;}
    size_t Inserts() const {return itsInserts;}
    void   Report(std::ostream&, const std::string& name) const;
private:

    typedef std::map<size_t,std::unique_ptr<const Cacheable3>> cache_3;
    typedef std::map<size_t,cache_3> cache_2;
    typedef std::map<size_t,cache_2> cache_t;

    mutable cache_t  cache;
    mutable cache_2* i1_cache;
    mutable cache_3* i2_cache;
    mutable size_t   i1,i2,i3; //Current indexes
    mutable size_t   itsLookups=0, itsInserts=0; //hit/miss stats
};
