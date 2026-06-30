// File: Cache2.C  Cache object based on two unsigned integer indices.
//
// The molecular analogue of Cache4.  Atoms cache the expensive 4-index Slater integrals Rk (Cache4);
// molecules cache the expensive 2-index charge distributions Omega_ab -- the Hermite expansion for a
// primitive pair.  2-centre integrals are trivial for atoms (the angular part factors out) but not
// for molecules, which is why molecules want a Cache2 and atoms never did.
//
// As with Cache4, the string RadialType key is paid once (in the global cache, keyed by RadialType);
// within a Cache2 the two indices are plain integers descended by loop_1/loop_2, so the inner loop is
// a single-index lookup -- no per-pair string keys.
module;
#include <map>
#include <memory>
#include <string>
#include <iosfwd>
#include <functional>
export module qchem.BasisSet.Internal.Cache2;

namespace qchem {


export class Cache2; //fwd

//
//  Derive from this class if you know how to create a Cache2 derived object.
//
export class Cache2_Client
{
public:
    virtual ~Cache2_Client() {};
    virtual std::string RadialType() const=0; //Used as the map key in the database.
    virtual Cache2*     MakeCache2() const=0;
};
//
//  Abstract base for the type being cached.  For molecules these are the charge distributions
//  Omega_ab (GaussianCD).
//
export class Cacheable2
{
public:
    virtual ~Cacheable2() {};
    virtual bool   isSupported(const Cache2_Client*) const=0;
    virtual size_t RAMsize() const=0;
};

//
//  Derive from this class if you need to run two index loops to build 2-centre data.  Cached objects
//  are keyed on two unsigned integer indices.  The two index lookup is executed as each index is
//  looped over (loop_1 descends, loop_2 returns/creates), leaving a one index lookup in the inner
//  loop.  Use covariant return types for the loop_2 overload.
//
export class Cache2
{
public:
    virtual ~Cache2();
    virtual void Register(Cache2_Client*);     //base body erases unsupported; subclass may extend

    // Facade form: re-entrant find-or-make on (i1,i2).  On a miss make() is invoked (it may perform
    // nested cached gets on OTHER Cache2s), the result is stored and a reference to the STORED object
    // is returned.  Use this when the caller already has the data to build the entry.  const because
    // the storage is mutable -- so a `const Cache2*` from GetCache2 can still populate.
    //
    // make is a TEMPLATE parameter (not std::function): findΩ/findRNLM call get() from the innermost
    // contraction loops, and the cache is ~99.9% hits, so a std::function would heap/construct a wrapper
    // every call only to (almost) never invoke it -- profiling put that at several % of MnD integral time.
    template <class Make>
    const Cacheable2& get(size_t i1, size_t i2, Make&& make) const
    {
        ++itsLookups;
        cache_2& sub = cache[i1];                  // outer map: find-or-create the i1 sub-map
        auto it = sub.find(i2);
        if (it==sub.end())
        {
            ++itsInserts;
            it = sub.emplace(i2, std::unique_ptr<const Cacheable2>(make())).first;
        }
        return *it->second;
    }

    // Descent form (for hot loops): loop_1 descends to the i1 sub-map, loop_2 returns/creates via
    // Create.  Override Create in a subclass that knows how to build from indices alone.
    void                      loop_1(size_t i1) const;
    virtual const Cacheable2* loop_2(size_t i2) const;

    virtual size_t            RAMsize() const;  //Optional override
    virtual const Cacheable2* Create(size_t i1,size_t i2) const; //default: unused by the facade form

    // Hit/miss stats (high reuse = the cache is shared, e.g. Omega across SALC irreps).  Report takes
    // the cache's name (its RadialType key, which the cache itself does not store).
    size_t Lookups() const {return itsLookups;}
    size_t Inserts() const {return itsInserts;}
    void   Report(std::ostream&, const std::string& name) const;
private:

    typedef std::map<size_t,std::unique_ptr<const Cacheable2>> cache_2;
    typedef std::map<size_t,cache_2> cache_t;

    mutable cache_t  cache;
    mutable cache_2* i1_cache;
    mutable size_t   i1,i2; //Current indexes
    mutable size_t   itsLookups=0, itsInserts=0; //hit/miss stats
};

} // namespace qchem