// File: Cache4.C Cache object based on four unsigned integer indices.
module;
#include <stdexcept>   // RequireEngine (the survey fix, 2026-09-08)
#include <string>
#include <map>
#include <memory>
#include <string>
#include <iosfwd>
export module qchem.BasisSet.Internal.Cache4;

namespace qchem {


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
export class Cacheable4
{
public:
    virtual ~Cacheable4() {};
    virtual bool isSupported(const Cache4_Client*) const=0;
    virtual size_t RAMsize() const=0;
};

//! \brief Cross-cast a \c Cacheable4 to the engine type that made it, THROWING with the type names when it
//! is not that type.
//!
//! WHY IT EXISTS (the \c dynamic_cast survey, 2026-09-08; CLAUDE.md's standing TODO: *"throw custom
//! exceptions full of relevant information in the event they fail"*).  Six call sites -- the
//! \c direct/\c exchange function-pointer pair in each of the Gaussian, Slater and BSpline atom
//! evaluators -- did this:
//! \code
//!     const RkEngine* cd = dynamic_cast<const RkEngine*>(c);
//!     return cd->DirectRk(la,lc,Ak);          // <-- no null check
//! \endcode
//! The cast is correct BY CONSTRUCTION (the \c Cacheable4 came from the same evaluator's \c MakeCache4),
//! but nothing enforced that, and the failure mode was a null dereference: a SEGFAULT in a static callback
//! several layers below the mistake, with no evaluator, no engine and no operand in the message.  This
//! turns that into a named exception, at the boundary, in Release as well as Debug -- an \c assert would be
//! compiled out exactly where the benchmark runs (the same lesson \c RequireSiteBlocks records in
//! src/Structure/Imp/UnitCell.C).
//!
//! \a who is the caller's own name, so the message says which engine and which operand.
export template <class E> const E& RequireEngine(const Cacheable4* c, const char* who)
{
    if (const E* e=dynamic_cast<const E*>(c)) return *e;
    throw std::runtime_error(std::string(who)+": this Cacheable4 is not the engine type this evaluator "
        "makes.  The cache entry and the evaluator that reads it must come from the same radial family -- "
        "a mismatch means a Cache4 built by one evaluator was handed to another, which is a composition "
        "error, not a recoverable condition.");
}

//
//  Derive from this class if you need to run four index loops in order to calculation HF Direct and Exchange
//  tables.  Cached object is based on four unsigned integer indices. Use covariant return types for the loop_4 overload. 
//  The four index lookup is executed as each index is looped over, leaving just a one index lookup in the innermost loop.
//
export class Cache4
{
public: 
    virtual ~Cache4();
    virtual void Register(Cache4_Client*)=0;

    void       loop_1(size_t i1) const;
    void       loop_2(size_t i2) const;
    void       loop_3(size_t i3) const;
    virtual const Cacheable4* loop_4(size_t i4) const;
//
//     template <class T> const T* Tloop_4(size_t i4) const
//     {
//         const Cacheable41* c=loop_4(i4);
//         assert(c);
//         const T* Tc=dynamic_cast<const T*>(c);
//         assert(Tc);
//         return Tc;
//     } 
    virtual size_t RAMsize() const; //Optional override
    virtual const Cacheable4* Create(size_t i1,size_t i2,size_t i3,size_t i4) const=0;

    // Hit/miss stats (see Cache2).  Report takes the cache's name (its RadialType key).
    size_t Lookups() const {return itsLookups;}
    size_t Inserts() const {return itsInserts;}
    void   Report(std::ostream&, const std::string& name) const;
private:

    typedef std::map<size_t,std::unique_ptr<const Cacheable4>> cache_4;
    typedef std::map<size_t,cache_4> cache_3;
    typedef std::map<size_t,cache_3> cache_2;
    typedef std::map<size_t,cache_2> cache_t;

    mutable cache_t cache;
    mutable cache_2* i1_cache;
    mutable cache_3* i2_cache;
    mutable cache_4* i3_cache;
    mutable size_t i1,i2,i3,i4; //Current indexes
    mutable size_t itsLookups=0, itsInserts=0; //hit/miss stats
};
} // namespace qchem