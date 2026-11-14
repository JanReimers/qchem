// File: BasisSet/Internal/PrimitiveStore.C  The primitive-level caches of ONE owner: charge distributions Ω, RNLM tables, Hermite3 blocks.
//
// doc/IntegralCacheOwnershipPlan.md S5 (user rulings 1 and 4, 2026-11-12; controller rulings 2026-11-14 on the two open choices):
//   - It belongs to an OWNER and dies with it.  A basis set's BasisCache holds one for its blocks (the within-geometry reuse: Ω across
//     irreps and quartets of one point, §4), and a fill whose keys reach outside its owner (a lattice image term, a block × fit 3-centre
//     fill on the per-component kernels) makes a scratch one on the stack and drops it at scope end.  That replaces the old process-wide
//     Cache2 "PG.Omega" / "PG.RNLM" / Cache3 "PG.H3" and the size-0 GeometryCacheBudget that kept them O(1) on the lattice paths.
//   - So the parent's store only ever holds keys over the parent's OWN primitives: no entry can outlive a partner, and nothing needs a
//     sweep.  A key is a primitive's (or a pair's) never-reused SERIAL (PrimGaussian / Ω take one from a process-wide atomic at
//     construction), so two owners' keys can never collide even when their primitives are equal.  (No content interning: measured
//     2026-11-14, a molecular water run constructs 43 primitives, all distinct, and a GPW Si run's 122,828 constructions are lattice image
//     clones whose entries the old size-0 budget evicted anyway: interning bought nothing.)
//   - It is reached EXPLICITLY: every M&D kernel that looks something up takes the store as a parameter (controller ruling (A)), so which
//     owner a lookup lands in is visible at the call, never a hidden rule.
//
// THREADING (plan §3, option c): single-thread by contract.  Only a fill's serial prepass touches a store (the threaded group kernels
// receive the Ω they need as a PrimitivePairTable built before the parallel region).  In an asserts-live build every access checks that
// no other thread is inside the same store.
module;
#include <atomic>
#include <cassert>
#include <thread>
export module qchem.BasisSet.Internal.PrimitiveStore;
export import qchem.BasisSet.Internal.Cache2;
export import qchem.BasisSet.Internal.Cache3;

export namespace qchem::BasisSet
{

class PrimitiveStore
{
public:
    PrimitiveStore() = default;
    PrimitiveStore(const PrimitiveStore&)            = delete;   // owned by one owner; a copy would be a second owner of the same keys
    PrimitiveStore& operator=(const PrimitiveStore&) = delete;

    //! RAII proof that the caller is the only thread inside this store (checked only where asserts are live).
    class Use
    {
    public:
        explicit Use(const PrimitiveStore& s) : itsStore(s)
        {
#ifndef NDEBUG
            const std::thread::id me=std::this_thread::get_id();
            std::thread::id none;
            if (!itsStore.itsUser.compare_exchange_strong(none, me)) assert(none==me && "PrimitiveStore used by two threads at once");
            else itsOwned=true;
#endif
        }
        ~Use()
        {
#ifndef NDEBUG
            if (itsOwned) itsStore.itsUser.store(std::thread::id());
#endif
        }
        Use(const Use&)=delete; Use& operator=(const Use&)=delete;
    private:
        const PrimitiveStore& itsStore;
        bool itsOwned=false;
    };

    const Cache2& Omega() const {return itsΩ;}      //!< Ω of a primitive pair, keyed (serial a, serial b)
    const Cache2& RNLM () const {return itsRNLM;}   //!< the 3/4-centre RNLM, keyed (Ω serial, primitive or Ω serial)
    const Cache3& H3   () const {return itsH3;}     //!< the 3-centre Hermite block, keyed by the primitive triple's serials

    size_t RAMsize() const {return itsΩ.RAMsize()+itsRNLM.RAMsize()+itsH3.RAMsize();}

    //! A never-reused serial for a primitive or a pair: the key component of every entry above (process-wide, so no two owners collide).
    static size_t NextSerial()
    {
        static std::atomic<size_t> theSerial{0};
        return ++theSerial;
    }

private:
    Cache2 itsΩ, itsRNLM;
    Cache3 itsH3;
    mutable std::atomic<std::thread::id> itsUser{};
};

} // namespace qchem::BasisSet
