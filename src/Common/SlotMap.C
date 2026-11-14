// File: Common/SlotMap.C  THE keyed cache of the code: find-or-create under a lock, a once-fill per slot, two disciplines.
//
// doc/IntegralCacheOwnershipPlan.md S1 (user rulings 2026-11-12, §8).  Every owner-held cache (a Hamiltonian term's blocks now; the
// basis-owned integral slots from S2 on) is one of these, so the locking and the fill rule are written ONCE.
//
// TWO DISCIPLINES, chosen by the type, so a misuse does not compile (CLAUDE.md "lazy caching only for what never changes"):
//   - \c SlotDiscipline::FillOnce: a pure function of fixed inputs (an integral, a static term's block).  \c GetOrFill builds on the first
//     ask and serves the same value forever.
//   - \c SlotDiscipline::Refresh: anything that depends on the density.  There is NO lazy fill: the only way to fill is through the
//     \c RefreshScope that \c BeginRefresh opens (it clears the map first), which is the "everybody, start using this density" line of
//     the SCF loop (V1.41).  Outside it a query is a pure \c Read.
//
// THREADING (the plan §3): the map lock is held only to find or create a node, never across a fill, so a fill that itself asks another
// slot of the same map does not deadlock, and fills of different slots run concurrently.  Each slot carries a \c std::once_flag, so two
// threads asking the SAME slot build it once (the other waits).  A fill that throws leaves its slot unfilled (the next ask retries).
// std::map nodes are stable, so a reference handed out stays valid until the map is cleared (a Refresh) or destroyed.
module;
#include <map>
#include <mutex>
#include <memory>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <cstdint>
#include <vector>
export module qchem.SlotMap;
import qchem.InstanceToken;

export namespace qchem
{

enum class SlotDiscipline
{
    FillOnce,   //!< fixed inputs: built on first ask, never cleared
    Refresh     //!< density-dependent: rebuilt only inside a RefreshScope, read purely otherwise
};

//! \brief A keyed cache of \a V: find-or-create under a lock, then a per-slot once-fill (see the file comment for the two disciplines).
//! Copyable: a copy holds the FILLED values (fresh locks and once-flags), so an owner that is copied keeps its built blocks.
template <class Key, class V, SlotDiscipline D> class SlotMap
{
    struct Slot
    {
        std::once_flag once;
        std::unique_ptr<V> value;    //!< null until filled (a value type need not have an "empty" sentinel)
    };
public:
    SlotMap() = default;
    SlotMap(const SlotMap& o) {CopyFilled(o);}
    SlotMap& operator=(const SlotMap& o)
    {
        if (this!=&o) {std::lock_guard<std::mutex> l(itsLock); itsSlots.clear(); CopyFilledLocked(o);}
        return *this;
    }

    //! FillOnce only: the slot for \a key, built by \a make() on the first ask (once, whichever thread asks first).
    template <class F> const V& GetOrFill(const Key& key, F&& make) const requires (D==SlotDiscipline::FillOnce)
    {
        return FillSlot(Node(key), std::forward<F>(make));
    }

    //! A PURE READ: the filled value for \a key, or null if there is none (absent, or not filled yet).  Never builds.
    const V* Read(const Key& key) const
    {
        std::lock_guard<std::mutex> l(itsLock);
        auto it=itsSlots.find(key);
        return it==itsSlots.end() ? nullptr : it->second.value.get();
    }

    //! \brief Refresh only: the one place a density-dependent slot is filled.  Opening it CLEARS the map; while it is open,
    //! \c Fill builds a slot that is not yet filled (so two keys a term folds together are built once).  One refresh at a time.
    class RefreshScope
    {
    public:
        ~RefreshScope() {std::lock_guard<std::mutex> l(itsMap.itsLock); itsMap.itsRefreshing=false;}
        RefreshScope(const RefreshScope&)=delete;
        RefreshScope& operator=(const RefreshScope&)=delete;
        template <class F> const V& Fill(const Key& key, F&& make) const
        {
            return itsMap.FillSlot(itsMap.Node(key), std::forward<F>(make));   // FillSlot is a member: it takes the map lock
        }
    private:
        friend class SlotMap;
        explicit RefreshScope(const SlotMap& m) : itsMap(m) {}
        const SlotMap& itsMap;
    };
    [[nodiscard]] RefreshScope BeginRefresh() const requires (D==SlotDiscipline::Refresh)
    {
        std::lock_guard<std::mutex> l(itsLock);
        if (itsRefreshing) throw std::logic_error("SlotMap: a refresh opened inside another refresh of the same map");
        itsRefreshing=true;
        itsSlots.clear();
        return RefreshScope(*this);
    }

    //! The bytes held by the filled slots, each measured by \a bytesOf(const V&) (the report's number, plan S6).
    template <class F> size_t Bytes(F&& bytesOf) const
    {
        std::lock_guard<std::mutex> l(itsLock);
        size_t n=0;
        for (const auto& [k,s] : itsSlots) if (s.value) n+=bytesOf(*s.value);
        return n;
    }
    //! Erase every FILLED slot whose key satisfies \a dead(key), and return how many went.  For an owner that KNOWS such a key can never
    //! be asked again (an entry keyed by a dead partner's never-reused serial, doc/IntegralCacheOwnershipPlan.md §3b): nobody can be
    //! building it or holding a reference to it, because any asker would be holding the partner alive.  An unfilled slot is left alone.
    template <class P> size_t EraseIf(P&& dead) const requires (D==SlotDiscipline::FillOnce)
    {
        std::lock_guard<std::mutex> l(itsLock);
        size_t n=0;
        for (auto it=itsSlots.begin(); it!=itsSlots.end();)
            if (it->second.value && dead(it->first)) {it=itsSlots.erase(it); n++;}
            else ++it;
        return n;
    }
    //! The number of FILLED slots.
    size_t Filled() const
    {
        std::lock_guard<std::mutex> l(itsLock);
        size_t n=0;
        for (const auto& [k,s] : itsSlots) if (s.value) n++;
        return n;
    }

private:
    Slot& Node(const Key& key) const
    {
        std::lock_guard<std::mutex> l(itsLock);
        return itsSlots.try_emplace(key).first->second;
    }
    //! The build runs OUTSIDE the map lock (a nested ask of another slot must not deadlock); only the pointer store takes it,
    //! so a concurrent \c Read of this slot sees null or the finished value, never a half-built one.
    template <class F> const V& FillSlot(Slot& s, F&& make) const
    {
        std::call_once(s.once, [&]
        {
            auto built=std::make_unique<V>(make());
            std::lock_guard<std::mutex> l(itsLock);
            s.value=std::move(built);
        });
        return *s.value;
    }
    void CopyFilled(const SlotMap& o) {std::lock_guard<std::mutex> l(itsLock); CopyFilledLocked(o);}
    void CopyFilledLocked(const SlotMap& o)
    {
        std::lock_guard<std::mutex> lo(o.itsLock);
        for (const auto& [k,s] : o.itsSlots)
            if (s.value)
            {
                Slot& mine=itsSlots.try_emplace(k).first->second;
                std::call_once(mine.once, [&]{mine.value=std::make_unique<V>(*s.value);});
            }
    }

    mutable std::mutex           itsLock;
    mutable std::map<Key,Slot>   itsSlots;
    mutable bool                 itsRefreshing=false;
};

//! \brief ONE cached value with the same promises as a \c FillOnce \c SlotMap slot (built once on the first ask, whichever thread
//! asks; a throwing build leaves it unfilled; a copy keeps a filled value).  For an owner's single integrals: a fit basis's metric,
//! its inverse, its charges (doc/IntegralCacheOwnershipPlan.md S2).
template <class V> class OnceSlot
{
public:
    template <class F> const V& GetOrFill(F&& make) const {return itsSlot.GetOrFill(0, std::forward<F>(make));}
    const V* Read() const {return itsSlot.Read(0);}
    template <class F> size_t Bytes(F&& bytesOf) const {return itsSlot.Bytes(std::forward<F>(bytesOf));}
private:
    SlotMap<int,V,SlotDiscipline::FillOnce> itsSlot;
};

//! \brief A \c FillOnce slot per PARTNER OBJECT, keyed by the partner's \c InstanceToken serial, that frees a dead partner's slot.
//! For an owner caching something over ANOTHER object it does not own (a δ fit's \f$\Phi_{gi}=\chi_i(r_g)\f$ table per orbital block:
//! doc/IntegralCacheOwnershipPlan.md S3c, user ruling 2026-11-14).  Never keyed by the partner's address (reused after it dies) or by a
//! content stand-in such as its symmetry (two partners with equal symmetry are still two partners: the latent defect S3c removed).
//! Every \c GetOrFill first sweeps out the slots of partners that no longer exist (\c EvictDead), as \c BasisCache does for ⟨ab|c⟩:
//! a dead partner's serial is never reused, so nobody can be building or holding its slot.
template <class V> class InstanceSlotMap
{
public:
    template <class F> const V& GetOrFill(const InstanceToken& partner, F&& make) const
    {
        EvictDead();
        {
            std::lock_guard<std::mutex> l(itsLock);
            itsAlive.try_emplace(partner.Serial(), partner.Alive());
        }
        return itsSlots.GetOrFill(partner.Serial(), std::forward<F>(make));
    }
    //! Free the slots of every partner that no longer exists; returns how many slots went.
    size_t EvictDead() const
    {
        std::vector<uint64_t> dead;
        {
            std::lock_guard<std::mutex> l(itsLock);
            for (auto it=itsAlive.begin(); it!=itsAlive.end();)
                if (it->second.expired()) {dead.push_back(it->first); it=itsAlive.erase(it);}
                else ++it;
        }
        if (dead.empty()) return 0;
        return itsSlots.EraseIf([&](uint64_t k){for (uint64_t d:dead) if (d==k) return true; return false;});
    }
    size_t Filled() const {return itsSlots.Filled();}
    template <class F> size_t Bytes(F&& bytesOf) const {return itsSlots.Bytes(std::forward<F>(bytesOf));}

    InstanceSlotMap() = default;
    InstanceSlotMap(const InstanceSlotMap& o) : itsSlots(o.itsSlots)
        {std::lock_guard<std::mutex> l(o.itsLock); itsAlive=o.itsAlive;}
    InstanceSlotMap& operator=(const InstanceSlotMap&) = delete;

private:
    SlotMap<uint64_t,V,SlotDiscipline::FillOnce>              itsSlots;
    mutable std::mutex                                       itsLock;    //!< itsAlive only, never held across a fill
    mutable std::map<uint64_t,std::weak_ptr<const uint64_t>> itsAlive;
};

} // namespace qchem
