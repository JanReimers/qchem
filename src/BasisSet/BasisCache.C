// File: BasisSet/BasisCache.C  The cache a BasisSet MAKES and hands to every block it builds: the tables over its own functions.
//
// doc/IntegralCacheOwnershipPlan.md S3a (the user's rulings 2026-11-12, §8: one local framework, nothing process-global; the handle is a
// shared_ptr; "the cache becomes a constructor argument of the irrep basis set.  The actual BasisSet will make a cache.").
//
// What lives here, and what is left of each key once the owner is implicit (§2):
//   - a block's own 2-centre tables (Overlap, ⟨p²⟩, RestMass, the RKB L/S cross kinetic): the block's index in this cache;
//   - ERI4 Direct / Exchange between two blocks of one parent: the UNORDERED pair of indices, stored canonical (i ≤ j);
//   - ⟨ab|c⟩ Overlap3C / Repulsion3C: (the block's index, the fit basis's INSTANCE serial).  The element type follows the FIT axis, so
//     there is one typed map per scalar, the way BasisSetImp holds blocks of both scalars.
// No key carries geometry: a new geometry is a new BasisSet, so a new cache (§6).
//
// A PARTNER can die before this cache (a fit basis rebuilt mid-run).  Its serial is never reused, so its entries can never be hit again:
// every 3-centre insert sweeps them out (the weak handle beside each serial says when), and says so on the console and in the run report
// under the caller's scope.  In a normal run the fits die WITH the orbital basis, so a non-zero count is itself news.
//
// Threading: each slot is a SlotMap once-fill (a short lock to find the slot, the build outside any lock, other askers of THAT slot wait),
// so a fill may itself ask the cache for another slot.
module;
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <format>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
export module qchem.BasisSet.BasisCache;
export import qchem.InstanceToken;
export import qchem.BasisSet.PairProjection;   // PairProjection<T>, the 3-centre record
import qchem.BasisSet.Internal.ERI4;           // ERI4, the 4-centre tensor (an Internal module: imported, never re-exported)
export import qchem.Types;
import qchem.SlotMap;
import qchem.Reporting;

export namespace qchem::BasisSet
{

class BasisCache
{
public:
    //! A block's own 2-centre tables (Hermitian).
    enum class Own   {Overlap, Kinetic, RestMass};
    //! The two ERI4 contractions between two blocks of this parent.
    enum class Pair  {Direct, Exchange};
    //! ⟨ab|c⟩ against a fit basis: the overlap and the Coulomb kinds.
    enum class Three {Overlap, Repulsion};

    BasisCache() = default;
    BasisCache(const BasisCache&)            = delete;   // shared by its parent and blocks through one shared_ptr, never copied
    BasisCache& operator=(const BasisCache&) = delete;

    //! A new block's index in this cache (its residual key).  Called once, by the block's constructor.
    size_t Enrol() {return itsBlocks++;}

    template <class T, class F> const hmat_t<T>& Get(Own op, size_t block, F&& make) const
    {
        return Own2<T>().GetOrFill(std::make_tuple(op, block), std::forward<F>(make));
    }
    //! The RKB L/S cross kinetic: a table over the large block and the small-component partner it OWNS, so keyed by the block alone.
    template <class T, class F> const mat_t<T>& GetRKBCross(size_t block, F&& make) const
    {
        return Cross<T>().GetOrFill(block, std::forward<F>(make));
    }
    //! ERI4 between blocks \a a and \a b: only the canonical orientation a ≤ b is stored, so the other one is refused (a caller asking
    //! for it would hold a second copy of the same integrals).
    template <class F> const ERI4& Get(Pair op, size_t a, size_t b, F&& make) const
    {
        if (a>b)
            throw std::logic_error(std::format("BasisCache: ERI4 asked for the non-canonical block pair ({},{}); ask ({},{}) "
                                               "and scatter both (Orbital_ERI4_IBS::ScatterBoth)", a, b, b, a));
        return itsERI4.GetOrFill(std::make_tuple(op, a, b), std::forward<F>(make));
    }
    //! ⟨ab|c⟩ of block \a block against the partner \a fit (named \a fitName for the report).  Sweeps dead partners' entries first.
    template <class TFit, class F> const PairProjection<TFit>& Get(Three op, size_t block, const InstanceToken& fit, const std::string& fitName,
                                                                   F&& make) const
    {
        EvictDeadPartners();
        {
            std::lock_guard<std::mutex> l(itsPartnerLock);
            itsPartners.try_emplace(fit.Serial(), Partner{fit.Alive(), fitName});
        }
        return ThreeC<TFit>().GetOrFill(std::make_tuple(op, block, fit.Serial()), std::forward<F>(make));
    }

    //! Sweep out every entry keyed by a partner that no longer exists, report each partner and the count, and return the count.
    size_t EvictDeadPartners() const
    {
        std::map<uint64_t,std::string> dead;
        {
            std::lock_guard<std::mutex> l(itsPartnerLock);
            for (auto it=itsPartners.begin(); it!=itsPartners.end();)
                if (it->second.alive.expired()) {dead.emplace(it->first, it->second.name); it=itsPartners.erase(it);}
                else ++it;
        }
        if (dead.empty()) return 0;
        auto isDead=[&](const auto& key){return dead.count(std::get<2>(key))>0;};
        size_t n=its3C_r.EraseIf(isDead)+its3C_c.EraseIf(isDead);
        std::string who;
        for (const auto& [serial,name] : dead) who+=std::format("{}{} #{}", who.empty() ? "" : ", ", name, serial);
        std::cout << std::format("[basis cache] evicted {} 3-centre table(s) of {} dead partner(s): {}", n, dead.size(), who) << std::endl;
        report::Set("evictedPartnerTables", n);
        report::Set("evictedPartners", who);
        return n;
    }

    //! Filled entries, for the unit tests (the eviction gate counts them).
    size_t NumThreeCentre() const {return its3C_r.Filled()+its3C_c.Filled();}

private:
    template <class V, class K> using once_t = SlotMap<K,V,SlotDiscipline::FillOnce>;
    using own_key   = std::tuple<Own,size_t>;
    using pair_key  = std::tuple<Pair,size_t,size_t>;
    using three_key = std::tuple<Three,size_t,uint64_t>;

    template <class T> const once_t<hmat_t<T>,own_key>& Own2() const
    {
        if constexpr (std::is_same_v<T,double>) return itsOwn_r; else return itsOwn_c;
    }
    template <class T> const once_t<mat_t<T>,size_t>& Cross() const
    {
        if constexpr (std::is_same_v<T,double>) return itsCross_r; else return itsCross_c;
    }
    template <class T> const once_t<PairProjection<T>,three_key>& ThreeC() const
    {
        if constexpr (std::is_same_v<T,double>) return its3C_r; else return its3C_c;
    }

    struct Partner
    {
        std::weak_ptr<const uint64_t> alive;
        std::string                   name;
    };

    std::atomic<size_t>                       itsBlocks{0};
    once_t<hmat_t<double>,own_key>            itsOwn_r;
    once_t<hmat_t<dcmplx>,own_key>            itsOwn_c;
    once_t<mat_t<double>,size_t>              itsCross_r;
    once_t<mat_t<dcmplx>,size_t>              itsCross_c;
    once_t<ERI4,pair_key>                     itsERI4;
    once_t<PairProjection<double>,three_key>  its3C_r;
    once_t<PairProjection<dcmplx>,three_key>  its3C_c;
    mutable std::mutex                        itsPartnerLock;
    mutable std::map<uint64_t,Partner>        itsPartners;
};

} // namespace qchem::BasisSet
