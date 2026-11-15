// File: BasisSet/BasisCache.C  The cache a BasisSet MAKES and hands to every block it builds: the tables over its own functions.
//
// doc/IntegralCacheOwnershipPlan.md S3a (the user's rulings 2026-11-12, §8: one local framework, nothing process-global; the handle is a
// shared_ptr; "the cache becomes a constructor argument of the irrep basis set.  The actual BasisSet will make a cache.").
//
// What lives here, and what is left of each key once the owner is implicit (§2):
//   - a block's own 2-centre tables (Overlap, ⟨p²⟩, RestMass, the RKB L/S cross kinetic): the block's index in this cache;
//   - ERI4 Direct / Exchange between two blocks of one parent: the UNORDERED pair of indices, stored canonical (i ≤ j);
//   - ⟨ab|c⟩ Overlap3C / Repulsion3C: (the block's index, the fit basis's INSTANCE serial).  The element type follows the FIT axis, so
//     there is one typed map per scalar, the way BasisSetImp holds blocks of both scalars;
//   - the Gaussian primitive caches (Ω, RNLM, H3) over this parent's own primitives: S5b.  Its blocks' evaluators share it, and a fill
//     whose keys reach outside this parent uses a scratch store of its own instead, so nothing here can outlive a partner;
//   - an atom basis's Rk Slater integrals (Cache4, its exponent grouper inside), one per radial family: S5a.  Its blocks (one per l)
//     register into it as they are built, so their exponent indices are THIS parent's, and an atom run's tables die with its basis.
// No key carries geometry: a new geometry is a new BasisSet, so a new cache (§6).
//
// A PARTNER can die before this cache (a fit basis rebuilt mid-run).  Its serial is never reused, so its entries can never be hit again:
// every 3-centre insert sweeps them out (the weak handle beside each serial says when), and says so on the console and in the run report
// under the caller's scope.  In a normal run the fits die WITH the orbital basis, so a non-zero count is itself news.
//
// THE REPORT (S6; the run-report rule, user 2026-11-02: the callee writes one level in, AS IT RUNS): every fill adds its own bytes and a count
// under `cache` at the caller's cursor -- `<kind>Bytes` / `<kind>Entries`, plus `primitiveBytes` / `rkBytes` for what the fill's integrals put
// in this parent's primitive store and atom Rk tables -- so the report shows where each table was committed, in run sequence, and nobody
// collects anything afterwards.
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
import qchem.BasisSet.Internal.Cache4;         // Cache4 / Cache4_Client: an atom basis's Rk tables (imported, never re-exported)
import qchem.BasisSet.Internal.PrimitiveStore; // the Gaussian primitive caches of this parent (S5; imported, never re-exported)
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
        return Own2<T>().GetOrFill(std::make_tuple(op, block), Reported(NameOf(op), std::forward<F>(make)));
    }
    //! The RKB L/S cross kinetic: a table over the large block and the small-component partner it OWNS, so keyed by the block alone.
    template <class T, class F> const mat_t<T>& GetRKBCross(size_t block, F&& make) const
    {
        return Cross<T>().GetOrFill(block, Reported("rkbCross", std::forward<F>(make)));
    }
    //! ERI4 between blocks \a a and \a b: only the canonical orientation a ≤ b is stored, so the other one is refused (a caller asking
    //! for it would hold a second copy of the same integrals).
    template <class F> const ERI4& Get(Pair op, size_t a, size_t b, F&& make) const
    {
        if (a>b)
            throw std::logic_error(std::format("BasisCache: ERI4 asked for the non-canonical block pair ({},{}); ask ({},{}) "
                                               "and scatter both (Orbital_ERI4_IBS::ScatterBoth)", a, b, b, a));
        return itsERI4.GetOrFill(std::make_tuple(op, a, b), Reported(op==Pair::Direct ? "eri4Direct" : "eri4Exchange", std::forward<F>(make)));
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
        return ThreeC<TFit>().GetOrFill(std::make_tuple(op, block, fit.Serial()),
                                        Reported(op==Three::Overlap ? "overlap3C" : "repulsion3C", std::forward<F>(make)));
    }

    //! S5b: the primitive caches of this parent's blocks (their evaluators co-own it, so the handle is shared).
    std::shared_ptr<const PrimitiveStore> Primitives() const {return itsPrimitives;}

    //! S5a: block \a eval of an atom basis joins its radial family's Rk table (made on the first block of that family).  Called by the
    //! block's constructor, which runs serially inside the parent's factory.
    void Register4(Cache4_Client* eval) const
    {
        Cache4* c=nullptr;
        {
            std::lock_guard<std::mutex> l(its4Lock);
            auto it=its4.find(eval->RadialType());
            if (it==its4.end())
            {
                it=its4.emplace(eval->RadialType(), std::unique_ptr<Cache4>(eval->MakeCache4())).first;
                std::cout << "Insert Cache4_Client RadialType=" << it->first << std::endl;
            }
            c=it->second.get();
        }
        c->Register(eval);
    }
    //! The Rk table of radial family \a radialType (the ERI4 builders' lookup).  THROWS if no block of that family registered here: an
    //! ERI4 pair is two blocks of ONE parent, so a missing family is a wiring error.
    const Cache4& GetCache4(const std::string& radialType) const
    {
        std::lock_guard<std::mutex> l(its4Lock);
        auto it=its4.find(radialType);
        if (it==its4.end())
            throw std::logic_error("BasisCache::GetCache4: no block of radial family '"+radialType+"' registered with this basis set");
        return *it->second;
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
    static const char* NameOf(Own op) {return op==Own::Overlap ? "overlap" : op==Own::Kinetic ? "kinetic" : "restMass";}
    // The bytes a filled value holds (the doubles or complex numbers it stores; blaze's scalar symmetric/Hermitian adaptors keep the full
    // square).  A matrix-free ⟨ab|c⟩ realization holds closures whose captures this cannot see: it reports its dense part only.
    template <class M> static double MatrixBytes(const M& m) {return double(m.rows())*double(m.columns())*sizeof(typename M::ElementType);}
    template <class M> requires requires (const M& m) {m.rows(); m.columns();}
    static double BytesOf(const M& m) {return MatrixBytes(m);}   // every dense table (hmat_t / mat_t are aliases: T is not deducible through them)
    static double BytesOf(const ERI4& e)   // one value per unordered outer pair (blaze's non-scalar adaptor shares (a,b) and (b,a)), each full
    {
        const double nab=double(e.Nab()), ncd=e.Nab() ? double(e.Ncd()) : 0.0;
        return 0.5*nab*(nab+1.0)*ncd*ncd*sizeof(double);
    }
    template <class T> static double BytesOf(const PairProjection<T>& p)
    {
        double b=0.0;
        for (const auto& d : p.dense) b+=MatrixBytes(d);
        return b;
    }
    double PrimitiveBytes() const {return double(itsPrimitives->RAMsize());}
    double RkBytes() const
    {
        std::lock_guard<std::mutex> l(its4Lock);
        double b=0.0;
        for (const auto& [k,c] : its4) b+=double(c->RAMsize());
        return b;
    }
    //! \a make wrapped so the fill reports itself as it runs (see the file comment).  The primitive store and the Rk tables are single-thread
    //! by contract (a fill's serial prepass), so their growth across one fill is that fill's own.
    template <class F> auto Reported(const char* kind, F&& make) const
    {
        return [this, kind, make=std::forward<F>(make)]() mutable
        {
            const double p0=PrimitiveBytes(), r0=RkBytes();
            auto v=make();
            report::Accumulate("cache", std::string(kind)+"Bytes", BytesOf(v));
            report::Accumulate("cache", std::string(kind)+"Entries", 1.0);
            if (const double dp=PrimitiveBytes()-p0; dp!=0.0) report::Accumulate("cache", "primitiveBytes", dp);
            if (const double dr=RkBytes()-r0;        dr!=0.0) report::Accumulate("cache", "rkBytes", dr);
            return v;
        };
    }

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
    std::shared_ptr<PrimitiveStore>                        itsPrimitives=std::make_shared<PrimitiveStore>();
    mutable std::mutex                                     its4Lock;   //!< its4 only, never held across a Cache4 fill
    mutable std::map<std::string,std::unique_ptr<Cache4>>  its4;       //!< the atom Rk tables, by radial family (S5a)
};

} // namespace qchem::BasisSet
