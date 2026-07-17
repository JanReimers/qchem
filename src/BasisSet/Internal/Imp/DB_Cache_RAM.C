// File: BasisSet/Internal/Imp/DB_Cache_RAM.C
module;
#include <map>
#include <cassert>
#include <vector>
#include <string>
#include <variant>
#include <iostream>
#include <fstream>
#include <format>
#include <chrono>
#include <functional>
#include <stdexcept>
#include <tuple>
module qchem.BasisSet.Internal.DB_Cache_RAM;

namespace std
{
    using I1C=qchem::BasisSet::IntegralsCache_Base::I1C;
    template <> struct formatter<I1C> : formatter<string_view> {  
    auto format(I1C c, std::format_context& ctx) const {  
        string_view name;  
        switch (c) { // Reuse switch-case logic, but integrate with format  
        case I1C::Charge:   name = "Charge"; break;  
        case I1C::Normalization: name = "Normalization"; break;  
        }  
        return formatter<string_view>::format(name, ctx);  
    }  
    }; 

    using I2C=qchem::BasisSet::IntegralsCache_Base::I2C;
    template <> struct formatter<I2C> : formatter<string_view> 
    {  
        auto format(I2C c, format_context& ctx) const 
        {  
            string_view name;  
            switch (c) { // Reuse switch-case logic, but integrate with format  
                case I2C::Overlap:      name = "Overlap"     ; break;  
                case I2C::Repulsion:    name = "Repulsion"   ; break;  
                case I2C::Kinetic:      name = "Kinetic"     ; break;  
                case I2C::Nuclear:      name = "Nuclear"     ; break;  
                case I2C::RestMass:     name = "RestMass"    ; break;  
                case I2C::InvOverlap:   name = "InvOverlap"  ; break;  
                case I2C::InvRepulsion: name = "InvRepulsion"; break;  
            }  
            return formatter<string_view>::format(name, ctx);  
        }
    };  
    using I2x=qchem::BasisSet::IntegralsCache_Base::I2x;
    template <> struct formatter<I2x> : formatter<string_view> 
    {  
        auto format(I2x c, format_context& ctx) const 
        {  
            string_view name;  
            switch (c) { // Reuse switch-case logic, but integrate with format  
                case I2x::Overlap:      name = "Overlap"     ; break;  
                case I2x::Repulsion:    name = "Repulsion"   ; break;  
                case I2x::Kinetic:      name = "Kinetic"     ; break;  
            }  
            return formatter<string_view>::format(name, ctx);  
        }
    };  

    using I2n=qchem::BasisSet::IntegralsCache_Base::I2n;
    template <> struct formatter<I2n> : formatter<string_view> 
    {  
        auto format(I2n c, format_context& ctx) const 
        {  
            string_view name;  
            switch (c) { // Reuse switch-case logic, but integrate with format
                case I2n::Nuclear:      name = "Nuclear"     ; break;
                case I2n::LocalPP:      name = "LocalPP"     ; break;
                case I2n::SeparablePP:  name = "SeparablePP" ; break;
            }
            return formatter<string_view>::format(name, ctx);  
        }
    };  

    using I3C=qchem::BasisSet::IntegralsCache_Base::I3C;
    template <> struct formatter<I3C> : formatter<string_view> 
    {  
        auto format(I3C c, format_context& ctx) const 
        {  
            string_view name;  
            switch (c) { // Reuse switch-case logic, but integrate with format  
                case I3C::Overlap:   name = "Overlap"     ; break;  
                case I3C::Repulsion: name = "Repulsion"     ; break;  
            }  
            return formatter<string_view>::format(name, ctx);  
        }
    };  

    using I4C=qchem::BasisSet::IntegralsCache_Base::I4C;
    template <> struct formatter<I4C> : formatter<string_view> 
    {  
        auto format(I4C c, format_context& ctx) const 
        {  
            string_view name;  
            switch (c) { // Reuse switch-case logic, but integrate with format  
                case I4C::Direct:      name = "Direct"     ; break;  
                case I4C::Exchange:      name = "Exchange"     ; break;  
            }  
            return formatter<string_view>::format(name, ctx);  
        }
    };  

}
namespace qchem::BasisSet
{

//
// Specialize std::formatter for I1C,I2C
//
 


template <class T> IntegralsCache_RAM<T>::IntegralsCache_RAM(bool makelog)
    : itsMakeLog(makelog)
    , itsMaxRAM(4096) // in MB
    , itsTotalRAM(0) // in MB
    {
        if  (itsMakeLog)
            itsLogger=std::ofstream("cache.log");
    };

template <class T>  IntegralsCache_RAM<T>::~IntegralsCache_RAM()
{
    ReportRAMUsage(std::cout);
}
template <class T> void IntegralsCache_RAM<T>::Register(Cache4_Client* eval)
{
    RadialTypeID_t key=eval->RadialType(); //This is really a type
    auto it=itsCache4s.find(key);
    if (it==itsCache4s.end())
    {
        const auto [iterator, success]=itsCache4s.insert({key,val_t(eval->MakeCache4())});
        assert(success);
        it=iterator;
        std::cout << "Insert Cache4_Client RadialType=" << key << std::endl;
    }
    it->second->Register(eval);
}

template <class T> const Cache4* IntegralsCache_RAM<T>::GetCache4(const RadialTypeID_t& type) const
{
    auto it=itsCache4s.find(type);
    if (it==itsCache4s.end())
    {
        std::cerr << "Cache4 error: Cannot find radial type '" << type << "', known types are:" << std::endl;
        for (auto& i:itsCache4s)
            std::cerr << i.first << std::endl;
    }
    assert(it!=itsCache4s.end());
    return it->second.get();
}

template <class T> void IntegralsCache_RAM<T>::Register(Cache2_Client* eval)
{
    RadialTypeID_t key=eval->RadialType(); //This is really a type
    auto it=itsCache2s.find(key);
    if (it==itsCache2s.end())
    {
        const auto [iterator, success]=itsCache2s.insert({key,val2_t(eval->MakeCache2())});
        assert(success);
        it=iterator;
    }
    it->second->Register(eval);
}

template <class T> const Cache2* IntegralsCache_RAM<T>::GetCache2(const RadialTypeID_t& type) const
{
    auto it=itsCache2s.find(type);
    if (it==itsCache2s.end())
    {
        std::cerr << "Cache2 error: Cannot find radial type '" << type << "', known types are:" << std::endl;
        for (auto& i:itsCache2s)
            std::cerr << i.first << std::endl;
    }
    assert(it!=itsCache2s.end());
    return it->second.get();
}

template <class T> void IntegralsCache_RAM<T>::Register(Cache3_Client* eval)
{
    RadialTypeID_t key=eval->RadialType();
    auto it=itsCache3s.find(key);
    if (it==itsCache3s.end())
    {
        const auto [iterator, success]=itsCache3s.insert({key,val3_t(eval->MakeCache3())});
        assert(success);
        it=iterator;
    }
    it->second->Register(eval);
}

template <class T> const Cache3* IntegralsCache_RAM<T>::GetCache3(const RadialTypeID_t& type) const
{
    auto it=itsCache3s.find(type);
    if (it==itsCache3s.end())
    {
        std::cerr << "Cache3 error: Cannot find radial type '" << type << "', known types are:" << std::endl;
        for (auto& i:itsCache3s)
            std::cerr << i.first << std::endl;
    }
    assert(it!=itsCache3s.end());
    return it->second.get();
}

// Self-check guard: on a cache hit or insert the stored matrix MUST have the dimension the client
// expects (CacheDim() = number of basis functions).  A mismatch means the BasisSetID() key was not
// specific enough -- two differently-sized geometries collided on one key -- and we dump operator / ID /
// dims and throw right here at the cache boundary, instead of letting a wrong-sized matrix detonate
// (Cholesky segfault) several layers down.  (Same-size-different-content collisions still slip through;
// that is a rarer, separate problem.)
namespace
{
void CheckCacheDim(size_t cached,size_t expected,const std::string& op,const std::string& id)
{
    if (cached!=expected)
        throw std::runtime_error(std::format(
            "IntegralsCache dim mismatch: operator={} BasisSetID='{}' cached dim={} expected dim={}"
            " -- BasisSetID() is not specific enough (differently-sized geometries collide on one key).",
            op,id,cached,expected));
}
void CheckCacheDim2(size_t crows,size_t ccols,size_t erows,size_t ecols,
                    const std::string& op,const std::string& ida,const std::string& idb)
{
    if (crows!=erows || ccols!=ecols)
        throw std::runtime_error(std::format(
            "IntegralsCache cross dim mismatch: operator={} a='{}' b='{}' cached={}x{} expected={}x{}"
            " -- a BasisSetID() is not specific enough.",
            op,ida,idb,crows,ccols,erows,ecols));
}
} //anon

//
// Self-contained lookup-or-compute.  No shared itsLastKeyX / iterator state, so these are
// re-entrant: make() may perform nested cached Get()s.  std::map is node-based, so the
// returned reference stays valid as later (or nested) inserts grow the map.  We always
// compute make() into a local *before* inserting, holding no find-iterator across it.
//
template <class T> const rvec_t& IntegralsCache_RAM<T>::Get(I1C i1c,const DBCacheClient* bs,std::function<rvec_t()> make)
{
    IBS_ID_t id=bs->BasisSetID();
    key1_t key(i1c,id);
    if (auto i=itsVecs.find(key); i!=itsVecs.end()) return i->second;
    if (itsMakeLog) itsLogger << "I1C " << std::format("{:<12}",i1c) << " compute a=" << id << std::endl;
    auto v=make();
    const auto [it,ok]=itsVecs.insert({key,std::move(v)});
    assert(ok);
    return it->second;
}

template <class T> const hmat_t<T>& IntegralsCache_RAM<T>::Get(I2C i2c,const DBCacheClient* bs,std::function<hmat_t<T>()> make)
{
    IBS_ID_t id=bs->BasisSetID();
    key2_t key(i2c,id);
    if (auto i=itsSMats.find(key); i!=itsSMats.end())
    {
        CheckCacheDim(i->second.rows(),bs->CacheDim(),std::format("I2C {}",i2c),id);
        return i->second;
    }
    if (itsMakeLog) itsLogger << "I2C " << std::format("{:<12}",i2c) << " compute a=" << id << std::endl;
    auto v=make();
    const auto [it,ok]=itsSMats.insert({key,std::move(v)});
    assert(ok);
    CheckCacheDim(it->second.rows(),bs->CacheDim(),std::format("I2C {}",i2c),id);
    return it->second;
}

template <class T> const hmat_t<T>& IntegralsCache_RAM<T>::Get(I2n i2n,const DBCacheClient* bs,const Structure_ID_t& st,std::function<hmat_t<T>()> make)
{
    IBS_ID_t id=bs->BasisSetID();
    keyn_t key(i2n,id,st);
    if (auto i=itsNMats.find(key); i!=itsNMats.end())
    {
        CheckCacheDim(i->second.rows(),bs->CacheDim(),std::format("I2n {}",i2n),id);
        return i->second;
    }
    if (itsMakeLog) itsLogger << "I2n " << std::format("{:<12}",i2n) << " compute a=" << id << " structure=" << st << std::endl;
    auto v=make();
    const auto [it,ok]=itsNMats.insert({key,std::move(v)});
    assert(ok);
    CheckCacheDim(it->second.rows(),bs->CacheDim(),std::format("I2n {}",i2n),id);
    return it->second;
}

template <class T> const mat_t<T>& IntegralsCache_RAM<T>::Get(I2x i2x,const DBCacheClient* a,const DBCacheClient* b,std::function<mat_t<T>()> make)
{
    IBS_ID_t ida=a->BasisSetID(), idb=b->BasisSetID();
    keyx_t key(i2x,ida,idb);
    if (auto i=itsMats.find(key); i!=itsMats.end())
    {
        CheckCacheDim2(i->second.rows(),i->second.columns(),a->CacheDim(),b->CacheDim(),
                       std::format("I2x {}",i2x),ida,idb);
        return i->second;
    }
    if (itsMakeLog) itsLogger << "I2x " << std::format("{:<12}",i2x) << " compute a=" << ida << " b=" << idb << std::endl;
    auto v=make();
    const auto [it,ok]=itsMats.insert({key,std::move(v)});
    assert(ok);
    CheckCacheDim2(it->second.rows(),it->second.columns(),a->CacheDim(),b->CacheDim(),
                   std::format("I2x {}",i2x),ida,idb);
    return it->second;
}

template <class T> const ERI3<T>& IntegralsCache_RAM<T>::Get(I3C i3c,const DBCacheClient* a,const DBCacheClient* b,std::function<ERI3<T>()> make)
{
    IBS_ID_t ida=a->BasisSetID(), idb=b->BasisSetID();
    key3_t key(i3c,ida,idb);
    if (auto i=itsERI3s.find(key); i!=itsERI3s.end()) return i->second;
    if (itsMakeLog) itsLogger << "I3C " << std::format("{:<12}",i3c) << " compute a=" << ida << " b=" << idb << std::endl;
    auto v=make();
    const auto [it,ok]=itsERI3s.insert({key,std::move(v)});
    assert(ok);
    return it->second;
}

// Reciprocal-space (plane-wave) 3-centre tensor: keyed on the same (op, orbital-ID, fit-ID) axis as ERI3, but
// stored in its own map (values are G_ERI3, not ERI3<T>).  The delta support is orbital-{G} intrinsic, so the
// fit-basis half of the key just records which fitter declared coverage.  No CacheDim cross-check: G_ERI3 is
// not a leading-dimension 2-centre matrix (like the ERI3 variant, it self-describes via its columns).
template <class T> const G_ERI3& IntegralsCache_RAM<T>::Get(I3C i3c,const DBCacheClient* a,const DBCacheClient* b,std::function<G_ERI3()> make)
{
    IBS_ID_t ida=a->BasisSetID(), idb=b->BasisSetID();
    key3_t key(i3c,ida,idb);
    if (auto i=itsG_ERI3s.find(key); i!=itsG_ERI3s.end()) return i->second;
    if (itsMakeLog) itsLogger << "I3C(G) " << std::format("{:<12}",i3c) << " compute a=" << ida << " b=" << idb << std::endl;
    auto v=make();
    const auto [it,ok]=itsG_ERI3s.insert({key,std::move(v)});
    assert(ok);
    return it->second;
}

template <class T> const ERI4& IntegralsCache_RAM<T>::Get(I4C i4c,const DBCacheClient* a,const DBCacheClient* b,std::function<ERI4()> make)
{
    IBS_ID_t ida=a->BasisSetID(), idb=b->BasisSetID();
    // Canonical-only invariant (doc/ERI4Rework.md §5.2): only the a<=b (by BasisSetID) block of each
    // unordered irrep pair may ever be built or stored.  The bra-ket partner J(b,a)=J(a,b)^T is derived on
    // the fly by ERI4::ScatterBoth from the canonical block, so requesting it here is a programming error --
    // fail loudly (rigid: a THROW, live in Release too, not a debug-only assert) rather than silently
    // re-materialize a redundant orientation and quietly double the ERI cache.  Diagonal (a==b) is canonical.
    if (ida>idb)
        throw std::runtime_error(std::format(
            "ERI4 cache: non-canonical {} request a='{}' b='{}' -- the bra-ket partner J(b,a)=J(a,b)^T must "
            "never be built/stored; request the canonical (a<=b by BasisSetID) block and derive the transpose "
            "via ERI4::ScatterBoth (doc/ERI4Rework.md §5.2).", i4c, ida, idb));
    map4_t& m = (i4c==IntegralsCache<T>::I4C::Direct) ? Jac : Kab;
    if (auto ia=m.find(ida); ia!=m.end())
    {
        if (auto ib=ia->second.find(idb); ib!=ia->second.end())
        {
            ERI4_timestamps[std::make_pair(ida,idb)]=std::chrono::system_clock::now();
            return ib->second;
        }
    }
    if (itsMakeLog) itsLogger << "I4C " << std::format("{:<12}",i4c) << " compute a=" << ida << " b=" << idb << std::endl;
    auto v=make();
    const auto [it,ok]=m[ida].insert({idb,std::move(v)});
    assert(ok);
    const id_pair_t key=std::make_pair(ida,idb);
    ERI4_timestamps[key]=std::chrono::system_clock::now();
    itsTotalRAM+=it->second.size()*sizeof(T)/1024/1024;
    if (itsTotalRAM>itsMaxRAM) RunGarbageCollector(key); // protect the entry just inserted
    return it->second;
}

template <class T> const rvec_t& IntegralsCache_RAM<T>::Get(I1C i1c,const DBCacheClient* bs,const Mesh_ID_t& mid,std::function<rvec_t()> make)
{
    IBS_ID_t id=bs->BasisSetID();
    key1m_t key(i1c,id,mid);
    if (auto i=itsmVecs.find(key); i!=itsmVecs.end()) return i->second;
    if (itsMakeLog) itsLogger << "1Cm " << std::format("{:<12}",i1c) << " compute a=" << id << " mesh=" << mid << std::endl;
    auto v=make();
    const auto [it,ok]=itsmVecs.insert({key,std::move(v)});
    assert(ok);
    return it->second;
}

template <class T> const rmat_t& IntegralsCache_RAM<T>::Get(I2x i2x,const DBCacheClient* a,const DBCacheClient* b,const Mesh_ID_t& mid,std::function<rmat_t()> make)
{
    IBS_ID_t ida=a->BasisSetID(), idb=b->BasisSetID();
    key2xm_t key(i2x,ida,idb,mid);
    if (auto i=itsmMats.find(key); i!=itsmMats.end()) return i->second;
    if (itsMakeLog) itsLogger << "2xm " << std::format("{:<12}",i2x) << " compute a=" << ida << " b=" << idb << " mesh=" << mid << std::endl;
    auto v=make();
    const auto [it,ok]=itsmMats.insert({key,std::move(v)});
    assert(ok);
    return it->second;
}

// --- Troubleshooting / DBCache-unit-test Clear() hooks (see the declarations in the interface partition).
// kCacheTestHooks is the compile-time switch -- the constexpr-bool idiom, NO preprocessor.  Default true so the DBCache unit tests run in the normal Release UTMain; flip to false
// to drop the bodies from a hardened production build (they are unreachable from production code anyway,
// which only holds the abstract IntegralsCache<T> face).  Each erases the entries for one operator; the
// enum value comparison is on std::get<0>(key) (the operator slot) for every keyed map.
constexpr bool kCacheTestHooks = true;

template <class T> void IntegralsCache_RAM<T>::Clear(I1C op)
{
    if constexpr (kCacheTestHooks)
    {
        std::erase_if(itsVecs,  [&](const auto& kv){ return std::get<0>(kv.first)==op; });
        std::erase_if(itsmVecs, [&](const auto& kv){ return std::get<0>(kv.first)==op; });
    }
}
template <class T> void IntegralsCache_RAM<T>::Clear(I2C op)
{
    if constexpr (kCacheTestHooks)
        std::erase_if(itsSMats, [&](const auto& kv){ return std::get<0>(kv.first)==op; });
}
template <class T> void IntegralsCache_RAM<T>::Clear(I2n op)
{
    // keyn_t = (I2n, IBS_ID_t, Structure_ID_t): erase only the requested operator's entries (Nuclear vs the
    // LocalPP/SeparablePP pseudopotential matrices), mirroring the per-operator Clear(I2C)/Clear(I4C) hooks.
    if constexpr (kCacheTestHooks)
        std::erase_if(itsNMats, [&](const auto& kv){ return std::get<0>(kv.first)==op; });
}
template <class T> void IntegralsCache_RAM<T>::Clear(I2x op)
{
    if constexpr (kCacheTestHooks)
    {
        std::erase_if(itsMats,  [&](const auto& kv){ return std::get<0>(kv.first)==op; });
        std::erase_if(itsmMats, [&](const auto& kv){ return std::get<0>(kv.first)==op; });
    }
}
template <class T> void IntegralsCache_RAM<T>::Clear(I3C op)
{
    if constexpr (kCacheTestHooks)
    {
        std::erase_if(itsERI3s,   [&](const auto& kv){ return std::get<0>(kv.first)==op; });
        std::erase_if(itsG_ERI3s, [&](const auto& kv){ return std::get<0>(kv.first)==op; });
    }
}
template <class T> void IntegralsCache_RAM<T>::Clear(I4C op)
{
    // Jac/Kab are nested id->id->ERI4 maps with no operator in the key; Direct lives in Jac, Exchange in
    // Kab.  ERI4_timestamps / itsTotalRAM are advisory GC bookkeeping (stale entries are tolerated and
    // self-correct), so a troubleshooting clear leaves them be.
    if constexpr (kCacheTestHooks)
        (op==I4C::Direct ? Jac : Kab).clear();
}

template struct IntegralsCache<double>;
template struct IntegralsCache_RAM<double>;
// Complex (plane-wave) path: the PW lineage only needs the Hermitian I2C Overlap, but the explicit
// instantiation pulls in every member, so all the dcmplx storage/Get variants must compile.
template struct IntegralsCache<dcmplx>;
template struct IntegralsCache_RAM<dcmplx>;

} //namespace
