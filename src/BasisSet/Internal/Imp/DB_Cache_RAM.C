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
module qchem.BasisSet.Internal.DB_Cache_RAM;

namespace std
{
    using I1C=BasisSet::IntegralsCache_Base::I1C;
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

    using I2C=BasisSet::IntegralsCache_Base::I2C;
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
    using I2x=BasisSet::IntegralsCache_Base::I2x;
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

    using I2n=BasisSet::IntegralsCache_Base::I2n;
    template <> struct formatter<I2n> : formatter<string_view> 
    {  
        auto format(I2n c, format_context& ctx) const 
        {  
            string_view name;  
            switch (c) { // Reuse switch-case logic, but integrate with format  
                case I2n::Nuclear:      name = "Nuclear"     ; break;  
            }  
            return formatter<string_view>::format(name, ctx);  
        }
    };  

    using I3C=BasisSet::IntegralsCache_Base::I3C;
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

    using I4C=BasisSet::IntegralsCache_Base::I4C;
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
namespace BasisSet
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

template <class T> const smat_t<T>& IntegralsCache_RAM<T>::Get(I2C i2c,const DBCacheClient* bs,std::function<smat_t<T>()> make)
{
    IBS_ID_t id=bs->BasisSetID();
    key2_t key(i2c,id);
    if (auto i=itsSMats.find(key); i!=itsSMats.end()) return i->second;
    if (itsMakeLog) itsLogger << "I2C " << std::format("{:<12}",i2c) << " compute a=" << id << std::endl;
    auto v=make();
    const auto [it,ok]=itsSMats.insert({key,std::move(v)});
    assert(ok);
    return it->second;
}

template <class T> const smat_t<T>& IntegralsCache_RAM<T>::Get(I2n i2n,const DBCacheClient* bs,const Cluster_ID_t& cl,std::function<smat_t<T>()> make)
{
    IBS_ID_t id=bs->BasisSetID();
    keyn_t key(id,cl);
    if (auto i=itsNMats.find(key); i!=itsNMats.end()) return i->second;
    if (itsMakeLog) itsLogger << "I2n " << std::format("{:<12}",i2n) << " compute a=" << id << " cluster=" << cl << std::endl;
    auto v=make();
    const auto [it,ok]=itsNMats.insert({key,std::move(v)});
    assert(ok);
    return it->second;
}

template <class T> const mat_t<T>& IntegralsCache_RAM<T>::Get(I2x i2x,const DBCacheClient* a,const DBCacheClient* b,std::function<mat_t<T>()> make)
{
    IBS_ID_t ida=a->BasisSetID(), idb=b->BasisSetID();
    keyx_t key(i2x,ida,idb);
    if (auto i=itsMats.find(key); i!=itsMats.end()) return i->second;
    if (itsMakeLog) itsLogger << "I2x " << std::format("{:<12}",i2x) << " compute a=" << ida << " b=" << idb << std::endl;
    auto v=make();
    const auto [it,ok]=itsMats.insert({key,std::move(v)});
    assert(ok);
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

template <class T> const ERI4& IntegralsCache_RAM<T>::Get(I4C i4c,const DBCacheClient* a,const DBCacheClient* b,std::function<ERI4()> make)
{
    IBS_ID_t ida=a->BasisSetID(), idb=b->BasisSetID();
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

template struct IntegralsCache<double>;
template struct IntegralsCache_RAM<double>;

} //namespace
