// File: BasisSet/Imp/DB_Cache1.C Global integrals cache allow data sharing between separate runs.
module;
#include <map>
#include <cassert>
#include <vector>
#include <string>
#include <variant>
#include <iostream>
#include <fstream>
#include <format>
module qchem.BasisSet1.DB_Cache;

// import qchem.BasisSet.Internal.ERI4;
// import qchem.BasisSet.Internal.ERI3;
// import qchem.BasisSet.Internal.IntegralEnums;

namespace std
{
    using I1C=BasisSet1::IntegralsCache_Base::I1C;
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

    using I2C=BasisSet1::IntegralsCache_Base::I2C;
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
    using I2x=BasisSet1::IntegralsCache_Base::I2x;
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

    using I2n=BasisSet1::IntegralsCache_Base::I2n;
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

    using I3C=BasisSet1::IntegralsCache_Base::I3C;
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

    using I4C=BasisSet1::IntegralsCache_Base::I4C;
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
namespace BasisSet1
{
// template  <> IntegralsCache<double>* IntegralsCache<double>::theGlobalCache;

// size_t Cache41::Register(const std::string& bf_id)
// {
//     size_t index;
//     auto i=itsUniqueBFs.find(bf_id);
//     if (i==itsUniqueBFs.end())
//         itsUniqueBFs[bf_id]=index=itsUniqueBFs.size();
//     else
//         index=i->second;
//     return index;
// }

//
// Specialize std::formatter for I1C,I2C
//
 


std::string Hit(bool hit)
{
    return hit ? "Hit " : "Miss";
}
std::ostream& operator<<(std::ostream& os,IntegralsCache_Base::IBS_ID_t id)
{
    return os << std::get<0>(id) << " " << std::get<1>(id);
}
template <class T> IntegralsCache_RAM<T>::IntegralsCache_RAM(bool makelog) 
    : itsMakeLog(makelog) 
    {
        if  (itsMakeLog)
            itsLogger=std::ofstream("cache.log");
    };

template <class T> bool IntegralsCache_RAM<T>::Has(Ix1 ix,const IBS_ID_t& id) const
{
    bool ret=false;
    std::string type;
    if (std::holds_alternative<I1C>(ix))
    {
        auto i1c=std::get<I1C>(ix);
        itsLastKey1=key1_t(i1c,id);
        auto i=itsVecs.find(itsLastKey1);
        its1CIterator=i;
        ret=i!=itsVecs.end(); 
        type=std::format("I1C {:<12}", i1c);
    }
    else if (std::holds_alternative<I2C>(ix))
    {
        auto i2c=std::get<I2C>(ix);
        itsLastKey2=key2_t(i2c,id);
        auto i=itsSMats.find(itsLastKey2);
        its2CIterator=i;
        ret=i!=itsSMats.end(); 
        type=std::format("I2C {:<12}", i2c);
    }
    else
    {
        std::cerr << "IntegralsCache_RAM: Unhandled integral type alternative " << std::endl;
        assert(false);
        exit(-1);
    }
    if  (itsMakeLog && !ret)
        itsLogger << "Ix1 " << type << " cache " << Hit(ret) << " a=" << id << std::endl;

    return ret;
}

template <class T> bool IntegralsCache_RAM<T>::Has(Ix2 ix,const IBS_ID_t& ida,const IBS_ID_t& idb) const
{
    bool ret=false;
    std::string type;
    if (std::holds_alternative<I2x>(ix))
    {
        auto i2x=std::get<I2x>(ix);
        itsLastKeyx=keyx_t(i2x,ida,idb);
        its2xIterator=itsMats.find(itsLastKeyx);
        ret=its2xIterator!=itsMats.end(); 
        type=std::format("I2x {:<12}", i2x);
    }
    else if (std::holds_alternative<I3C>(ix))
    {
        auto i3c=std::get<I3C>(ix);
        itsLastKey3=key3_t(i3c,ida,idb);
        its3CIterator=itsERI3s.find(itsLastKey3);
        ret=its3CIterator!=itsERI3s.end(); 
        type=std::format("I3C {:<12}", i3c);
    }
    else if (std::holds_alternative<I4C>(ix))
    {
        itsLastKey4a=ida;
        itsLastKey4b=idb;
        auto i4c=std::get<I4C>(ix);
        type=std::format("I4C {:<12}", i4c);
        switch (i4c)
        {
            case IntegralsCache<T>::I4C::Direct:
            {
                auto ia=Jac.find(ida);
                if (ia!=Jac.end())
                {
                    its4CIterator=ia->second.find(idb);
                    ret=its4CIterator!=ia->second.end(); 
                }
                break;
            }
            case IntegralsCache<T>::I4C::Exchange:
            {
                auto ia=Kab.find(ida);
                if (ia!=Kab.end())
                {
                    its4CIterator=ia->second.find(idb);
                    ret=its4CIterator!=ia->second.end(); 
                }
                break;
            }
        } //switch
        
    }
    else
    {
        std::cerr << "IntegralsCache_RAM: Unhandled integral type alternative " << std::endl;
        assert(false);
        exit(-1);
    }
    if  (itsMakeLog && !ret)
        itsLogger << "Ix2 " << type << " cache " << Hit(ret) << " a=" << ida << " b=" << idb << std::endl;

    return ret;
}

template <class T> bool IntegralsCache_RAM<T>::Has(I2n i2n,const IBS_ID_t& IBS_id,const Cluster_ID_t& cluster_id) const
{
    itsLastKeyn=keyn_t(IBS_id,cluster_id);
    auto i=itsNMats.find(itsLastKeyn);
    its2CIterator=i;
    bool ret=i!=itsNMats.end();
    std::string type=std::format("    {:<12}", i2n);
    if  (itsMakeLog && !ret)
        itsLogger << "I2n " << type << " cache " << Hit(ret) << " a=" << IBS_id << " cluster=" << cluster_id << std::endl;
    return ret; 
}

template <class T>  bool IntegralsCache_RAM<T>::Has(I1C i1c,const IBS_ID_t& id,const Mesh_ID_t& mid) const
{
    itsLastKey1m=key1m_t(i1c,id,mid);
    auto i=itsmVecs.find(itsLastKey1m);
    its1CIterator=i;
    bool ret=i!=itsmVecs.end();
    std::string type=std::format("    {:<12}", i1c);
    if  (itsMakeLog && !ret)
        itsLogger << "1Cm " << type << " cache " << Hit(ret) << " a=" << id << " mesh=" << mid << std::endl;
    return ret; 
}
template <class T>  bool IntegralsCache_RAM<T>::Has(I2x i2x,const IBS_ID_t& a,const IBS_ID_t& b,const Mesh_ID_t& mid) const
{
    itsLastKey2xm=key2xm_t(i2x,a,b,mid);
    auto i=itsmMats.find(itsLastKey2xm);
    its2xmIterator=i;
    bool ret=i!=itsmMats.end(); 
    std::string type=std::format("    {:<12}", i2x);
    if  (itsMakeLog && !ret)
        itsLogger << "2xm " << type << " cache " << Hit(ret) << " a=" << a << " b=" << b << " mesh=" << mid << std::endl;
    return ret;
}


template <class T>  const rvec_t& IntegralsCache_RAM<T>::GetVec() const
{
    if (std::holds_alternative<typename map1_t::const_iterator>(its1CIterator))
        return std::get<typename  map1_t::const_iterator>(its1CIterator)->second;
    else //
        return std::get<typename map1m_t::const_iterator>(its1CIterator)->second;

}
template <class T>  const smat_t<T>& IntegralsCache_RAM<T>::GetSMat() const
{
    if (std::holds_alternative<typename map2_t::const_iterator>(its2CIterator))
        return std::get<typename map2_t::const_iterator>(its2CIterator)->second;
    else //if (std::holds_alternative<typename mapn_t::const_iterator>(its2CnIterator))
        return std::get<typename mapn_t::const_iterator>(its2CIterator)->second;
}
template <class T>  const mat_t<T>& IntegralsCache_RAM<T>::GetMat() const
{
    return its2xIterator->second;
}
template <class T>  const ERI3<T>& IntegralsCache_RAM<T>::GetERI3() const
{
    return its3CIterator->second;
}
template <class T>  const ERI4& IntegralsCache_RAM<T>::GetERI4() const
{
    return its4CIterator->second;
}



template <class T>  const rvec_t&  IntegralsCache_RAM<T>::Set(const rvec_t& v)
{
    if (std::holds_alternative<typename map1_t::const_iterator>(its1CIterator))
    {
        const auto [iterator, success]=itsVecs.insert({itsLastKey1,v});  //Non-nuclear
        assert(success);
        its1CIterator=iterator;
        return iterator->second; 
    }
    else //if (std::holds_alternative<typename map1m_t::const_iterator>(its1CmIterator))
    {
        const auto [iterator, success]=itsmVecs.insert({itsLastKey1m,v});  //Non-nuclear
        assert(success);
        its1CIterator=iterator;
        return iterator->second; 
    }
}

template <class T> const smat_t<T>& IntegralsCache_RAM<T>::Set(const smat_t<T>& m)
{
    if (std::holds_alternative<typename map2_t::const_iterator>(its2CIterator))
    {
        const auto [iterator, success]=itsSMats.insert({itsLastKey2,m});  //Non-nuclear
        assert(success);
        its2CIterator=iterator;
        return iterator->second; 
    }
    else //if (std::holds_alternative<typename mapn_t::const_iterator>(its2CnIterator))
    {
        const auto [iterator, success]=itsNMats.insert({itsLastKeyn,m}); //Nuclear
        assert(success);
        its2CIterator=iterator;
        return iterator->second;
    }
    
}

template <class T> const  mat_t<T>& IntegralsCache_RAM<T>::Set(const  mat_t<T>& m)
{
    const auto [iterator, success]=itsMats.insert({itsLastKeyx,m});  //Non-nuclear
    assert(success);
    its2xIterator=iterator;
    return iterator->second; 
}

template <class T> const   ERI3<T>& IntegralsCache_RAM<T>::Set(const   ERI3<T>& eri3)
{
    const auto [iterator, success]=itsERI3s.insert({itsLastKey3,eri3});  
    assert(success);
    its3CIterator=iterator;
    return iterator->second; 
}

template <class T> const ERI4     &  IntegralsCache_RAM<T>::SetDirect(const   ERI4   & eri4)
{
    const auto [iterator, success]=Jac[itsLastKey4a].insert({itsLastKey4b,eri4});
    assert(success);
    its4CIterator=iterator;
    return iterator->second;
}
template <class T> const ERI4     &  IntegralsCache_RAM<T>::SetExchange(const   ERI4   & eri4)
{
    const auto [iterator, success]=Kab[itsLastKey4a].insert({itsLastKey4b,eri4});
    assert(success);
    its4CIterator=iterator;
    return iterator->second;
}

template struct IntegralsCache<double>;
template struct IntegralsCache_RAM<double>;

} //namespace