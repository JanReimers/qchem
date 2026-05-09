// File: BasisSet/Imp/DB_Cache1.C Global integrals cache allow data sharing between separate runs.
module;
#include <map>
#include <cassert>
#include <vector>
#include <string>
#include <variant>
#include <iostream>
#include <fstream>
module qchem.BasisSet1.DB_Cache;

// import qchem.BasisSet.Internal.ERI4;
// import qchem.BasisSet.Internal.ERI3;
// import qchem.BasisSet.Internal.IntegralEnums;

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
        if (itsMakeLog)
            itsLogger=std::ofstream("cache.log");
    };

template <class T> bool IntegralsCache_RAM<T>::Has(Ix1 ix,const IBS_ID_t& id) const
{
    bool ret=false;
    std::string type;
    if (std::holds_alternative<I1C>(ix))
    {
        itsLastKey1=key1_t(std::get<I1C>(ix),id);
        auto i=itsVecs.find(itsLastKey1);
        its1CIterator=i;
        ret=i!=itsVecs.end(); 
        type="1C ";
    }
    else if (std::holds_alternative<I2C>(ix))
    {
        itsLastKey2=key2_t(std::get<I2C>(ix),id);
        auto i=itsSMats.find(itsLastKey2);
        its2CIterator=i;
        ret=i!=itsSMats.end(); 
        type="2C ";
    }
    else
    {
        std::cerr << "IntegralsCache_RAM: Unhandled integral type alternative " << std::endl;
        assert(false);
        exit(-1);
    }
    if (itsMakeLog)
        itsLogger << "Ix1 " << type << " cache " << Hit(ret) << " " << id << std::endl;

    return ret;
}

template <class T> bool IntegralsCache_RAM<T>::Has(Ix2 ix,const IBS_ID_t& ida,const IBS_ID_t& idb) const
{
    bool ret=false;
    std::string type;
    if (std::holds_alternative<I2x>(ix))
    {
        itsLastKeyx=keyx_t(std::get<I2x>(ix),ida,idb);
        its2xIterator=itsMats.find(itsLastKeyx);
        ret=its2xIterator!=itsMats.end(); 
        type="2Cx";
    }
    else if (std::holds_alternative<I3C>(ix))
    {
        itsLastKey3=key3_t(std::get<I3C>(ix),ida,idb);
        its3CIterator=itsERI3s.find(itsLastKey3);
        ret=its3CIterator!=itsERI3s.end(); 
        type="3C ";
    }
    else if (std::holds_alternative<I4C>(ix))
    {
        itsLastKey4a=ida;
        itsLastKey4b=idb;
        switch (std::get<I4C>(ix))
        {
            case IntegralsCache<T>::I4C::Direct:
            {
                auto ia=Jac.find(ida);
                if (ia!=Jac.end())
                {
                    its4CIterator=ia->second.find(idb);
                    ret=its4CIterator!=ia->second.end(); 
                }
                type="Jac";
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
                type="Kab";
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
    assert(type.size()==3);
    if (itsMakeLog)
        itsLogger << "Ix2 " << type << " cache " << Hit(ret) << " a=" << ida << " b=" << idb << std::endl;

    return ret;
}

template <class T> bool IntegralsCache_RAM<T>::Has(I2n,const IBS_ID_t& IBS_id,const Cluster_ID_t& cluster_id) const
{
    itsLastKeyn=keyn_t(IBS_id,cluster_id);
    auto i=itsNMats.find(itsLastKeyn);
    its2CIterator=i;
    bool ret=i!=itsNMats.end();
    if (itsMakeLog)
        itsLogger << "I2n     cache " << Hit(ret) << " " << IBS_id << " cluster=" << cluster_id << std::endl;
    return ret; 
}

template <class T>  bool IntegralsCache_RAM<T>::Has(I1C ix,const IBS_ID_t& id,const Mesh_ID_t& mid) const
{
    itsLastKey1m=key1m_t(ix,id,mid);
    auto i=itsmVecs.find(itsLastKey1m);
    its1CIterator=i;
    bool ret=i!=itsmVecs.end();
    if (itsMakeLog)
        itsLogger << "1Cm     cache " << Hit(ret) << " " << id << " mesh=" << mid << std::endl;
    return ret; 
}
template <class T>  bool IntegralsCache_RAM<T>::Has(I2x ix,const IBS_ID_t& a,const IBS_ID_t& b,const Mesh_ID_t& mid) const
{
    itsLastKey2xm=key2xm_t(ix,a,b,mid);
    auto i=itsmMats.find(itsLastKey2xm);
    its2xmIterator=i;
    bool ret=i!=itsmMats.end(); 
    if (itsMakeLog)
        itsLogger << "I2xm    cache " << Hit(ret) << " a=" << a << " b=" << b << " mesh=" << mid << std::endl;
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