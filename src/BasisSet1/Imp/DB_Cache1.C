// File: BasisSet/Imp/DB_Cache1.C Global integrals cache allow data sharing between separate runs.
module;
#include <map>
#include <cassert>
#include <vector>
#include <string>
#include <variant>
#include <iostream>
module qchem.BasisSet1.DB_Cache;
// import qchem.BasisSet.Internal.ERI4;
// import qchem.BasisSet.Internal.ERI3;
// import qchem.BasisSet.Internal.IntegralEnums;

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

template <class T> bool IntegralsCache_RAM<T>::Has(Ix1 ix,const IBS_ID_t& id) const
{
    bool ret=false;
    if (std::holds_alternative<I1C>(ix))
    {
        itsLastKey1=key1_t(std::get<I1C>(ix),id);
        its1CIterator=itsVecs.find(itsLastKey1);
        ret=its1CIterator!=itsVecs.end(); 
    }
    else if (std::holds_alternative<I2C>(ix))
    {
        itsLastKey2=key2_t(std::get<I2C>(ix),id);
        auto i=itsSMats.find(itsLastKey2);
        its2CnIterator=i;
        ret=i!=itsSMats.end(); 
    }
    else
    {
        std::cerr << "IntegralsCache_RAM: Unhandled integral type alternative " << std::endl;
        assert(false);
        exit(-1);
    }
    return ret;
}

template <class T> bool IntegralsCache_RAM<T>::Has(Ix2 ix,const IBS_ID_t& ida,const IBS_ID_t& idb) const
{
    bool ret=false;
    if (std::holds_alternative<I2x>(ix))
    {
        itsLastKeyx=keyx_t(std::get<I2x>(ix),ida,idb);
        its2xIterator=itsMats.find(itsLastKeyx);
        ret=its2xIterator!=itsMats.end(); 
    }
    else if (std::holds_alternative<I3C>(ix))
    {
        itsLastKey3=key3_t(std::get<I3C>(ix),ida,idb);
        its3CIterator=itsERI3s.find(itsLastKey3);
        ret=its3CIterator!=itsERI3s.end(); 
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
    return ret;
}

template <class T> bool IntegralsCache_RAM<T>::Has(I2n,const IBS_ID_t& IBS_id,const Cluster_ID_t& cluster_id) const
{
    itsLastKeyn=keyn_t(IBS_id,cluster_id);
    auto i=itsNMats.find(itsLastKeyn);
    its2CnIterator=i;
    return i!=itsNMats.end(); 
}

template <class T>  const rvec_t& IntegralsCache_RAM<T>::GetVec() const
{
    assert(its1CIterator!=itsVecs.end());
    return its1CIterator->second;
}
template <class T>  const smat_t<T>& IntegralsCache_RAM<T>::GetSMat() const
{
    if (std::holds_alternative<typename map2_t::const_iterator>(its2CnIterator))
        return std::get<typename map2_t::const_iterator>(its2CnIterator)->second;
    else //if (std::holds_alternative<typename mapn_t::const_iterator>(its2CnIterator))
        return std::get<typename mapn_t::const_iterator>(its2CnIterator)->second;
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



template <class T>  void IntegralsCache_RAM<T>::Set(const rvec_t& v)
{
    itsVecs[itsLastKey1]=v;
}

template <class T> const smat_t<T>& IntegralsCache_RAM<T>::Set(const smat_t<T>& m)
{
    // auto iterator= std::holds_alternative<typename map2_t::const_iterator>(its2CnIterator)
    //     ? itsSMats.insert({itsLastKey2,m})  //Non-nuclear
    //     : itsNMats.insert({itsLastKeyn,m}); //Nuclear
    // its2CnIterator=iterator.first;
    // return iterator.second;
    if (std::holds_alternative<typename map2_t::const_iterator>(its2CnIterator))
    {
        const auto [iterator, success]=itsSMats.insert({itsLastKey2,m});  //Non-nuclear
        assert(success);
        its2CnIterator=iterator;
        return iterator->second; 
    }
    else //if (std::holds_alternative<typename mapn_t::const_iterator>(its2CnIterator))
    {
        const auto [iterator, success]=itsNMats.insert({itsLastKeyn,m}); //Nuclear
        assert(success);
        its2CnIterator=iterator;
        return iterator->second;
    }
    
}

template <class T> void IntegralsCache_RAM<T>::Set(const  mat_t<T>& m)
{
    itsMats[itsLastKeyx]=m;
}

template <class T> void IntegralsCache_RAM<T>::Set(const   ERI3<T>& eri3)
{
    itsERI3s[itsLastKey3]=eri3;
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
