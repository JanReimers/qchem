// File: BasisSet/Internal/DB_Cache.C Cache integrals for a full basis set.
module;
#include <map>
#include <cassert>
#include <vector>
export module qchem.BasisSet.Internal.DB_Cache;
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Internal.ERI3;
import qchem.BasisSet.Internal.IntegralEnums;
import qchem.Orbital_HF_IBS;
import Common.UniqueID;
 
//
//  Each (full) Basis set has one of these integral caches built in through inheritance.
//  Each irrep basis set (IBS) hold a pointer to this cache.  Cached integrals are distinguised based
//  on the unique ID of each IBS.  
//
export template  <class T> struct DB_cache  
{
    typedef UniqueID::IDtype IDType;
    typedef std::map<IDType,std::map<IDType,ERI4> > erij_t;
    typedef std::tuple<qchem::IType2C,IDType> id2c_t; //Integral key for one IBS.
    typedef std::tuple<qchem::IType2C,IDType,IDType> idx_t; //Cross Integral key for two IBSs.
    typedef std::tuple<qchem::IType3C,IDType,IDType> id3c_t; //Ingteral key for 3 center integrals.
    
    mutable std::map<id2c_t ,rsmat_t> itsbSMats; 
    mutable std::map< idx_t , rmat_t> itsbMats; 
    mutable std::map<id2c_t , rvec_t> itsbVecs; 
    mutable std::map<id3c_t ,ERI3<T>> itsERI3s; 
    mutable erij_t Jac,Kab; //4 center, 2 electron integrals.  Direct: Jac and Exchange: Kab
};

export template <class T> class DB_BS_2E : public virtual Integrals_BS_2E<T>, public DB_cache<T>
{
    typedef UniqueID::IDtype IDType;   
public: 
    virtual ERI4 Direct  (IDType a,IDType c) const;
    virtual ERI4 Exchange(IDType a,IDType b) const;
protected:
    void Append(const Orbital_HF_IBS<T>*);
    virtual ERI4 MakeDirect  (const Orbital_HF_IBS<T>* a, const Orbital_HF_IBS<T>* c) const=0;
    virtual ERI4 MakeExchange(const Orbital_HF_IBS<T>* a, const Orbital_HF_IBS<T>* b) const=0;
private:
    //! Internally called once to build direct and exchange supermatrix tables.
    virtual void MakeDirect  () const; 
    virtual void MakeExchange() const; 
    std::vector<const Orbital_HF_IBS<T>*> itsIrreps; //Used for 2-electron integrals.
    using DB_cache<T>::Jac;
    using DB_cache<T>::Kab;
};
