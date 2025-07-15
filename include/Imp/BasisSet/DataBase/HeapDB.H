// File: HeapDB.H  Implement an integral data base that uses heap for storage.
#ifndef _HeapDB_H_
#define _HeapDB_H_



#include <map>
#include "../ERI4.H"
#include "oml/vector.h"
#include <BasisSet/Integrals.H>
#include "../IEClient.H"

template  <class T> class DB_cache  : virtual public Integrals_Base<T>
{
    typedef UniqueID::IDtype IDType;
    typedef Integrals_Base<T> Base;
    typedef typename Base::SMat SMat;    
    typedef typename Base::Mat Mat;    
    typedef typename Base::Vec Vec;    
    typedef typename Base::ERI3 ERI3;    
public:
    typedef std::map<IDType,std::map<IDType,ERI4> > erij_t;
    typedef std::tuple<qchem::IType2C,IDType> id2c_t;
    typedef std::tuple<qchem::IType2C,IDType,IDType> idx_t;
    typedef std::tuple<qchem::IType3C,IDType,IDType> id3c_t;
    
    mutable std::map<id2c_t ,SMat> itsSMats; 
    mutable std::map< idx_t , Mat> itsMats; 
    mutable std::map<id2c_t , Vec> itsVecs; 
    mutable std::map<id3c_t ,ERI3> itsERI3s; 
    mutable erij_t Jac,Kab;
};
 
template <class T> class DB_Common 
    : virtual public Integrals_Base<T>
    , virtual public UniqueID
{
protected:
    DB_Common(const DB_cache<T>* db) : itsCache(db) {assert(itsCache);}
    const DB_cache<T>* itsCache;
};

// It would be nice use virtual inheretance from DB_Common<T>, but that opens a contructor hornets nest.
template <class T> class DB_Overlap  : public DB_Common<T>, public virtual Integrals_Overlap<T>
{    
protected:
    DB_Overlap(const DB_cache<T>* db) : DB_Common<T>(db) {};
    virtual typename Integrals_Base<T>::SMat_ref Overlap() const;
    virtual typename Integrals_Base<T>::SMat MakeOverlap() const=0;
};
template <class T> class DB_Kinetic    : public DB_Common<T>, public virtual Integrals_Kinetic<T>
{    
protected:
    DB_Kinetic(const DB_cache<T>* db) : DB_Common<T>(db) {};
    virtual typename Integrals_Base<T>::SMat_ref Kinetic() const;
    virtual typename Integrals_Base<T>::SMat MakeKinetic() const=0;
};
template <class T> class DB_Nuclear  : public DB_Common<T>, public virtual Integrals_Nuclear<T>
{    
protected:
    DB_Nuclear(const DB_cache<T>* db) : DB_Common<T>(db) {};
    virtual typename Integrals_Base<T>::SMat_ref Nuclear(const Cluster*) const;
    virtual typename Integrals_Base<T>::SMat MakeNuclear(const Cluster*) const=0;
};
template <class T> class DB_XKinetic   : public DB_Common<T>, public virtual Integrals_XKinetic<T>
{    
protected:
    DB_XKinetic(const DB_cache<T>* db) : DB_Common<T>(db) {};
    virtual typename Integrals_Base<T>::Mat_ref Kinetic(const Orbital_RKBS_IBS<T>* rkbs) const;
    virtual typename Integrals_Base<T>::Mat MakeKinetic(const Orbital_RKBS_IBS<T>* rkbs) const=0;
};
template <class T> class DB_RestMass : public DB_Common<T>, public virtual Integrals_RestMass<T>
{    
protected:
    DB_RestMass(const DB_cache<T>* db) : DB_Common<T>(db) {};
    virtual typename Integrals_Base<T>::SMat_ref RestMass() const;
    virtual typename Integrals_Base<T>::SMat MakeRestMass() const=0;
};




#endif
