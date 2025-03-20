// File: HeapDB.H  Implement an integral data base that uses heap for storage.
#ifndef _HeapDB_H_
#define _HeapDB_H_



#include <IEClient.H>
#include <AnalyticIE.H>
#include "Imp/Containers/ERI4.H"
#include <map>
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
 
template <class T> class DB_Common : virtual public Integrals_Base<T>
{
protected:
    typedef UniqueID::IDtype IDType;
    virtual IDType GetID() const=0;
    DB_Common(const DB_cache<T>* db) : itsCache(db) {assert(itsCache);}
    const DB_cache<T>* itsCache;
};

// It would be nice use virtual inheretance from DB_Common<T>, but that opens a contructor hornets nest.
template <class T> class DB_Overlap : public DB_Common<T>, public virtual Integrals_Overlap<T>
{    
protected:
    DB_Overlap(const DB_cache<T>* db) : DB_Common<T>(db) {};
    virtual typename Integrals_Base<T>::SMat_ref Overlap() const;
    virtual typename Integrals_Base<T>::SMat MakeOverlap() const=0;
};
template <class T> class DB_Kinetic : public DB_Common<T>, public virtual Integrals_Kinetic<T>
{    
protected:
    DB_Kinetic(const DB_cache<T>* db) : DB_Common<T>(db) {};
    virtual typename Integrals_Base<T>::SMat_ref Kinetic() const;
    virtual typename Integrals_Base<T>::SMat MakeKinetic() const=0;
};
template <class T> class DB_Nuclear : public DB_Common<T>, public virtual Integrals_Nuclear<T>
{    
protected:
    DB_Nuclear(const DB_cache<T>* db) : DB_Common<T>(db) {};
    virtual typename Integrals_Base<T>::SMat_ref Nuclear(const Cluster*) const;
    virtual typename Integrals_Base<T>::SMat MakeNuclear(const Cluster*) const=0;
};
template <class T> class DB_XKinetic : public DB_Common<T>, public virtual Integrals_XKinetic<T>
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

template <class T> class DB_DFT 
    : virtual public Integrals_DFT<T>
    , public DB_Common<T>
{
    typedef typename Integrals_Base<T>::fbs_t fbs_t;
    typedef typename Integrals_Base<T>::ERI3 ERI3;
protected:
    DB_DFT(const DB_cache<T>* db) : DB_Common<T>(db) {}
    virtual const ERI3& Overlap3C  (const fbs_t& c) const; //<ab|c>
    virtual const ERI3& Repulsion3C(const fbs_t& c) const; //<a(1)b(1)|1/r12|c(2)>
    virtual ERI3 MakeOverlap3C  (const fbs_t& c) const=0;
    virtual ERI3 MakeRepulsion3C(const fbs_t& c) const=0;
};


template <class T> class DB_BS_2E : public virtual Integrals_BS_2E<T>, public DB_cache<T>
{
    typedef UniqueID::IDtype IDType;   
public: 
    virtual ERI4 Direct  (IDType a,IDType c) const;
    virtual ERI4 Exchange(IDType a,IDType b) const;
protected:
    void Append(const IrrepIEClient*);
    virtual ERI4 MakeDirect  (const IrrepIEClient* a, const IrrepIEClient* c) const=0;
    virtual ERI4 MakeExchange(const IrrepIEClient* a, const IrrepIEClient* b) const=0;
private:
    //! Internally called once to build direct and exchange supermatrix tables.
    virtual void MakeDirect  () const; 
    virtual void MakeExchange() const; 
    std::vector<const IrrepIEClient*> itsIrreps; //Used for 2-electron integrals.
    using DB_cache<T>::Jac;
    using DB_cache<T>::Kab;
};
template <class T> class DB_2E :  virtual public Integrals_2E<T>
{
    typedef Integrals_Base<T> Base;
    typedef typename Base::obs_t obs_t;
public:
    typedef UniqueID::IDtype IDType;
    virtual ERI4 Direct  (const obs_t& c) const;
    virtual ERI4 Exchange(const obs_t& b) const;
protected:
    DB_2E() : itsDB_BS_2E(0) {};
    DB_2E(const DB_BS_2E<T>* db);
private:
    virtual IDType GetID() const=0;
    const DB_BS_2E<T>* itsDB_BS_2E; //Database of all supermatrix tables.
};

class DB_Fit 
    : public virtual FitIntegrals_1E
    , public DB_Common<double>
{
protected:
    DB_Fit(const DB_cache<double>* db) : DB_Common<double>(db) {};

    using Integrals_Overlap<double>::Overlap; 
    virtual Vec_ref  Charge   () const;   
    virtual SMat_ref Repulsion() const;
    virtual  Mat_ref Repulsion(const fbs_t&) const;
    virtual SMat_ref InvOverlap(const LAParams&) const;
    virtual SMat_ref InvRepulsion(const LAParams&) const;
    virtual  Vec_ref Norm   (const Mesh*        ) const; //Numerical .
    virtual  Vec_ref Charge (const Mesh*        ) const; //Numerical .
    virtual  Mat_ref Overlap(const Mesh*,const fbs_t& b) const; //Numerical X overlap.


private:
    // One time calls to un-buffered integral calculations.
    using FitIntegrals_1E::MakeCharge;
    virtual Vec  MakeCharge() const=0;
    // virtual SMat MakeOverlap() const=0;
    virtual SMat MakeRepulsion() const=0;
    virtual  Mat MakeRepulsion(const fbs_t&) const=0;
    
    //! \brief Return the Penrose inverse of a symmetric matrix using SVD decomposition
    //! If \f$ S=UsV^{\dagger} \f$, then \f$ S^{-1}=V\frac{1}{s}U^{\dagger} \f$
    static  SMat MakeInverse  (SMat_ref,const LAParams&); //Numerically stable algo required.

};

template <class T> class DB_RKB 
    : virtual public Integrals_RKB<T>
    , public DB_Overlap<T>
    , public DB_Kinetic<T>
    , public DB_Nuclear<T>
    , public DB_RestMass<T>
{    
protected:
    DB_RKB(const DB_cache<T>* db) :  DB_Overlap<T>(db), DB_Kinetic<T>(db), DB_Nuclear<T>(db), DB_RestMass<T>(db) {}; 
   
};
template <class T> class DB_RKBL 
    : virtual public Integrals_RKBL<T>
    , public DB_Overlap<T>
    , public DB_XKinetic<T>
    , public DB_Nuclear<T>
{
protected:
    DB_RKBL(const DB_cache<T>* db) : DB_Overlap<T>(db), DB_XKinetic<T>(db), DB_Nuclear<T>(db) {};
};
template <class T> class DB_RKBS 
    : virtual public Integrals_RKBS<T>
    , public DB_Kinetic<T>
    , public DB_Nuclear<T>
{
protected:
    DB_RKBS(const DB_cache<T>* db) :  DB_Kinetic<T>(db), DB_Nuclear<T>(db) {};  
   
};

#endif
