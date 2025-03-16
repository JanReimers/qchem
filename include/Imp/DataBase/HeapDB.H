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
public:
    typedef std::map<IDType,std::map<IDType,ERI4> > erij_t;
    typedef std::tuple<qchem::IType2C,IDType> id2c_t;
    typedef std::tuple<qchem::IType2C,IDType,IDType> idx_t;
    
    mutable std::map<id2c_t ,SMat> itsSMats; 
    mutable std::map< idx_t , Mat> itsMats; 
    mutable std::map<id2c_t , Vec> itsVecs; 
    mutable erij_t Jac,Kab;
};
 
template <class T> class DB_1E :
    virtual public Integrals_1E<T>
{    
    typedef Integrals_Base<T> Base;
    typedef typename Base::SMat_ref SMat_ref;
public:
    //  This is all buffered.
    virtual SMat_ref Overlap() const;
    virtual SMat_ref Kinetic() const;
    virtual SMat_ref Nuclear(const Cluster*) const;
protected:
    typedef typename Base::SMat SMat;
    DB_1E(const DB_cache<T>* db) : itsCache(db) {assert(itsCache);}
private:
    typedef typename DB_cache<T>::id2c_t id2c_t;
    typedef UniqueID::IDtype IDType;
    virtual IDType GetID() const=0;
    // One time calls to un-buffered integral calculations.
    virtual SMat MakeOverlap() const=0;
    virtual SMat MakeKinetic() const=0;
    virtual SMat MakeNuclear(const Cluster*) const=0;
    const DB_cache<T>* itsCache;


};

template <class T> class DB_DFT :
    virtual public Integrals_DFT<T>
{
    typedef Integrals_Base<T> Base;
    typedef typename Base::fbs_t fbs_t;
    typedef typename Base::SMat SMat;
    typedef typename Integrals_DFT<T>::ERI3 ERI3;
public:
    virtual const ERI3& Overlap3C  (const fbs_t& c) const; //<ab|c>
    virtual const ERI3& Repulsion3C(const fbs_t& c) const; //<a(1)b(1)|1/r12|c(2)>
    protected:
    DB_DFT(const DB_cache<T>* db) : itsCache(db) {assert(itsCache);}
private:
    typedef UniqueID::IDtype IDType;
    virtual IDType GetID() const=0;
    virtual ERI3 MakeOverlap3C  (const fbs_t& c) const=0;
    virtual ERI3 MakeRepulsion3C(const fbs_t& c) const=0;
    typedef std::tuple<qchem::IType3C,IDType,IDType> id3c_t;
    mutable std::map<id3c_t ,ERI3> itsBuffer; 
    const DB_cache<T>* itsCache;
};


template <class T> class DB_BS_2E : public virtual Integrals_BS_2E<T>, public DB_cache<T>
{
    typedef UniqueID::IDtype IDType;    
public:
    virtual    ERI4 Direct  (IDType a,IDType c) const;
    virtual    ERI4 Exchange(IDType a,IDType b) const;
protected:
    virtual void Append(const IrrepIEClient*);
private:
    //! Internally called once to build direct and exchange supermatrix tables.
    virtual void MakeDirect  () const; 
    virtual void MakeExchange() const; 
    virtual ERI4 MakeDirect  (const IrrepIEClient* a, const IrrepIEClient* c) const=0;
    virtual ERI4 MakeExchange(const IrrepIEClient* a, const IrrepIEClient* b) const=0;

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

class DB_Fit : public virtual FitIntegrals_1E<double>
{
    typedef Integrals_Base<double> Base;
    typedef typename Base::Mat_ref Mat_ref;
    typedef typename Base::SMat_ref SMat_ref;
    typedef typename Base::Vec_ref Vec_ref;
    typedef typename Base::fbs_t fbs_t;
public:
    virtual Vec_ref  Charge   () const;   
    virtual SMat_ref Overlap  () const;
    virtual SMat_ref Repulsion() const;
    virtual  Mat_ref Repulsion(const fbs_t&) const;
    virtual SMat_ref InvOverlap(const LAParams&) const;
    virtual SMat_ref InvRepulsion(const LAParams&) const;
    virtual  Vec_ref Norm   (const Mesh*        ) const; //Numerical .
    virtual  Vec_ref Charge (const Mesh*        ) const; //Numerical .
    virtual  Mat_ref Overlap(const Mesh*,const fbs_t& b) const; //Numerical X overlap.

protected:
    DB_Fit(const DB_cache<double>* db) : itsCache(db) {assert(itsCache);}

    typedef Base::SMat SMat;
    typedef Base::Mat  Mat;
    typedef Base::Vec  Vec;
private:
    typedef UniqueID::IDtype IDType;
    virtual IDType GetID() const=0;
    // One time calls to un-buffered integral calculations.
    using FitIntegrals_1E<double>::MakeOverlap;
    using FitIntegrals_1E<double>::MakeCharge;
    virtual Vec  MakeCharge() const=0;
    virtual SMat MakeOverlap() const=0;
    virtual SMat MakeRepulsion() const=0;
    virtual  Mat MakeRepulsion(const fbs_t&) const=0;
    
    //! \brief Return the Penrose inverse of a symmetric matrix using SVD decomposition
    //! If \f$ S=UsV^{\dagger} \f$, then \f$ S^{-1}=V\frac{1}{s}U^{\dagger} \f$
    static  SMat MakeInverse  (SMat_ref,const LAParams&); //Numerically stable algo required.

    typedef typename DB_cache<double>::id2c_t id2c_t;
    typedef typename DB_cache<double>::idx_t idx_t;
    const DB_cache<double>* itsCache;
};

template <class T> class DB_RKB :
    virtual public Integrals_RKB<T>
{    
public:
    typedef UniqueID::IDtype IDType;
    typedef Integrals_RKB<T> Base;
    typedef Matrix<T> Mat;
    typedef SMatrix<T> SMat;
    typedef typename Base::Mat_ref Mat_ref;
    typedef typename Base::SMat_ref SMat_ref;
    virtual SMat_ref Kinetic() const;
    virtual SMat_ref RestMass() const;
private:
    virtual IDType GetID() const=0;
    virtual SMat MakeKinetic() const=0;
    virtual SMat MakeRestMass() const=0;
    typedef std::tuple<qchem::IType2C,IDType> id2c_t;
    mutable std::map<id2c_t ,SMat> itsBuffer;
};

template <class T> class DB_RKBL :
    virtual public Integrals_RKBL<T>
{
public:
    typedef UniqueID::IDtype IDType;
    typedef Integrals_RKBL<T> Base;
    typedef typename Base::SMat_ref SMat_ref;
    typedef typename Base::Mat_ref Mat_ref;
    typedef typename Base::SMat SMat;
    typedef typename Base::Mat  Mat;

    virtual SMat_ref Overlap() const;
    virtual Mat_ref  Kinetic(const Orbital_RKBS_IBS<T>* rkbs) const;
    virtual SMat_ref Nuclear(const Cluster*) const;
private:

    virtual IDType GetID() const=0;
    virtual SMat MakeOverlap() const;
    virtual Mat  MakeKinetic(const Orbital_RKBS_IBS<T>* rkbs) const;
    virtual SMat MakeNuclear(const Cluster*) const;
    SMat MakeIntegrals(qchem::IType,const Cluster* cl=0) const;
    virtual T    Integral(qchem::IType,double ea , double eb,size_t l) const=0;

    typedef std::tuple<qchem::IType2C,IDType> id2c_t;
    mutable std::map<id2c_t ,SMat> itsBuffer; 
    mutable std::map<id2c_t ,Mat> itsBufferX; 
};

template <class T> class DB_RKBS :
    virtual public Integrals_RKBS<T>
{
public:
    typedef UniqueID::IDtype IDType;
    typedef Integrals_RKBS<T> Base;
    typedef typename Base::SMat_ref SMat_ref;
    typedef typename Base::SMat SMat;
   
    virtual SMat_ref Overlap() const;
    virtual SMat_ref Nuclear(const Cluster*) const;
    virtual SMat_ref RestMass() const;
private:
    virtual IDType GetID() const=0;
    virtual SMat MakeOverlap() const;
    virtual SMat MakeNuclear(const Cluster*) const;
    virtual SMat MakeRestMass() const;
    SMat MakeIntegrals(qchem::IType,const Cluster* cl=0) const;
    virtual T    Integral(qchem::IType,double ea , double eb,size_t l) const=0;

    typedef std::tuple<qchem::IType2C,IDType> id2c_t;
    mutable std::map<id2c_t ,SMat> itsBuffer; 
};

#endif
