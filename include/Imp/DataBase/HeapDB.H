// File: HeapDB.H  Implement an integral data base that uses heap for storage.
#ifndef _HeapDB_H_
#define _HeapDB_H_



#include "IntegralDataBase.H"
#include <IEClient.H>
#include <AnalyticIE.H>
#include "Imp/Containers/ERI4.H"
#include <map>

//-----------------------------------------------------------------------------
//
//  This class implements an IntegralDataBase that simply stores all integrals
//  on the heap.  This means RAM and/or vitual memory in a real OS.
//
template <class T> class HeapDB :
    virtual public IntegralDataBase<T>
{
protected:
    typedef IntegralDataBase<T> Base;
    typedef typename Base::Vec   Vec;
    typedef typename Base::Mat   Mat;
    typedef typename Base::RVec  RVec;
    typedef typename Base::SMat  SMat;
    typedef typename Base::ERI3 ERI3;
    typedef UniqueID::IDtype IDType;
    typedef typename Base::bs_t bs_t;
    typedef typename Base::iec_t iec_t;
    typedef typename Base::Rf Rf;
    typedef typename Base::Vf Vf;
public:

    HeapDB();
    virtual ~HeapDB();
    HeapDB(AnalyticIE<T>* ie);

#ifndef UT_COVERAGE_ONLY
    virtual const RVec&    GetNormalization(bs_t& a);
    virtual const SMat&    GetOverlap      (bs_t& a);
    virtual const Mat&     GetRepulsion    (bs_t& ,bs_t& b);
#endif

  //
//  These are used for charge and Vxc fitting.  They change with iterations
//  So they MUST not be cached.
//
    virtual const Vec      GetOverlap    (const Mesh*,bs_t&,Rf&   f); //Numerical .
    virtual const Vec      GetRepulsion  (const Mesh*,bs_t&,Rf&   f); //Numerical .
    virtual const RVec&    GetCharge     (const Mesh*,bs_t&        ); //Numerical .
    virtual const Mat&     GetOverlap    (const Mesh*,bs_t&,bs_t& b); //Numerical .
//
//  These are all done analytically.
//
    virtual const RVec&    GetCharge          (iec_t* a);
    //virtual const SMat&    GetOverlap         (iec_t* a); 
    // virtual const SMat&    GetRepulsion       (iec_t* a);
    // virtual const SMat&    GetInverseOverlap  (iec_t* a,const LAParams&); //uses analytic overlap
    // virtual const SMat&    GetInverseRepulsion(iec_t* a,const LAParams&); //Uses analytic overlap
    // virtual const SMat&    GetKinetic         (iec_t* a);
    // virtual const SMat&    GetNuclear         (iec_t* a,const Cluster&);
    //virtual const SMat&    GetRestMass        (iec_t* a);
    
    virtual const Mat&     GetRepulsion  (iec_t* a ,iec_t* b) ; //<a(1)|1/r12|b(2)>
    virtual const ERI3&    GetOverlap3C  (iec_t* ab,iec_t* c) ; //<ab|c>
    virtual const ERI3&    GetRepulsion3C(iec_t* ab,iec_t* c) ; //<a(1)b(1)|1/r12|c(2)>

    // virtual    ERI4      GetDirect__4C(bs_t& a,bs_t& b);
    // virtual    ERI4      GetExchange4C(bs_t& a,bs_t& b);

    virtual void Report(std::ostream&) const;

private:
    const RVec& GetNumericalNormalization(const Mesh* m,bs_t& a);


    AnalyticIE<T>*     itsAnalyticIE; //Integral engine for analytic integrals.
    
    typedef std::tuple<qchem::IType2C,IDType> id2c_t;
    typedef std::tuple<qchem::IType2C,IDType,IDType> id2cx_t; //Cross basis set integrals.
    typedef std::tuple<qchem::IType3C,IDType,IDType> id3c_t; //Cross basis set integrals.
    std::map<id2c_t ,RVec> its1C; //Norm, Charge 
    std::map<id2c_t ,Vec> its1Cx; //Norm, Charge Overlap(function)
    std::map<id2c_t ,SMat> its2C; //Overlap, Rep, Kinetic, InvOverlap, InvRepulsion
    std::map<id2cx_t,SMat> its2CNuc; //Nuclear IDed by cluster
    std::map<id2cx_t, Mat> its2Cx; //Overlap, Rep, Kinetic, InvOverlap, InvRepulsion
    std::map<id3c_t ,ERI3> its3C; //Overlap, Rep, Kinetic, InvOverlap, InvRepulsion
    // AnalyticIE<double>::erij_t Jac; 
    // AnalyticIE<double>::erij_t Kab;

};

template <class T> class DB_1E :
    virtual public Integrals_1E<T>
{    
public:
    typedef Integrals_1E<T> Base;
    typedef typename Base::SMat_ref SMat_ref;
    //  This is all buffered.
    virtual SMat_ref Overlap() const;
    virtual SMat_ref Kinetic() const;
    virtual SMat_ref Nuclear(const Cluster*) const;
protected:
    typedef SMatrix<T> SMat;
private:
    typedef UniqueID::IDtype IDType;
    virtual IDType GetID() const=0;
    // One time calls to un-buffered integral calculations.
    virtual SMat MakeOverlap() const=0;
    virtual SMat MakeKinetic() const=0;
    virtual SMat MakeNuclear(const Cluster*) const=0;

    typedef std::tuple<qchem::IType2C,IDType> id2c_t;
    mutable std::map<id2c_t ,SMat> itsBuffer; 
};

template <class T> class DB_BS_2E : public virtual Integrals_BS_2E<T>
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
    typedef std::map<IDType,std::map<IDType,ERI4> > erij_t;
    mutable erij_t Jac,Kab;
};

template <class T> class DB_2E :  virtual public Integrals_2E<T>
{
public:
    typedef UniqueID::IDtype IDType;
    typedef Integrals_2E<T> Base;
    typedef typename Base::bs_t bs_t;
    virtual ERI4 Direct  (const bs_t& c) const;
    virtual ERI4 Exchange(const bs_t& b) const;
protected:
    DB_2E() : itsDB_BS_2E(0) {};
    DB_2E(const DB_BS_2E<T>* db);
private:
    virtual IDType GetID() const=0;
    const DB_BS_2E<T>* itsDB_BS_2E; //Database of all supermatrix tables.
};

class DB_Fit : public virtual FitIntegrals_1E<double>
{
public:
    typedef FitIntegrals_1E<double> Base;
    typedef typename Base::SMat_ref SMat_ref;
    virtual SMat_ref Overlap() const;
    virtual SMat_ref Repulsion() const;
    virtual SMat_ref InvOverlap(const LAParams&) const;
    virtual SMat_ref InvRepulsion(const LAParams&) const;
    protected:
    typedef SMatrix<double> SMat;
private:
    typedef UniqueID::IDtype IDType;
    virtual IDType GetID() const=0;
    // One time calls to un-buffered integral calculations.
    virtual SMat MakeOverlap() const=0;
    virtual SMat MakeRepulsion() const=0;
    //! \brief Return the Penrose inverse of a symmetric matrix using SVD decomposition
    //! If \f$ S=UsV^{\dagger} \f$, then \f$ S^{-1}=V\frac{1}{s}U^{\dagger} \f$
    static  SMat MakeInverse  (SMat_ref,const LAParams&); //Numerically stable algo required.

    typedef std::tuple<qchem::IType2C,IDType> id2c_t;
    mutable std::map<id2c_t ,SMat> itsBuffer; 
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
