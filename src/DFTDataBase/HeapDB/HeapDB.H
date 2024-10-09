// File: HeapDB.H  Implement an integral data base that uses heap for storage.
#ifndef _HeapDB_H_
#define _HeapDB_H_



#include "IntegralDataBase.H"
#include "Misc/Persistent/IDRef.H"
#include "Misc/ERI4.H"
#include "Misc/rc_ptr.H"
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
    typedef typename Base::iecv_t iecv_t;
    typedef typename Base::Rf Rf;
    typedef typename Base::Vf Vf;
public:

    HeapDB();
    virtual void  Insert(const NumericalIE<T>*);
    virtual void  Insert(const AnalyticIE<T>*);
    virtual void  Insert(const BasisGroup*);
    virtual void  Insert(const ERI4& J,const ERI4& K);
    virtual bool  operator==(const IntegralDataBase<T>&) const;

#ifndef UT_COVERAGE_ONLY
    virtual const RVec&    GetNormalization(bs_t& a);
    virtual const SMat&    GetOverlap      (bs_t& a);
    virtual const Mat&     GetRepulsion    (bs_t& ,bs_t& b);
#endif

  //
//  These are used for charge and Vxc fitting.  They change with iterations
//  So they MUST not be cached.
//
    virtual const Vec      GetOverlap    (bs_t&,Rf&   f); //Numerical .
    virtual const Vec      GetRepulsion  (bs_t&,Rf&   f); //Numerical .
    virtual const RVec&    GetCharge     (bs_t&        ); //Numerical .
    virtual const Mat&     GetOverlap    (bs_t&,bs_t& b); //Numerical .
//
//  These are all done analytically.
//
    virtual const RVec&    GetCharge          (iec_t* a);
    virtual const SMat&    GetOverlap         (iec_t* a); 
    virtual const SMat&    GetRepulsion       (iec_t* a);
    virtual const SMat&    GetInverseOverlap  (iec_t* a); //uses analytic overlap
    virtual const SMat&    GetInverseRepulsion(iec_t* a); //Uses analytic overlap
    virtual const SMat&    GetKinetic         (iec_t* a);
    virtual const SMat&    GetNuclear         (iec_t* a,const Cluster&);

    virtual const Mat&     GetRepulsion  (iec_t* a ,iec_t* b) ; //<a(1)|1/r12|b(2)>
    virtual const ERI3&    GetOverlap3C  (iec_t* ab,iec_t* c) ; //<ab|c>
    virtual const ERI3&    GetRepulsion3C(iec_t* ab,iec_t* c) ; //<a(1)b(1)|1/r12|c(2)>

    virtual       ERI4&    GetRepulsion4C(iecv_t&) ;
    virtual       ERI4&    GetExchange4C (iecv_t&)  ;

    virtual std::ostream&  Write(std::ostream&) const;
    virtual std::istream&  Read (std::istream&)      ;
    virtual IntegralDataBase<T>* Clone() const;

private:
    using StreamableObject::Binary;
    using StreamableObject::Ascii;
    using StreamableObject::Pretty;

    void WipeCleanAllData   ()      ;
    void BuildERIs(iecv_t& abcd);
    const RVec& GetNumericalNormalization(bs_t& a);

    typedef std::vector<UniqueID::IDtype>      IDList;

    const NumericalIE<T>* itsNumericalIE;
    const AnalyticIE<T>*     itsAnalyticIE; //Integral engine for analytic integrals.
    const BasisGroup*        itsBasisGroup;
 
    typedef std::tuple<qchem::IType2C,IDType> id2c_t;
    typedef std::tuple<qchem::IType2C,IDType,IDType> id2cx_t; //Cross basis set integrals.
    typedef std::tuple<qchem::IType3C,IDType,IDType> id3c_t; //Cross basis set integrals.
    std::map<id2c_t ,RVec> its1C; //Norm, Charge 
    std::map<id2c_t ,Vec> its1Cx; //Norm, Charge Overlap(function)
    std::map<id2c_t ,SMat> its2C; //Overlap, Rep, Kinetic, InvOverlap, InvRepulsion
    std::map<id2cx_t,SMat> its2CNuc; //Nuclear IDed by cluster
    std::map<id2cx_t, Mat> its2Cx; //Overlap, Rep, Kinetic, InvOverlap, InvRepulsion
    std::map<id3c_t ,ERI3> its3C; //Overlap, Rep, Kinetic, InvOverlap, InvRepulsion
    ERI4 itsJTable; //Coulomb integrals.
    ERI4 itsKTable; //Exchange integrals.
    
};

#endif
