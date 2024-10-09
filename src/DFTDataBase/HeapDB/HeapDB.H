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

    virtual const RVec&    GetNormalization(bs_t& a       );
    virtual const RVec&    GetCharge       (bs_t&  a      );
    virtual const Mat&     GetOverlap    (bs_t&, bs_t&          );
    virtual const SMat&    GetOverlap    (bs_t&   ); //Numerical
    virtual const Mat&     GetRepulsion  (bs_t& ,bs_t& b);
    
    virtual const SMat&    GetOverlap    (iec_t* a); //Analytic
    virtual const SMat&    GetRepulsion  (iec_t* a   ) ; //Analytic
    virtual const Mat&     GetRepulsion  (iec_t* a,iec_t* b);
    
    virtual const ERI3&    GetOverlap3C  (iec_t* ab,iec_t* c     );
    virtual const ERI3&    GetRepulsion3C(iec_t* ab,iec_t* c);
    virtual       ERI4&    GetRepulsion4C(iecv_t&);
    virtual       ERI4&    GetExchange4C (iecv_t&);

    virtual const RVec&    GetCharge       (iec_t* a      );

    virtual const SMat&    GetKinetic         (iec_t* a     );
    virtual const SMat&    GetNuclear         (iec_t* a,const Cluster&);
    virtual const SMat&    GetInverseOverlap  (iec_t* a      );
    virtual const SMat&    GetInverseRepulsion(iec_t* a     );
//
//  These are used for charge and Vxc fitting.  They change with iterations
//  So they MUST not be cached.
//
    virtual const Vec      GetOverlap    (Vf&,Rf&);
    virtual const Vec      GetRepulsion  (Vf&,Rf&);


    virtual const NumericalIE<T>*  GetIntegralEngine(              ) const;

    virtual std::ostream&  Write(std::ostream&) const;
    virtual std::istream&  Read (std::istream&)      ;
    virtual IntegralDataBase<T>* Clone() const;

private:
    using StreamableObject::Binary;
    using StreamableObject::Ascii;
    using StreamableObject::Pretty;

    void WipeCleanAllData   ()      ;
    void BuildERIs(iecv_t& abcd);

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
