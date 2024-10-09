// File: IntegralDataBase.H  Interface for a database of integrals over basis functions.
#ifndef _IntegralDataBase_H_
#define _IntegralDataBase_H_



#include "types.H"
#include "Misc/Persistent/Persistent.H"
#include "Misc/ERI4.H"

namespace qchem
{
    // Integral types.
    enum IType3C {Overlap3C, Repulsion3C}; // <ab|c> and <ar|1/r12|c>
    enum IType2C {Overlap2C, Repulsion2C,Kinetic,Nuclear, InvOverlap, InvRepulsion, Charge, Normalization};
}
//-----------------------------------------------------------------------------
//
//  Once supplied with an integral engine, this class manages storage
//  of all matrix elements.  All calculations are done by the integral
//  engine.  A smart implementation will request a calculation of the
//  data does not yet exist in the data base.
//
template <class T> class IntegralDataBase
    : public virtual Persistent
{
public:
    typedef Matrix<T>      Mat;
    typedef SMatrix<T>     SMat;
    typedef Vector<T>      Vec;
    typedef Vector<double> RVec;
    typedef std::vector<SMat> ERI3;
    
    typedef const ScalarFunction<double> Rf;
    typedef const ScalarFunction<T>      Sf;
    typedef const VectorFunction<T>      Vf;
    
    typedef const TIrrepBasisSet<T> bs_t;
    typedef const IrrepIEClient iec_t;
    typedef std::vector<iec_t*> iecv_t;

    virtual ~IntegralDataBase() {};
    virtual void  Insert(const NumericalIE<T>*)=0;
    virtual void  Insert(const AnalyticIE<T>*)=0;
    virtual void  Insert(const BasisGroup*)=0;
    virtual void  Insert(const ERI4& J,const ERI4& K)=0;

    virtual bool  operator==(const IntegralDataBase&) const=0;
//
//  Integrals that get calculated once and then cached.
//
    virtual const RVec&    GetNormalization(bs_t& a         )=0;
    virtual const RVec&    GetCharge       (bs_t&  a      )=0;
    virtual const SMat&    GetOverlap    (bs_t& a  )=0; //Numerical
    virtual const Mat&     GetOverlap    (bs_t& a,bs_t& b) =0; //Numerical only.

    virtual const SMat&    GetOverlap    (iec_t* a)=0; //Analytic
    virtual const SMat&    GetRepulsion  (iec_t* a   ) =0;
    virtual const Mat&     GetRepulsion  (iec_t* a,iec_t* b) =0; //Analytic only.
    
    virtual const ERI3&    GetOverlap3C  (iec_t* ab,iec_t* c   ) =0; //<ab|c>
    virtual const ERI3&    GetRepulsion3C(iec_t* ab,iec_t* c  ) =0; //<a(1)b(1)|1/r12|c(2)>

    virtual       ERI4&    GetRepulsion4C(iecv_t&) =0;
    virtual       ERI4&    GetExchange4C (iecv_t&) =0;

    virtual const RVec&    GetCharge       (iec_t* a      )=0;
    
    virtual const SMat&    GetKinetic         (iec_t* a      )=0;
    virtual const SMat&    GetNuclear         (iec_t* a,const Cluster&)=0;
    virtual const SMat&    GetInverseOverlap  (iec_t* a     )=0;
    virtual const SMat&    GetInverseRepulsion(iec_t* a      )=0;
//
//  These are used for charge and Vxc fitting.  They change with iterations
//  So they MUST not be cached.
//
    virtual const Vec      GetOverlap    (Vf&,Rf&) =0;
    virtual const Vec      GetRepulsion  (Vf&,Rf&) =0;

    
    virtual const NumericalIE<T>* GetIntegralEngine(               ) const=0;

    virtual IntegralDataBase*      Clone            (               ) const=0;
    static  IntegralDataBase*      Factory          (std::istream&       )       ;
};

#endif //_IntegralDataBase_H_
