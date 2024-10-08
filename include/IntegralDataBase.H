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
    typedef const TIrrepBasisSet<T> bs_t;
    typedef const IrrepIEClient iec_t;

    virtual ~IntegralDataBase() {};
    virtual void  Insert(bs_t*, const NumericalIE<T>*)=0;
    virtual void  Insert(bs_t*, const AnalyticIE<T>*)=0;
    virtual void  Insert(const BasisGroup*)=0;
    virtual void  Insert(const ERI4& J,const ERI4& K)=0;

    virtual bool  operator==(const IntegralDataBase&) const=0;
//
//  Integrals that get calculated once and then cached.
//
    virtual const SMat&    GetOverlap    (bs_t&   )=0; //Numerical
    virtual const SMat&    GetOverlap    (iec_t* a)=0; //Analytic
    virtual const Mat&     GetOverlap    (bs_t& a,bs_t& b) =0; //Numerical only.
    
    virtual const ERI3&    GetOverlap3C  (bs_t&          ) =0;

    virtual const SMat&    GetRepulsion  (                             ) =0;
    //virtual const Mat&     GetRepulsion  (bs_t& ,bs_t& b ) =0;
    virtual const Mat&     GetRepulsion  (iec_t* a,iec_t* b) =0; //Analytic only.
    virtual const ERI3&    GetRepulsion3C(bs_t&          ) =0;
    virtual       ERI4view GetRepulsion4C(bs_t* bs_cd    ) =0;
    virtual       ERI4view GetExchange4C (bs_t* bs_cd    ) =0;

    virtual const RVec&    GetNormalization(              )=0;
    virtual const RVec&    GetCharge       (              )=0;

    virtual const SMat&    GetKinetic         (              )=0;
    virtual const SMat&    GetNuclear         (const Cluster&)=0;
    virtual const SMat&    GetInverseOverlap  (              )=0;
    virtual const SMat&    GetInverseRepulsion(              )=0;
//
//  These are used for charge and Vxc fitting.  They change with iterations
//  So they MUST not be cached.
//
    virtual const Vec      GetOverlap    (const ScalarFunction<double>&) =0;
    virtual const Vec      GetRepulsion  (const ScalarFunction<double>&) =0;

    
    virtual const NumericalIE<T>* GetIntegralEngine(               ) const=0;

    virtual IntegralDataBase*      Clone            (               ) const=0;
    static  IntegralDataBase*      Factory          (std::istream&       )       ;
};

#endif //_IntegralDataBase_H_
