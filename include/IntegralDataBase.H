// File: IntegralDataBase.H  Interface for a database of integrals over basis functions.
#ifndef _IntegralDataBase_H_
#define _IntegralDataBase_H_



#include "types.H"
#include "Misc/Persistent/Persistent.H"
#include "Misc/ERIProxy.H"

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
    typedef std::vector<SMat> MList;

    virtual ~IntegralDataBase() {};
    virtual void  Insert(const TBasisSet<T>*, const NumericalIE<T>*)=0;
    virtual void  Insert(const TBasisSet<T>*, const AnalyticIE<T>*)=0;
    virtual void  Insert(const BasisGroup*)=0;
    virtual void  Insert(const ERIList1& J,const ERIList1& K)=0;

    virtual bool  operator==(const IntegralDataBase&) const=0;

    virtual const SMat&    GetOverlap    (                             ) =0;
    virtual const Mat&     GetOverlap    (const TBasisSet<T>&          ) =0;
    virtual const Vec&     GetOverlap    (const ScalarFunction<double>&) =0;
    virtual const MList&   GetOverlap3C  (const TBasisSet<T>&          ) =0;

    virtual const SMat&    GetRepulsion  (                             ) =0;
    virtual const Mat&     GetRepulsion  (const TBasisSet<T>&          ) =0;
    virtual const Vec&     GetRepulsion  (const ScalarFunction<double>&) =0;
    virtual const MList&   GetRepulsion3C(const TBasisSet<T>&          ) =0;
    virtual ERIProxy1 GetRepulsion4C_1(const TBasisSet<T>* bs_cd    ) =0;
    virtual ERIProxy1 GetExchange4C_1 (const TBasisSet<T>* bs_cd    ) =0;

    virtual const RVec&    GetNormalization(              )=0;
    virtual const RVec&    GetCharge       (              )=0;

    virtual const SMat&    GetKinetic         (              )=0;
    virtual const SMat&    GetNuclear         (const Cluster&)=0;
    virtual const SMat&    GetInverseOverlap  (              )=0;
    virtual const SMat&    GetInverseRepulsion(              )=0;

    virtual const NumericalIE<T>* GetIntegralEngine(               ) const=0;

    virtual IntegralDataBase*      Clone            (               ) const=0;
    static  IntegralDataBase*      Factory          (std::istream&       )       ;
};

#endif //_IntegralDataBase_H_
