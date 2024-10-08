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
public:

    HeapDB();
    virtual void  Insert(const TBasisSet<T>*,const NumericalIE<T>*);
    virtual void  Insert(const TBasisSet<T>*,const AnalyticIE<T>*);
    virtual void  Insert(const BasisGroup*);
    virtual void  Insert(const ERI4& J,const ERI4& K);
    virtual bool  operator==(const IntegralDataBase<T>&) const;

    virtual const SMat&    GetOverlap    (                             );
    virtual const Mat&     GetOverlap    (const TBasisSet<T>&          );
    virtual const Vec      GetOverlap    (const ScalarFunction<double>&);
    virtual const ERI3&    GetOverlap3C  (const TBasisSet<T>&          );

    virtual const SMat&    GetRepulsion  (                             );
    virtual const Mat&     GetRepulsion  (const TBasisSet<T>&          );
    virtual const Vec      GetRepulsion  (const ScalarFunction<double>&);
    virtual const ERI3&    GetRepulsion3C(const TBasisSet<T>&          );
    virtual       ERI4view GetRepulsion4C(const TBasisSet<T>* bs_cd    );
    virtual       ERI4view GetExchange4C (const TBasisSet<T>* bs_cd    );

    virtual const RVec&    GetNormalization (                );
    virtual const RVec&    GetCharge        (                );

    virtual const SMat&    GetKinetic         (              );
    virtual const SMat&    GetNuclear         (const Cluster&);
    virtual const SMat&    GetInverseOverlap  (              );
    virtual const SMat&    GetInverseRepulsion(              );

    virtual const NumericalIE<T>*  GetIntegralEngine(              ) const;

    virtual std::ostream&  Write(std::ostream&) const;
    virtual std::istream&  Read (std::istream&)      ;
    virtual IntegralDataBase<T>* Clone() const;

private:
    using StreamableObject::Binary;
    using StreamableObject::Ascii;
    using StreamableObject::Pretty;

    void WipeCleanAllData   ()      ;
    void BuildERIs();

    typedef std::vector<UniqueID::IDtype>      IDList;

    const NumericalIE<T>* itsNumericalIE;
    const AnalyticIE<T>*     itsAnalyticIE; //Integral engine for analytic integrals.
    const IrrepBasisSet*          itsBasisSet;
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
