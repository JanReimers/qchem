// File: TIrrepCommon.H
#ifndef _TIrrepCommon_H_
#define _TIrrepCommon_H_

#include <BasisSet.H>
#include <LAParams.H>
//---------------------------------------------------------------------
//
//  This class also owns the IntegralDataBase.
//
template <class T> class TIrrepBasisSetCommon
    : public virtual TIrrepBasisSet<T>
{
protected:
    typedef IrrepBasisSet Base;
    typedef typename Base::RVec    RVec;
    typedef typename Base::Mat     Mat;  //Matrix.
    typedef typename Base::SMat    SMat; //Symmetrix matrix.
    typedef typename Base::ERI3 ERI3;
    typedef typename VectorFunction<T>::Vec     Vec;  //Vector of scalars.
    typedef typename VectorFunction<T>::Vec3Vec Vec3Vec;//vector of 3 space vectors.
    typedef typename VectorFunction<T>::Vec3Mat Vec3Mat;//matrix of 3 space vectors.
  
public:
    TIrrepBasisSetCommon(                              );
    TIrrepBasisSetCommon(const LAParams& lap, IntegralDataBase<T>*);
    TIrrepBasisSetCommon(const TIrrepBasisSetCommon&);
    ~TIrrepBasisSetCommon(                              );

    using TIrrepBasisSet<T>::GetVectorSize;
//
//  Make a gen/ EV solver that already has the overlap S factorized.
//
//    virtual LASolver<double>* CreateSolver() const;
//
//  Used by fitting routines.
// 
    virtual RVec GetCharge          () const;
    virtual SMat GetInverseRepulsion(const LAParams&) const;
    virtual SMat GetInverseOverlap  (const LAParams&) const;
    virtual  Mat GetRepulsion  (                 const IrrepBasisSet* ff) const;
    virtual  Mat GetOverlap    (const Mesh*,     const IrrepBasisSet* ff) const;
    virtual RVec GetOverlap    (const Mesh*,     const ScalarFunction<double>* sf) const;
    virtual RVec GetRepulsion  (const Mesh*,     const ScalarFunction<double>* sf) const;
    virtual RVec GetOverlap3C  (const SMat& Dcd, const IrrepBasisSet* ff) const;
    virtual RVec GetRepulsion3C(const SMat& Dcd, const IrrepBasisSet* ff) const;
    virtual const ERI3& GetOverlap3C  (const IrrepBasisSet* ff) const;
    virtual const ERI3& GetRepulsion3C  (const IrrepBasisSet* ff) const;
//
//  Used by HamiltonianTerm::CalculateHamiltonianMatrix
//
    virtual SMat GetOverlap () const;
    virtual SMat GetKinetic (              ) const;
    virtual SMat GetNuclear (const Cluster*) const;
    virtual SMat GetRestMass(              ) const;
    virtual SMat GetRepulsion(const SMat& Dcd, const TIrrepBasisSet<T>* bs_cd) const;
    virtual SMat GetExchange (const SMat& Dcd, const TIrrepBasisSet<T>* bs_cd) const;

    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

protected:
    IntegralDataBase<T>* GetDataBase () const;
   
//    virtual void EvalGrad(const Mesh&, Vec3Mat&) const;
protected:
    LAParams          itsLAParams; //Numerical control of general eigen solution.
private:
    mutable IntegralDataBase<T>*  itsDataBase;
};

#include "Imp/DataBase/HeapDB.H"

template <class T> class Orbital_IBS_Common
    : public virtual Orbital_IBS<T>
    , protected  TIrrepBasisSetCommon<T> 
{
    public:
    Orbital_IBS_Common() {};
    Orbital_IBS_Common(const LAParams& lap, IntegralDataBase<T>* db) 
        : TIrrepBasisSetCommon<T>(lap,db) {};

    //
    //  Make a gen/ EV solver that already has the overlap S factorized.
    //
    virtual LASolver<double>* CreateSolver() const;

    using TIrrepBasisSetCommon<T>::GetDataBase;
};

#endif //_TIrrepCommon_H_
