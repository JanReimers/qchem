// File: TBasisSetImplementation.H
#ifndef _TBasisSetImplementation_H_
#define _TBasisSetImplementation_H_

#include "BasisSet.H"
#include "FunctionsImp/VectorFunctionBuffer.H"
#include "SCFIterator/IterationParams.H"
#include "Misc/Persistent/PerRef.H"
#include "Misc/rc_ptr.H"

//---------------------------------------------------------------------
//
//  This class implements functionality common to
//  all basis sets.  It stores a list of BasisFunction*'s.
//  A basis set owns these BasisFunction*'s and deletes them in
//  the destructor.  Hence copy operations are not allowed.
//  This class also owns the IntegralDataBase.
//
template <class T> class TIrrepBasisSetCommon
    : public virtual TIrrepBasisSet<T>
    , private VectorFunctionBuffer<T>
{
protected:
    typedef VectorFunctionBuffer<T> Base;
    typedef typename Base::Mat     Mat;  //Matrix.
    typedef typename IrrepBasisSet::SMat    SMat; //Symmetrix matrix.
    typedef typename Base::Vec     Vec;  //Vector of scalars.
    typedef typename Base::Vec3    Vec3;   //3 vector (possibly complex).
    typedef typename Base::Vec3Vec Vec3Vec;//vector of 3 space vectors.
    typedef typename Base::RVec3   RVec3;  //Real space vector.
    typedef typename Base::RVec    RVec;
    typedef typename Base::Vec3Mat Vec3Mat;//matrix of 3 space vectors.
    typedef typename TIrrepBasisSet<T>::ERI3 ERI3;
  
public:
    TIrrepBasisSetCommon(                              );
    TIrrepBasisSetCommon(const LinearAlgebraParams& lap, IntegralDataBase<T>*);
    TIrrepBasisSetCommon(const TIrrepBasisSetCommon&);
    ~TIrrepBasisSetCommon(                              );

    using TIrrepBasisSet<T>::GetVectorSize;
//
//  This is where the real work gets done.
//
    virtual OrbitalGroup* CreateOrbitals(const rc_ptr<const IrrepBasisSet>&,const Hamiltonian*, const Spin&) const;
//
//  Used by fitting routines.
// 
    virtual RVec GetCharge          () const;
    virtual SMat GetInverseRepulsion() const;
    virtual SMat GetInverseOverlap  () const;
    virtual  Mat GetRepulsion  (                 const IrrepBasisSet* ff) const;
    virtual  Mat GetOverlap    (const Mesh*,     const IrrepBasisSet* ff) const;
    virtual RVec GetOverlap3C  (const SMat& Dcd, const IrrepBasisSet* ff) const;
    virtual RVec GetRepulsion3C(const SMat& Dcd, const IrrepBasisSet* ff) const;
//
//  Used by HamiltonianTerm::CalculateHamiltonianMatrix
//
    virtual SMat GetOverlap         () const;
    virtual SMat GetKinetic  (                     ) const;
    virtual SMat GetNuclear  (const Cluster*       ) const;
    virtual SMat GetOverlap  (const FittedFunction*) const;
    virtual SMat GetRepulsion(const FittedFunction*) const;
    virtual SMat GetRepulsion(const SMat& Dcd, const TIrrepBasisSet<T>* bs_cd) const;
    virtual SMat GetExchange (const SMat& Dcd, const TIrrepBasisSet<T>* bs_cd) const;
//  Load overlap (or repulsion) of this basis set with a scalar
//  funciton into a fitted function.
//
    virtual void SetFitOverlap  (FittedFunction* ff,const ScalarFunction<double>&) const;
    virtual void SetFitRepulsion(FittedFunction* ff,const ScalarFunction<double>&) const;


    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;


    virtual LASolver<T>*         GetLASolver () const;

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

protected:
    void  Insert(AnalyticIE<T>*);
private:
    IntegralDataBase<T>* GetDataBase () const;
   
    virtual void Eval    (const Mesh&,     Mat&) const;
    virtual void EvalGrad(const Mesh&, Vec3Mat&) const;
protected:
    LinearAlgebraParams          itsLAParams;
private:
    mutable IntegralDataBase<T>*  itsDataBase;
    mutable LASolver<T>*                 itsLASolver;
};

#endif //_TBasisSetImplementation_H_
