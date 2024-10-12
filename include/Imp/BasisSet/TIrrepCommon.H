// File: TBasisSetImplementation.H
#ifndef _TBasisSetImplementation_H_
#define _TBasisSetImplementation_H_

#include "BasisSet.H"
#include "Imp/LAParams.H"
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
//  This is where the real work gets done.
//
    virtual OrbitalGroup* CreateOrbitals(const Hamiltonian*, const Spin&) const;
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
    virtual SMat GetOverlap         () const;
    virtual SMat GetKinetic  (                     ) const;
    virtual SMat GetNuclear  (const Cluster*       ) const;
    virtual SMat GetRepulsion(const SMat& Dcd, const TIrrepBasisSet<T>* bs_cd) const;
    virtual SMat GetExchange (const SMat& Dcd, const TIrrepBasisSet<T>* bs_cd) const;

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
    LAParams          itsLAParams; //Numerical control of general eigen solution.
private:
    mutable IntegralDataBase<T>*  itsDataBase;
    mutable LASolver<T>*          itsLASolver;
};

#endif //_TBasisSetImplementation_H_
