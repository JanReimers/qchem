// File: BasisSet.H  Abstract interface for a basis set.
#ifndef _BasisSet_H_
#define _BasisSet_H_

#include "types.H"
#include "Misc/Persistent/Persistent.H"



//--------------------------------------------------------------
//
//  This interface of the basis function provides no notion
//  of weather the function is real or complex (or vector!)
//  valued.  Normalization and total charge should always
//  be real.
//
class BasisFunction
    : public Persistent
{
public:
    virtual ~BasisFunction()  {};

    virtual bool   operator==(const BasisFunction&) const=0;

    virtual BasisFunction* Clone  (             ) const=0;
    static  BasisFunction* Factory(std::istream&)        ;
};

#include "Functions/ScalarFunction.H"
//---------------------------------------------------------------
//
//  Here the basis becomes a real or complex valued function
//  depending of T.
//
template <class T> class TBasisFunction
    : public virtual BasisFunction
    , public virtual ScalarFunction<T>
{
};


#include <vector>
//----------------------------------------------------------------------------
//
//  Interface for an irreducible representation basis sets.  H is block diagonal with one
//  block for  IrrepBasisSet,  For atoms each L get.s an IrrepBasisSet and an H  block. 
//  All matrices of integrals are supplied by the
//  IntegralDataBase, which only calculates the integrals once and buffers them.
//  The Factory allows one to read basis sets from a stream polymorphically.
//  Clone(RVec&) creates a similar set with a new center.
//  The quantum number could be L for atoms, Irreducable rep for molecules, or
//  the wave vector k for solids.
//  Persistent is derived from StreamableObject which has a virtual destructor
//
#include "IEClient.H"

class IrrepBasisSet
    : public virtual IrrepIEClient
    , public UniqueID
{
public:
    typedef std::vector<BasisFunction*>::      iterator       iterator;
    typedef std::vector<BasisFunction*>::const_iterator const_iterator;

    typedef Vector <double> RVec;
    typedef SMatrix<double> SMat;
    typedef  Matrix<double>  Mat;
    typedef std::vector<SMat> ERI3;
    
    virtual bool                 operator==      (const IrrepBasisSet&) const=0;
    virtual int                  GetStartIndex   (               ) const=0;
    virtual size_t               GetNumFunctions (               ) const=0;
    virtual const QuantumNumber& GetQuantumNumber(               ) const=0;
//
//  This is where the real work gets done.
//
    virtual OrbitalGroup* CreateOrbitals(const rc_ptr<const IrrepBasisSet>&,const Hamiltonian*, const Spin&) const=0;
//
//  Used by fitting routines.
//
    virtual RVec GetCharge          () const=0;
    virtual SMat GetInverseRepulsion() const=0;
    virtual SMat GetInverseOverlap  () const=0;
    virtual  Mat GetOverlap    (const Mesh*,     const IrrepBasisSet* ff) const=0;
    virtual  Mat GetRepulsion  (                 const IrrepBasisSet* ff) const=0;
    virtual RVec GetOverlap    (const Mesh*,     const ScalarFunction<double>* sf) const=0;
    virtual RVec GetRepulsion  (const Mesh*,     const ScalarFunction<double>* sf) const=0;
    virtual RVec GetOverlap3C  (const SMat& Dcd, const IrrepBasisSet* ff) const=0;
    virtual RVec GetRepulsion3C(const SMat& Dcd, const IrrepBasisSet* ff) const=0;
    virtual const ERI3& GetOverlap3C  (const IrrepBasisSet* ff) const=0;
    virtual const ERI3& GetRepulsion3C  (const IrrepBasisSet* ff) const=0;
//
//  Used by HamiltonianTerm::CalculateHamiltonianMatrix
//
    virtual SMat GetOverlap  (                     ) const=0;
    virtual SMat GetKinetic  (                     ) const=0;
    virtual SMat GetNuclear  (const Cluster*       ) const=0;
     
    virtual const_iterator begin() const=0;
    virtual const_iterator end  () const=0;
    virtual       iterator begin()      =0;
    virtual       iterator end  ()      =0;
    
    
    virtual UniqueID::IDtype GetID() const {return UniqueID::GetID();}
//
//  StreamableObject stuff.
//
    virtual std::ostream&  Write(std::ostream&    ) const=0;
    virtual std::istream&  Read (std::istream&    )      =0;
    virtual IrrepBasisSet* Clone  (            ) const=0; //TODO do we need this?  It prevent usage of IDs.
    virtual IrrepBasisSet* Clone  (const RVec3&) const=0;
    static  IrrepBasisSet* Factory(std::istream&    )        ;
private:
    friend class BasisSetImp;
    virtual void SetStartIndex(int)=0;
};

#include "Functions/VectorFunction.H"
#include "Imp/Containers/ptr_vector.h"
//----------------------------------------------------------------------------
//
//  Extend basis to be a set of real or complex valued functions
//
template <class T> class TIrrepBasisSet
    : public virtual IrrepBasisSet
    , public virtual VectorFunction<T>
{
public:
    typedef dynamic_cast_iterator<BasisFunction*,TBasisFunction<T>* > const_iterator;
    typedef std::vector<SMat> ERI3;
    size_t GetVectorSize() const {return GetNumFunctions();}
    virtual  Mat GetRepulsion(const IrrepBasisSet* ff) const=0; //unhide
    virtual SMat GetRepulsion(const SMat& Dcd, const TIrrepBasisSet* bs_cd) const=0;
    virtual SMat GetExchange (const SMat& Dcd, const TIrrepBasisSet* bs_cd) const=0;

    virtual IntegralDataBase<T>* GetDataBase() const=0;
    
    const_iterator beginT() const {return const_iterator(begin());};
};



//----------------------------------------------------------------------------
//
//  Interface for all basis sets.  All matrices of integrals are supplied by the
//  IntegralDataBase, which only calculates the integrals once and buffers them.
//  The Factory allows one to read basis sets from a stream polymorphically.
//  Clone(RVec&) creates a similar set with a new center.
//  The quantum number could be L for atoms, Irreducable rep for molecules, or
//  the wave vector k for solids.
//  Persistent is derived from StreamableObject which has a virtual destructor
//
class BasisSet
    : public virtual IEClient
    , public UniqueID
{
public:
    typedef optr_vector1<IrrepBasisSet*>::const_iterator const_iterator;

    BasisSet() {};

    virtual ~BasisSet();
    virtual size_t GetNumFunctions() const=0;
    virtual size_t GetNumBasisSets() const=0;

    virtual const_iterator begin() const=0;
    virtual const_iterator end  () const=0;
    virtual IntegralDataBase<double>* GetDataBase() const =0;

    virtual UniqueID::IDtype GetID() const {return UniqueID::GetID();}
    
    static  BasisSet*  Factory(std::istream&    )      ;
};


#endif //_BasisSet_H_
