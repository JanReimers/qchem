// File: TIrrepCommon.H
#ifndef _TIrrepCommon_H_
#define _TIrrepCommon_H_

#include <Irrep_BS.H>
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
    typedef typename VectorFunction<T>::Vec     Vec;  //Vector of scalars.
    typedef typename VectorFunction<T>::Vec3Vec Vec3Vec;//vector of 3 space vectors.
    //typedef typename VectorFunction<T>::Vec3Mat Vec3Mat;//matrix of 3 space vectors.
  
public:
    TIrrepBasisSetCommon(                              );
    TIrrepBasisSetCommon(const TIrrepBasisSetCommon&);
    ~TIrrepBasisSetCommon(                              );

    virtual void Set(const LAParams&);

    using TIrrepBasisSet<T>::GetVectorSize;

    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;

protected:
    LAParams          itsLAParams; //Numerical control of general eigen solution.
};

#include "Imp/DataBase/HeapDB.H"

template <class T> class Orbital_IBS_Common
    : public virtual TOrbital_IBS<T>
    , public  TIrrepBasisSetCommon<T> 
{
    public:
    Orbital_IBS_Common() {};
    //
    //  Make a gen/ EV solver that already has the overlap S factorized.
    //
    virtual LASolver<double>* CreateSolver() const;

};

#endif //_TIrrepCommon_H_
