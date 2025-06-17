// File: TIBS_Common.H
#ifndef _TIBS_Common_H_
#define _TIBS_Common_H_

#include <BasisSet/Irrep_BS.H>
#include <LASolver/LAParams.H>

template <class T> class TIBS_Common
    : public virtual TIrrepBasisSet<T>
{
protected:
    typedef IrrepBasisSet Base;
    typedef typename VectorFunction<T>::Vec     Vec;  //Vector of scalars.
    typedef typename VectorFunction<T>::Vec3Vec Vec3Vec;//vector of 3 space vectors.
    //typedef typename VectorFunction<T>::Vec3Mat Vec3Mat;//matrix of 3 space vectors.
  
public:
    TIBS_Common(                              );
    TIBS_Common(const TIBS_Common&);
    ~TIBS_Common(                              );

    virtual void Set(const LAParams&);

    using TIrrepBasisSet<T>::GetVectorSize;

    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;

protected:
    LAParams          itsLAParams; //Numerical control of general eigen solution.
};

template <class T> class Orbital_IBS_Common
    : public virtual TOrbital_IBS<T>
    , public  TIBS_Common<T> 
{
    public:
    Orbital_IBS_Common() {};
    //
    //  Make a gen/ EV solver that already has the overlap S factorized.
    //
    virtual LASolver<double>* CreateSolver() const;

};

#endif //_TIBS_Common_H_
