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
    typedef typename VectorFunction<T>::Vec     Vec;  //Vector of scalars.
    typedef typename VectorFunction<T>::Vec3Vec Vec3Vec;//vector of 3 space vectors.
    //typedef typename VectorFunction<T>::Vec3Mat Vec3Mat;//matrix of 3 space vectors.
  
public:
    TIrrepBasisSetCommon(                              );
    TIrrepBasisSetCommon(const LAParams& lap);
    TIrrepBasisSetCommon(const TIrrepBasisSetCommon&);
    ~TIrrepBasisSetCommon(                              );

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
    Orbital_IBS_Common(const LAParams& lap) 
        : TIrrepBasisSetCommon<T>(lap) {};

    //
    //  Make a gen/ EV solver that already has the overlap S factorized.
    //
    virtual LASolver<double>* CreateSolver() const;

};

template <class T> class Orbital_DFT_IBS_Common
    : public virtual TOrbital_DFT_IBS<T>
{
    typedef typename Integrals_Base<T>::SMat SMat;
    typedef typename Integrals_Base<T>::fbs_t fbs_t;
    typedef typename Integrals_DFT<T>::Vec  Vec;
public:
    using Integrals_DFT<T>::Overlap3C; //Unhide
    using Integrals_DFT<T>::Repulsion3C; //Unhide
    virtual Vec Overlap3C  (const SMat& Dcd, const fbs_t* ff) const;
    virtual Vec Repulsion3C(const SMat& Dcd, const fbs_t* ff) const;
};

template <class T> class Orbital_HF_IBS_Common
    : public virtual TOrbital_HF_IBS<T>
{
    typedef typename Integrals_Base<T>::SMat SMat;
public:
    virtual SMat Direct  (const SMat& Dcd, const TIrrepBasisSet<T>* bs_cd) const;
    virtual SMat Exchange(const SMat& Dcd, const TIrrepBasisSet<T>* bs_cd) const;
};

class Fit_IBS_Common : public virtual Fit_IBS, public virtual FitIntegrals_1E<double>
{
    typedef Integrals_Base<double>::Vec Vec;
    typedef Integrals_Base<double>::Mat Mat;
    typedef Integrals_Base<double>::bs_t bs_t;
    public:
    virtual  Vec MakeNorm   (const Mesh*        ) const; //Numerical .
    virtual  Vec MakeCharge (const Mesh*        ) const; //Numerical .
    virtual  Mat MakeOverlap(const Mesh*,const bs_t& b) const; //Numerical X overlap.

    virtual const Vec Overlap  (const Mesh*,const Sf&) const; //Numerical  
    virtual const Vec Repulsion(const Mesh*,const Sf&) const; //Numerical 
};

#endif //_TIrrepCommon_H_
