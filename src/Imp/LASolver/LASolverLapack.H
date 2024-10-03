// File: LASolver.H  General eigen solver.
#ifndef _LASolverLapack_H_
#define _LASolverLapack_H_

#include "LASolver/LASolver.H"
#include "Imp/LASolver/LASolverImp.H"
#include "oml/numeric/EigenSolver.H"
#include "oml/numeric/SVDSolver.H"


template <class T> class LASolverLapackCommon
    : public virtual  LASolver<T>
    , protected LASolverCommon<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
    typedef typename Base::UdType UdType;
public:
    LASolverLapackCommon(const LinearAlgebraParams&);

    virtual UdType Solve          (const SMat& H) const;
    
    using  LASolverCommon<T>::MakeSymmetric;
    using  LASolverCommon<T>::itsParams;

    oml::EigenSolver<T>* itsLapackEigenSolver;
    oml::  SVDSolver<T>* itsLapackSVDSolver;
    using LASolverCommon<T>::V;
    using LASolverCommon<T>::Vd;
};

template <class T> class LASolverLapackEigen
    : public virtual  LASolver<T>
    , private LASolverLapackCommon<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
    typedef typename Base::RVec RVec;
    typedef typename Base::RSMat RSMat;
public:
    LASolverLapackEigen(const LinearAlgebraParams& lap) : LASolverLapackCommon<T>(lap) {};

    virtual void   SetBasisOverlap(const SMat& S);
    virtual RSMat   Inverse(const RSMat& S) const;

    using  LASolverLapackCommon<T>::itsLapackEigenSolver;
    using  LASolverCommon<T>::MakeSymmetric;
    using  LASolverCommon<T>::itsParams;

};

template <class T> class LASolverLapackSVD
    : public virtual  LASolver<T>
    , private LASolverLapackCommon<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
    typedef typename Base::RVec RVec;
    typedef typename Base::RSMat RSMat;
public:
    LASolverLapackSVD(const LinearAlgebraParams& lap) : LASolverLapackCommon<T>(lap) {};

    virtual void   SetBasisOverlap(const SMat& S);
    virtual RSMat   Inverse(const RSMat& S) const;

    using  LASolverLapackCommon<T>::itsLapackSVDSolver;
    using  LASolverCommon<T>::MakeSymmetric;
    using  LASolverCommon<T>::itsParams;
};

template <class T> class LASolverLapackCholsky
    : public virtual  LASolver<T>
    , private LASolverLapackCommon<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
    typedef typename Base::RSMat RSMat;
public:
    LASolverLapackCholsky(const LinearAlgebraParams& lap) : LASolverLapackCommon<T>(lap) {};
    virtual void   SetBasisOverlap(const SMat& S);
    virtual RSMat   Inverse(const RSMat& S) const;
};


#endif //_LASolver_H_
