// File: EigenSolver.H  General eigen solver.
#ifndef _EigenSolverLapack_H_
#define _EigenSolverLapack_H_

#include "Misc/EigenSolver.H"
#include "oml/numeric/EigenSolver.H"
#include "oml/numeric/SVDSolver.H"


template <class T> class EigenSolverLapackCommon
    : public virtual  EigenSolver<T>
    , protected EigenSolverCommon<T>
{
    typedef EigenSolver<T> Base;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
    typedef typename Base::UdType UdType;
public:
    EigenSolverLapackCommon();
    virtual UdType Solve(const SMat&) const;
    
    oml::EigenSolver<T>* itsLapackEigenSolver;
    oml::  SVDSolver<T>* itsLapackSVDSolver;
    double eps;
    using EigenSolverCommon<T>::V;
    using EigenSolverCommon<T>::Vd;
};

template <class T> class EigenSolverLapackEigen
    : public virtual  EigenSolver<T>
    , private EigenSolverLapackCommon<T>
{
    typedef EigenSolver<T> Base;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
public:
    EigenSolverLapackEigen(const SMat& S, double tolerance);
    using  EigenSolverLapackCommon<T>::itsLapackEigenSolver;
    using  EigenSolverLapackCommon<T>::eps;

};

template <class T> class EigenSolverLapackSVD
    : public virtual  EigenSolver<T>
    , private EigenSolverLapackCommon<T>
{
    typedef EigenSolver<T> Base;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
    typedef typename Base::RVec RVec;
public:
    EigenSolverLapackSVD(const SMat& S, double tolerance);
    using  EigenSolverLapackCommon<T>::itsLapackSVDSolver;
    using  EigenSolverLapackCommon<T>::eps;
};

template <class T> class EigenSolverLapackCholsky
    : public virtual  EigenSolver<T>
    , private EigenSolverLapackCommon<T>
{
    typedef EigenSolver<T> Base;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
public:
    EigenSolverLapackCholsky(const SMat& S, double tolerance);
};


#endif //_EigenSolver_H_
