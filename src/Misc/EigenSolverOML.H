// File: EigenSolver.H  General eigen solver.
#ifndef _EigenSolverOML_H_
#define _EigenSolverOML_H_

#include "Misc/EigenSolver.H"

template <class T> class EigenSolverOMLCommon
    : public virtual  EigenSolver<T>
    , protected EigenSolverCommon<T>
{
    typedef EigenSolver<T> Base;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
    typedef typename Base::UdType UdType;
public:
    virtual UdType Solve(const SMat&) const;
    
    using EigenSolverCommon<T>::V;
    using EigenSolverCommon<T>::Vd;
};

template <class T> class EigenSolverOMLEigen
    : public virtual  EigenSolver<T>
    , private EigenSolverOMLCommon<T>
{
    typedef EigenSolver<T> Base;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
public:
    EigenSolverOMLEigen(const SMat& S, double tolerance);
};

template <class T> class EigenSolverOMLSVD
    : public virtual  EigenSolver<T>
    , private EigenSolverOMLCommon<T>
{
    typedef EigenSolver<T> Base;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
public:
    EigenSolverOMLSVD(const SMat& S, double tolerance);
};

template <class T> class EigenSolverOMLCholsky
    : public virtual  EigenSolver<T>
    , private EigenSolverOMLCommon<T>
{
    typedef EigenSolver<T> Base;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
public:
    EigenSolverOMLCholsky(const SMat& S, double tolerance);
};


#endif //_EigenSolverOML_H_
