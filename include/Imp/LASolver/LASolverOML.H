// File: LASolver.H  General eigen solver.
#ifndef _LASolverOML_H_
#define _LASolverOML_H_

#include "Imp/LASolver/LASolverImp.H"

template <class T> class LASolverOMLCommon
    : public virtual  LASolver<T>
    , protected LASolverCommon<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
    typedef typename Base::UdType UdType;
public:
    LASolverOMLCommon(const LAParams& lap) : LASolverCommon<T>(lap) {};

    virtual UdType Solve(const SMat&) const;
    
    using LASolverCommon<T>::MakeSymmetric;
    using LASolverCommon<T>::V;
    using LASolverCommon<T>::Vd;
};

template <class T> class LASolverOMLEigen
    : public virtual  LASolver<T>
    , private LASolverOMLCommon<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
    typedef typename Base::RVec RVec;
    typedef typename Base::RSMat RSMat;
public:
    LASolverOMLEigen(const LAParams& lap) : LASolverOMLCommon<T>(lap) {};

    virtual void   SetBasisOverlap(const SMat& S);
    virtual RSMat   Inverse(const RSMat& S) const;

    using LASolverCommon<T>::MakeSymmetric;
    using  LASolverOMLCommon<T>::itsParams;
};

template <class T> class LASolverOMLSVD
    : public virtual  LASolver<T>
    , private LASolverOMLCommon<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
    typedef typename Base::RVec RVec;
    typedef typename Base::RSMat RSMat;
public:
    LASolverOMLSVD(const LAParams& lap) : LASolverOMLCommon<T>(lap) {};

    virtual void   SetBasisOverlap(const SMat& S);
    virtual RSMat   Inverse(const RSMat& S) const;

    using LASolverCommon<T>::MakeSymmetric;
    using  LASolverCommon<T>::itsParams;
};

template <class T> class LASolverOMLCholsky
    : public virtual  LASolver<T>
    , private LASolverOMLCommon<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
    typedef typename Base::RSMat RSMat;
public:
    LASolverOMLCholsky(const LAParams& lap) : LASolverOMLCommon<T>(lap) {};

    virtual void   SetBasisOverlap(const SMat& S);
    virtual RSMat   Inverse(const RSMat& S) const;

    using  LASolverCommon<T>::itsParams;
};


#endif //_LASolverOML_H_
