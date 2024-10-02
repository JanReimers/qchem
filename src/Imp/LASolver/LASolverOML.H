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
    LASolverOMLCommon(const LinearAlgebraParams& lap) : LASolverCommon<T>(lap) {};
    virtual UdType Solve(const SMat&) const;
    
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
public:
    LASolverOMLEigen(const LinearAlgebraParams& lap) : LASolverOMLCommon<T>(lap) {};

    virtual void   SetBasisOverlap(const SMat& S);

    using  LASolverOMLCommon<T>::itsParams;
};

template <class T> class LASolverOMLSVD
    : public virtual  LASolver<T>
    , private LASolverOMLCommon<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
public:
    LASolverOMLSVD(const LinearAlgebraParams& lap) : LASolverOMLCommon<T>(lap) {};

    virtual void   SetBasisOverlap(const SMat& S);

    using  LASolverCommon<T>::itsParams;
};

template <class T> class LASolverOMLCholsky
    : public virtual  LASolver<T>
    , private LASolverOMLCommon<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
public:
    LASolverOMLCholsky(const LinearAlgebraParams& lap) : LASolverOMLCommon<T>(lap) {};

    virtual void   SetBasisOverlap(const SMat& S);

    using  LASolverCommon<T>::itsParams;
};


#endif //_LASolverOML_H_
