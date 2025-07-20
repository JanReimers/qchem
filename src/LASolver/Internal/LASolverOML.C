// File: LASolverOML.C  Use OML/Numrical recipes linrary for all linear algebra ops.
export module qchem.LASolver.Internal.OML;
import qchem.LASolver.Internal.Common;

export template <class T> class LASolverOMLCommon
    : public virtual  LASolver<T>
    , protected LASolverCommon<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
    typedef typename Base::UdType UdType;
    typedef typename Base::UUdType UUdType;
public:
    LASolverOMLCommon(const LAParams& lap) : LASolverCommon<T>(lap) {};

    virtual UdType Solve(const SMat&) const;
    //! returns U,U',E  where U' has not been back transformed, U=V*Uprime.
    virtual UUdType SolveOrtho(const SMat&) const;
    
    using LASolverCommon<T>::MakeSymmetric;
    using LASolverCommon<T>::V;
    using LASolverCommon<T>::Vd;
};

export template <class T> class LASolverOMLEigen
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

export template <class T> class LASolverOMLSVD
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

export template <class T> class LASolverOMLCholsky
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
