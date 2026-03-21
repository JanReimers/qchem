// File: LASolver_blaze/Internal/LASolverLapack.C  Use blaze interface to Lapack linrary for all linear algebra ops.
export module qchem.LASolver_blaze.Internal.Lapack;
import qchem.LASolver_blaze.Internal.Common;


export template <class T> class LASolverLapackCommon_blaze
    : public virtual  LASolver_blaze<T>
    , protected LASolverCommon_blaze<T>
{
    typedef LASolver_blaze<T> Base;
    typedef typename Base:: rvec_t  rvec_t;
    typedef typename Base::  mat_t   mat_t;
    typedef typename Base:: smat_t  smat_t;
    typedef typename Base:: umat_t  umat_t;
    typedef typename Base:: lmat_t  lmat_t;
    typedef typename Base::rsmat_t rsmat_t;
    typedef typename Base:: dmat_t  dmat_t;
    typedef typename Base::   Ud_t    Ud_t;
    typedef typename Base::  UUd_t   UUd_t; //U,U',E  where U' has not been back transformed, U=V*Uprime.
public:
    LASolverLapackCommon_blaze(const LAParams&);
    ~LASolverLapackCommon_blaze();
    
    virtual  Ud_t Solve     (const smat_t& H) const;
    //! returns U,U',E  where U' has not been back transformed, U=V*Uprime.
    virtual UUd_t SolveOrtho(const smat_t& Hprime) const; //Hprime = Vd * H * V Hamiltonian/Fock matrix.
   
    using  LASolverCommon_blaze<T>::MakeSymmetric;
    using  LASolverCommon_blaze<T>::itsParams;

    using LASolverCommon_blaze<T>::V;
    using LASolverCommon_blaze<T>::Vd;
};

export template <class T> class LASolverLapackEigen_blaze
    : public virtual  LASolver_blaze<T>
    , private LASolverLapackCommon_blaze<T>
{
    typedef LASolver_blaze<T> Base;
    typedef typename Base:: rvec_t  rvec_t;
    typedef typename Base::  mat_t   mat_t;
    typedef typename Base:: smat_t  smat_t;
    typedef typename Base:: umat_t  umat_t;
    typedef typename Base:: lmat_t  lmat_t;
    typedef typename Base::rsmat_t rsmat_t;
    typedef typename Base:: dmat_t  dmat_t;
    typedef typename Base::   Ud_t    Ud_t;
    typedef typename Base::  UUd_t   UUd_t; //U,U',E  where U' has not been back transformed, U=V*Uprime.public:

    LASolverLapackEigen_blaze(const LAParams& lap) : LASolverLapackCommon_blaze<T>(lap) {};

    virtual void    SetBasisOverlap(const  smat_t& S);
    virtual rsmat_t         Inverse(const rsmat_t& S) const;

    using  LASolverCommon_blaze<T>::MakeSymmetric;
    using  LASolverCommon_blaze<T>::itsParams;

};

export template <class T> class LASolverLapackSVD_blaze
    : public virtual  LASolver_blaze<T>
    , private LASolverLapackCommon_blaze<T>
{
    typedef LASolver_blaze<T> Base;
    typedef typename Base:: rvec_t  rvec_t;
    typedef typename Base::  mat_t   mat_t;
    typedef typename Base:: smat_t  smat_t;
    typedef typename Base:: umat_t  umat_t;
    typedef typename Base:: lmat_t  lmat_t;
    typedef typename Base::rsmat_t rsmat_t;
    typedef typename Base:: dmat_t  dmat_t;
    typedef typename Base::   Ud_t    Ud_t;
    typedef typename Base::  UUd_t   UUd_t; //U,U',E  where U' has not been back transformed, U=V*Uprime.public:
public:
    LASolverLapackSVD_blaze(const LAParams& lap) : LASolverLapackCommon_blaze<T>(lap) {};

    virtual void    SetBasisOverlap(const  smat_t& S);
    virtual rsmat_t         Inverse(const rsmat_t& S) const;

    using  LASolverCommon_blaze<T>::MakeSymmetric;
    using  LASolverCommon_blaze<T>::itsParams;
};

export template <class T> class LASolverLapackCholsky_blaze
    : public virtual  LASolver_blaze<T>
    , private LASolverLapackCommon_blaze<T>
{
    typedef LASolver_blaze<T> Base;
    typedef typename Base:: rvec_t  rvec_t;
    typedef typename Base::  mat_t   mat_t;
    typedef typename Base:: smat_t  smat_t;
    typedef typename Base:: umat_t  umat_t;
    typedef typename Base:: lmat_t  lmat_t;
    typedef typename Base::rsmat_t rsmat_t;
    typedef typename Base:: dmat_t  dmat_t;
    typedef typename Base::   Ud_t    Ud_t;
    typedef typename Base::  UUd_t   UUd_t; //U,U',E  where U' has not been back transformed, U=V*Uprime.public:
public:
    LASolverLapackCholsky_blaze(const LAParams& lap) : LASolverLapackCommon_blaze<T>(lap) {};
    virtual void    SetBasisOverlap(const  smat_t& S);
    virtual rsmat_t         Inverse(const rsmat_t& S) const;
};
