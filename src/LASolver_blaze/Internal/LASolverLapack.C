// File: LASolver_blaze/Internal/LASolverLapack.C  Use blaze interface to Lapack linrary for all linear algebra ops.
module;
#include "blaze/Math.h" 
export module qchem.LASolver_blaze.Internal.Lapack;
import qchem.LASolver.Internal.Common;


export template <class T> class LASolverLapackCommon_blaze
    : public virtual  LASolver<T>
    , protected LASolverCommon<T>
{
    typedef LASolver<T> Base;
    typedef blaze::   UpperMatrix< mat_t<T>> umat_t;
    typedef blaze::   LowerMatrix< mat_t<T>> lmat_t;
    typedef blaze::DiagonalMatrix<rmat_t   > dmat_t; //For eigen and singular values.
    typedef typename Base::   Ud_t    Ud_t;
    typedef typename Base::  UUd_t   UUd_t; //U,U',E  where U' has not been back transformed, U=V*Uprime.
public:
    LASolverLapackCommon_blaze(double truncationTolerance);
    ~LASolverLapackCommon_blaze();
    
    virtual  Ud_t Solve     (const smat_t<T>& H) const;
    //! returns U,U',E  where U' has not been back transformed, U=V*Uprime.
    virtual UUd_t SolveOrtho(const smat_t<T>& Hprime) const; //Hprime = Vd * H * V Hamiltonian/Fock matrix.
   
    using  LASolverCommon<T>::MakeSymmetric;
    using  LASolverCommon<T>::itsTruncationTolerance;

    using LASolverCommon<T>::V;
    using LASolverCommon<T>::Vd;
};

export template <class T> class LASolverLapackEigen_blaze
    : public virtual  LASolver<T>
    , private LASolverLapackCommon_blaze<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::   Ud_t    Ud_t;
    typedef typename Base::  UUd_t   UUd_t; //U,U',E  where U' has not been back transformed, U=V*Uprime.public:
public:
    LASolverLapackEigen_blaze(double truncationTolerance) : LASolverLapackCommon_blaze<T>(truncationTolerance) {};

    virtual void    SetBasisOverlap(const  smat_t<T>& S);
    virtual rsmat_t         Inverse(const rsmat_t& S) const;

    using  LASolverCommon<T>::MakeSymmetric;
    using  LASolverCommon<T>::itsTruncationTolerance;

};

export template <class T> class LASolverLapackSVD_blaze
    : public virtual  LASolver<T>
    , private LASolverLapackCommon_blaze<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::   Ud_t    Ud_t;
    typedef typename Base::  UUd_t   UUd_t; //U,U',E  where U' has not been back transformed, U=V*Uprime.public:
public:
    LASolverLapackSVD_blaze(double truncationTolerance) : LASolverLapackCommon_blaze<T>(truncationTolerance) {};

    virtual void    SetBasisOverlap(const  smat_t<T>& S);
    virtual rsmat_t         Inverse(const rsmat_t& S) const;

    using  LASolverCommon<T>::MakeSymmetric;
    using  LASolverCommon<T>::itsTruncationTolerance;
};

export template <class T> class LASolverLapackCholsky_blaze
    : public virtual  LASolver<T>
    , private LASolverLapackCommon_blaze<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::   Ud_t    Ud_t;
    typedef typename Base::  UUd_t   UUd_t; //U,U',E  where U' has not been back transformed, U=V*Uprime.public:
public:
    LASolverLapackCholsky_blaze(double truncationTolerance) : LASolverLapackCommon_blaze<T>(truncationTolerance) {};
    virtual void    SetBasisOverlap(const  smat_t<T>& S);
    virtual rsmat_t         Inverse(const rsmat_t& S) const;
};
