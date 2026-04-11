// File: LASolver_blaze/Internal/LASolverLapack.C  Use blaze interface to Lapack linrary for all linear algebra ops.
module;
#include "blaze/Math.h" 
export module qchem.LASolver_blaze.Internal.Lapack;
import qchem.LASolver.Internal.Common;

export template <class T> class LASolverEigen
    : public virtual  LASolver<T>
    , private LASolverCommon<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::   Ud_t    Ud_t;
    typedef typename Base::  UUd_t   UUd_t; //U,U',E  where U' has not been back transformed, U=V*Uprime.public:
public:
    LASolverEigen(double truncationTolerance) : LASolverCommon<T>(truncationTolerance) {};

    virtual void    SetBasisOverlap(const  smat_t<T>& S);
    virtual rsmat_t         Inverse(const rsmat_t& S) const;

    using  LASolverCommon<T>::MakeSymmetric;
    using  LASolverCommon<T>::itsTruncationTolerance;
private:
    static void Rescale (mat_t<T>& V, const rvec_t& w);
    static void Truncate(mat_t<T>& U, rvec_t& w, double tol);

};

export template <class T> class LASolverSVD
    : public virtual  LASolver<T>
    , private LASolverCommon<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::   Ud_t    Ud_t;
    typedef typename Base::  UUd_t   UUd_t; //U,U',E  where U' has not been back transformed, U=V*Uprime.public:
public:
    LASolverSVD(double truncationTolerance) : LASolverCommon<T>(truncationTolerance) {};

    virtual void    SetBasisOverlap(const  smat_t<T>& S);
    virtual rsmat_t         Inverse(const rsmat_t& S) const;

    using  LASolverCommon<T>::MakeSymmetric;
    using  LASolverCommon<T>::itsTruncationTolerance;
private:    
    static void Rescale (mat_t<T>& U, const rvec_t& w, mat_t<T>& Vt);
    static void Truncate(mat_t<T>& U, rvec_t& s, mat_t<T>& V , double tol);

};

export template <class T> class LASolverCholsky
    : public virtual  LASolver<T>
    , private LASolverCommon<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::   Ud_t    Ud_t;
    typedef typename Base::  UUd_t   UUd_t; //U,U',E  where U' has not been back transformed, U=V*Uprime.public:
public:
    LASolverCholsky(double truncationTolerance) : LASolverCommon<T>(truncationTolerance) {};
    virtual void    SetBasisOverlap(const  smat_t<T>& S);
    virtual rsmat_t         Inverse(const rsmat_t& S) const;

};
