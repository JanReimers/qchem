// File: LASolver_blaze/Internal/LASolverImp.C  Linear algebra for Lowden orthogonalization and eigen solutions.
module;
#include <string>
#include "blaze/Math.h" 
export module qchem.LASolver_blaze.Internal.Common;
export import qchem.LASolver_blaze;

//#################################################################################
//
//  The goal here is to provide solutions to the generalized eigen problem:
//
//                H*U = e*S*U
//
//  Where H is the Hamiltonia matrix, U is a matrix of eigen vectors, e is a diagonal matrix
//  of eigen values, S is the orbital basis overlap matrix.  S is fixed at construction time
//  And U nd E are then calculated for various matrices H.  We need to find a decomposition
//  of S such that the general eigen problem is reduced to a normal eigen problem:
//
//                  H'*U'=e*S*U'
//
//  There are numerous ways of doing this: Eigen, SVD or Cholsky decomposition of S. These can 
//  be selected using the enums defined in LAPArams.  This version of the LASolver system is
//  hard coded to use blaze vectors and matrices.  Blaze calls to LAPACK are used for the linear
//  algebra work.
//


export template <class T> class LASolverCommon_blaze 
    : public virtual  LASolver_blaze<T>
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

    virtual rsmat_t  Transform(const rsmat_t& M) const;
    virtual mat_t    Transform(const   mat_t& M) const;
    virtual rvec_t   Get_BS_Diagonal() const {return Diag;}

protected:
    LASolverCommon_blaze(double truncationTolerance) : itsTruncationTolerance(truncationTolerance) {};
    ~LASolverCommon_blaze() {};
    //
    //  Helper functions used by derived classes.
    //
    static void Rescale (mat_t& V, const rvec_t& w);
    static void Rescale (mat_t& U, const rvec_t& w, mat_t& Vt);
    static void Truncate(mat_t& U, rvec_t& w, double tol);
    static void Truncate(mat_t& U, rvec_t& s, mat_t& V , double tol);
    static smat_t MakeSymmetric(mat_t&,std::string name);
    
    void AssignVs(const mat_t& _V, const mat_t& _Vd) {V=_V;Vd=_Vd;}

    double itsTruncationTolerance;
    mat_t V;
    mat_t Vd;  //Basis Overlap S = V*Vd
    rvec_t Diag; //s for SVD, e for eigen, diag for Cholsky.

};

