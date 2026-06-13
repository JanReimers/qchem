// File: LASolver/Internal/LASolverImp.C  Linear algebra for Lowden orthogonalization and eigen solutions.
module;
#include <string>
#include "blaze/Math.h" 
export module qchem.LASolver.Internal.Common;
export import qchem.LASolver;

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


export template <class T> class LASolverCommon 
    : public virtual  LASolver<T>
{

    virtual rsmat_t  Transform(const rsmat_t& M) const;
    virtual mat_t<T>    Transform(const   mat_t<T>& M) const;
    virtual mat_t<T>    BackTransform(const mat_t<T>& Uprime) const {return V*Uprime;}
    virtual rvec_t   Get_BS_Diagonal() const {return Diag;}
    typedef LASolver<T> Base;
    typedef blaze::   UpperMatrix< mat_t<T>> umat_t;
    typedef blaze::   LowerMatrix< mat_t<T>> lmat_t;
    typedef blaze::DiagonalMatrix<rmat_t   > dmat_t; //For eigen and singular values.
    typedef typename Base::   Ud_t    Ud_t;
    typedef typename Base::  UUd_t   UUd_t; //U,U',E  where U' has not been back transformed, U=V*Uprime.

protected:
    LASolverCommon(double truncationTolerance) : itsTruncationTolerance(truncationTolerance) {};
    ~LASolverCommon() {};
    virtual  Ud_t Solve     (const smat_t<T>& H) const;
    //! returns U,U',E  where U' has not been back transformed, U=V*Uprime.
    virtual UUd_t SolveOrtho(const smat_t<T>& Hprime) const; //Hprime = Vd * H * V Hamiltonian/Fock matrix.
    //
    //  Helper functions used by derived classes.
    //
    static smat_t<T> MakeSymmetric(mat_t<T>&,std::string name);
    
    void AssignVs(const mat_t<T>& _V, const mat_t<T>& _Vd) {V=_V;Vd=_Vd;}

    double itsTruncationTolerance;
    mat_t<T> V;
    mat_t<T> Vd;  //Basis Overlap S = V*Vd
    rvec_t Diag; //s for SVD, e for eigen, diag for Cholsky.

};

