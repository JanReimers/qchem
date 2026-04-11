// File: LASolver_blaze/Internal/LASolverImp.C  Linear algebra for Lowden orthogonalization and eigen solutions.
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
    virtual rvec_t   Get_BS_Diagonal() const {return Diag;}

protected:
    LASolverCommon(double truncationTolerance) : itsTruncationTolerance(truncationTolerance) {};
    ~LASolverCommon() {};
    //
    //  Helper functions used by derived classes.
    //
    static void Rescale (mat_t<T>& V, const rvec_t& w);
    static void Rescale (mat_t<T>& U, const rvec_t& w, mat_t<T>& Vt);
    static void Truncate(mat_t<T>& U, rvec_t& w, double tol);
    static void Truncate(mat_t<T>& U, rvec_t& s, mat_t<T>& V , double tol);
    static smat_t<T> MakeSymmetric(mat_t<T>&,std::string name);
    
    void AssignVs(const mat_t<T>& _V, const mat_t<T>& _Vd) {V=_V;Vd=_Vd;}

    double itsTruncationTolerance;
    mat_t<T> V;
    mat_t<T> Vd;  //Basis Overlap S = V*Vd
    rvec_t Diag; //s for SVD, e for eigen, diag for Cholsky.

};

