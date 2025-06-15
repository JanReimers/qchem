// File: LASolverImp.H  General eigen solver.
#ifndef _LASolverImp_H_
#define _LASolverImp_H_

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
//  be selected using the enums defined below.  Also one can chose the linear algebra library
//  to use: OML (derived from Numerical Recipes), or Lapack.
//

#include <LASolver.H>
#include <LAParams.H>
#include "oml/matrix.h"

template <class T> class LASolverCommon 
    : public virtual  LASolver<T>
{
    typedef LASolver<T> Base;
    typedef typename Base::RVec RVec;
    typedef typename Base::Mat  Mat;
    typedef typename Base::SMat SMat;
    typedef typename Base::RSMat RSMat;
    typedef DiagonalMatrix<T>  DMat;

    virtual RSMat  Transform(const RSMat& M) const;
    virtual Mat    Transform(const   Mat& M) const;
    virtual RVec   Get_BS_Diagonal() const {return Diag;}

protected:
    LASolverCommon(const LAParams& lap) : itsParams(lap) {};
    ~LASolverCommon() {};
    //
    //  Helper functions used by derived classes.
    //
    static void Rescale (Mat& V, const RVec& w);
    static void Rescale (Mat& U, const RVec& w, Mat& Vt);
    static void Truncate(Mat& U, RVec& w, double tol);
    static void Truncate(Mat& U, RVec&, Mat& V , double tol);
    static SMat MakeSymmetric(Mat&,c_str name);
    
    void AssignVs(const Mat& _V, const Mat& _Vd) {V=_V;Vd=_Vd;}

    LAParams itsParams;
    Mat  V,Vd;  //Basis Overlap S = V*Vd
    RVec Diag; //s for SVD, e for eigen, diag for Cholsky.
};


#endif //_LASolverImp_H_
