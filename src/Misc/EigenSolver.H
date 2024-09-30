// File: EigenSolver.H  General eigen solver.
#ifndef _EigenSolver_H_
#define _EigenSolver_H_



//#################################################################################
//
//  The goal here is to provide solutions to
//
//                HU = ESU
//
//  Where H is the Hamiltonia matrix, U is a matrix of eigen vectors, E is a diaginal matrix
//  of eigen values, S is the orbital basis overlap matrix.  S is fixed at construction time
//  And U nd E are then calculated for various matricies H.  If we diagonalize S = s+ lambda s
//  and define V=sqrt(lamnda) S,  then we can write S=V+ V.  so V is a square root of S and
//  represents linear comniations of the basis functions that are ortho-normal.
//
//

#include "oml/matrix.h"
#include "oml/vector.h"

template <class T> class TBasisSet;
template <class T> class SMatrix;

template <class T> class EigenSolver
{
public:
    typedef Matrix<T>      Mat;
    typedef SMatrix<T>     SMat;
    typedef Vector<T>      Vec;
    typedef Vector<double> RVec;
    
    //
    //  OML code is derived from the Numerical Recipes FORTRAN code
    //    adapted for OML Matrix and Vector containers.
    //  Lapack is wrapped to work with OML containers.
    //
    enum Pkg {OML,Lapack}; //Use linear algebra code from which package
    //
    //  Cholsky S = L*U
    //  Eigen   S = U*e*U_dagger
    //  SVD     S = U*s*V_dagger with the option for truncation singular values below 
    //                           some threshold.
    //
    enum Ovp {Cholsky,Eigen,SVD}; //Basis set orthogonalization techniquw

    EigenSolver(const TBasisSet<T>&, double Tolerance);

    void  Solve          (const SMat&)      ;
    const Mat&  GetEigenVectors(           ) const;
    const RVec& GetEigenValues (           ) const;

private:
    Mat  V;
    Mat  EigenVectors;
    RVec EigenValues ;
};

#endif //_EigenSolver_H_
