// File: LAParams.H
#ifndef _LAParams_H_
#define _LAParams_H_

namespace qchem {
    //
    //  OML code is derived from the Numerical Recipes FORTRAN code
    //  adapted for OML Matrix and Vector containers.
    //  Lapack is wrapped to work with OML containers.
    //
    enum Pkg {OML,Lapack}; //Use linear algebra code from which package
    //
    //  Cholsky S = L*U
    //  Eigen   S = U*e*U_dagger
    //  SVD     S = U*s*V_dagger 
    //
    enum Ortho {Cholsky,Eigen,SVD}; //Basis set orthogonalization technique
    
} //Namespace qchem

struct LAParams
{
    //LAParams() : LinearAlgebraPackage(qchem::Lapack), BasisOrthoAlgorithm(qchem::SVD), TruncationTolerance(1e-10), abstol(1e-12) {};
    qchem::Pkg    LinearAlgebraPackage; //{OML,Lapack}
    qchem::Ortho  BasisOrthoAlgorithm;  //{Cholsky,Eigen,SVD}
    //For Ortho=SVD/Eigen we can truncation solutions with small
    //Eigen/singular values.  THis remove linear dependency.
    double        TruncationTolerance;  
    //Lapack eigen convergence tolerance.   See details below.  Not used for SVD.                               
    double abstol;               
};
//  Pasted form the Lapack documentation https://netlib.org/lapack//explore-html/d4/de0/group__heevx_gaa457d630a97fca067f18f73de0dfce01.html

/*
   ABSTOL is DOUBLE PRECISION
          The absolute error tolerance for the eigenvalues.
          An approximate eigenvalue is accepted as converged
          when it is determined to lie in an interval [a,b]
          of width less than or equal to

                  ABSTOL + EPS *   max( |a|,|b| ) ,

          where EPS is the machine precision.  If ABSTOL is less than
          or equal to zero, then  EPS*|T|  will be used in its place,
          where |T| is the 1-norm of the tridiagonal matrix obtained
          by reducing A to tridiagonal form.

          Eigenvalues will be computed most accurately when ABSTOL is
          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
          If this routine returns with INFO>0, indicating that some
          eigenvectors did not converge, try setting ABSTOL to
          2*DLAMCH('S').

          See "Computing Small Singular Values of Bidiagonal Matrices
          with Guaranteed High Relative Accuracy," by Demmel and
          Kahan, LAPACK Working Note #3.
*/


#endif // _LAParams_H_
