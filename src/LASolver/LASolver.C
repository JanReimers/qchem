// File: LASolver/LASolver.C  Linear algebra for Lowden orthogonalization and eigen solutions.
module;
#include <tuple>
export module qchem.LASolver;
export import qchem.Types;
import qchem.Blaze;

export namespace qchem
{
    enum Ortho {Cholesky, Eigen, SVD};
}

//#################################################################################
//
//  The goal here is to provide solutions to the generalized eigen problem:
//
//                H*C = e*S*C
//
//  Where H is the Hamiltonian matrix, C is a matrix of eigen vectors, e is a diagonal matrix
//  of eigen values, S is the orbital basis overlap matrix.  S is fixed with a call to SetBasisOverlap
//  And C and e are then calculated for various matrices H.  We need to find a decomposition
//  of S such that the general eigen problem is reduced to a normal eigen problem:
//
//                  H'*C'=e*C'
//
//  There are numerous ways of doing this: Eigen, SVD or Cholesky decomposition of S. These can
//  be selected using the enum qchem::Ortho.  This interface is hard coded to use blaze matrix
//  classes.  Under the hood blaze uses the Lapack and Blas libraries for linear algebra operations.
//
//  hmat_t<double> = blaze::SymmetricMatrix, hmat_t<dcmplx> = blaze::HermitianMatrix.
//

export template <class T> class LASolver
{
public:
    virtual ~LASolver() {};

    typedef std::tuple<mat_t<T>         ,rvec_t> Ud_t;
    typedef std::tuple<mat_t<T>,mat_t<T>,rvec_t> UUd_t; //U,U',E  where U' has not been back transformed, U=V*Uprime.

    static  LASolver* Factory(qchem::Ortho ortho, double TruncationTolerance);

    //! Factor S=U*Ud then invert U, V=U^-1  and store V and Vd.
    virtual void       SetBasisOverlap(const hmat_t<T>& S)=0;
    virtual rvec_t     Get_BS_Diagonal() const=0;
    virtual Ud_t       Solve(const hmat_t<T>& H) const=0;       //! H' = Vd*H*V, eigen-solve, back-transform.
    //! returns U,U',E  where U' has not been back transformed, U=V*Uprime.
    virtual UUd_t      SolveOrtho(const hmat_t<T>& Hprime) const=0; //! Hprime = Vd * H * V
    virtual hmat_t<T>  Transform(const hmat_t<T>& M) const=0;  //! M' = Vd * M * V
    //! Back-transform orthonormal-basis coefficients to the AO basis: U = V * Uprime.
    virtual mat_t<T>   BackTransform(const mat_t<T>& Uprime) const=0;
};
