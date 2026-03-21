// File: LASolver_blaze/LASolver.C  Linear algebra for Lowden orthogonalization and eigan solutions.
module;
#include <tuple>
#include "blaze/Math.h" 
export module qchem.LASolver_blaze;
export import qchem.LAParams;
export import oml;
export import qchem.Types;

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
//  For fitting charge density and exchange potential profile in DFT calculation we require
//  the inverse of wither the overlap and repulsion matrices.  This can be very unstable
//  for a high precision basis set that is close to being linearly dependent.  SVD + Truncate
//  can yield the so called Penrose inverse  S = U*s*Vd, Inv(S) = V * 1/s * Ud .
//  Similarly we can do an eigen inverse S=U*w*Ud, Inv(S)=U*1/w*Ud
//

export template <class T> class LASolver_blaze
{
protected:
    typedef blaze::DynamicMatrix<double,blaze::columnMajor> mat_t; 
    typedef blaze::SymmetricMatrix<mat_t> smat_t;
    typedef blaze::    UpperMatrix<mat_t> umat_t;
    typedef blaze::    LowerMatrix<mat_t> lmat_t;

    typedef blaze::  DynamicVector<double> rvec_t;
    typedef blaze::  DynamicMatrix<double,blaze::columnMajor> rmat_t;
    typedef blaze::SymmetricMatrix<rmat_t> rsmat_t;
    typedef blaze::DiagonalMatrix <rmat_t>  dmat_t; //For eigen and singular values.

public:
    virtual ~LASolver_blaze() {};
    
    typedef std::tuple<mat_t,rvec_t> Ud_t;
    typedef std::tuple<mat_t,mat_t,rvec_t> UUd_t; //U,U',E  where U' has not been back transformed, U=V*Uprime.

    static  LASolver_blaze* Factory(const LAParams&);  
      
    //! Factor S=U*Ud then invert U, V=U^-1  and store V and Vd.
    virtual void    SetBasisOverlap(const smat_t& S)=0;
    virtual rvec_t  Get_BS_Diagonal() const=0;
    virtual Ud_t    Solve(const smat_t& H) const=0; //H=Hamiltonian/Fock matrix.  H'=
    //! returns U,U',E  where U' has not been back transformed, U=V*Uprime.
    virtual UUd_t   SolveOrtho(const smat_t& Hprime) const=0; //Hprime = Vd * H * V Hamiltonian/Fock matrix.
    virtual rsmat_t Inverse(const rsmat_t& S) const=0;
    // For unit testing
    virtual rsmat_t Transform(const rsmat_t& M) const=0; // M' = Vd * M * V, where S = V * Vd
    virtual   mat_t Transform(const   mat_t& M) const=0; // M' = Vd * M * V, where S = V * Vd
};
