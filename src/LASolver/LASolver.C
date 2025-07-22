// File: LASolver.C  Linear algebra for Lowden orthogonalization and eigne solutions.
module;
#include <tuple>
export module qchem.LASolver;
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


export template <class T> class LASolver
{
public:
    virtual ~LASolver() {};
    
    typedef Matrix<T>      Mat;
    typedef SMatrix<T>     SMat;
    typedef SMatrix<double> RSMat;
    typedef std::tuple<Mat,RVec> UdType;
    typedef std::tuple<Mat,Mat,RVec> UUdType; //U,U',E  where U' has not been back transformed, U=V*Uprime.

    static  LASolver* Factory(const LAParams&);  
      
    //! Factor S=U*Ud then invert U, V=U^-1  and store V and Vd.
    virtual void   SetBasisOverlap(const SMat& S)=0;
    virtual RVec   Get_BS_Diagonal() const=0;
    virtual UdType Solve(const SMat& H) const=0; //H=Hamiltonian/Fock matrix.  H'=
    //! returns U,U',E  where U' has not been back transformed, U=V*Uprime.
    virtual UUdType SolveOrtho(const SMat& Hprime) const=0; //Hprime = Vd * H * V Hamiltonian/Fock matrix.
    virtual RSMat  Inverse(const RSMat& S) const=0;
    // For unit testing
    virtual RSMat  Transform(const RSMat& M) const=0; // M' = Vd * M * V, where S = V * Vd
    virtual Mat    Transform(const   Mat& M) const=0; // M' = Vd * M * V, where S = V * Vd
};
