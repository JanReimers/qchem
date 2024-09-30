// File: EigenSolver.H  General eigen solver.
#ifndef _EigenSolver_H_
#define _EigenSolver_H_

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

#include "types.H"
#include "oml/matrix.h"

template <class T> class EigenSolver
{
public:
    virtual ~EigenSolver() {};
    
    typedef Matrix<T>      Mat;
    typedef SMatrix<T>     SMat;
    typedef Vector<T>      Vec;
    typedef Vector<double> RVec;
    typedef std::tuple<Mat,RVec> UdType;
    typedef std::tuple<Mat,Mat>  VVType;
    
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
    enum Ortho {Cholsky,Eigen,SVD}; //Basis set orthogonalization techniquw

    static EigenSolver* Factory(Pkg,Ortho,const SMat& S, double Tolerance);    

    virtual UdType Solve(const SMat& H) const=0;
};

template <class T> class EigenSolverCommon 
    : public virtual  EigenSolver<T>
{
    typedef EigenSolver<T> Base;
    typedef typename Base::RVec RVec;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
    typedef typename Base::VVType VVType;

protected:
    static void Rescale(Mat& V, const RVec& w);
    void AssignVs(const Mat& _V, const Mat& _Vd) {V=_V;Vd=_Vd;}
    Mat  V,Vd;  //Basis Overlap S = V*Vd
};

template <class T> class EigenSolverOMLCommon
    : public virtual  EigenSolver<T>
    , protected EigenSolverCommon<T>
{
    typedef EigenSolver<T> Base;
    typedef typename Base::Vec Vec;
    typedef typename Base::RVec RVec;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
    typedef typename Base::UdType UdType;
public:
    virtual UdType Solve(const SMat&) const;
    
    using EigenSolverCommon<T>::V;
    using EigenSolverCommon<T>::Vd;
};

template <class T> class EigenSolverOMLEigen
    : public virtual  EigenSolver<T>
    , private EigenSolverOMLCommon<T>
{
    typedef EigenSolver<T> Base;
    typedef typename Base::Vec Vec;
    typedef typename Base::RVec RVec;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
    typedef typename Base::UdType UdType;
public:
    EigenSolverOMLEigen(const SMat& S, double tolerance);
};

template <class T> class EigenSolverOMLSVD
    : public virtual  EigenSolver<T>
    , private EigenSolverOMLCommon<T>
{
    typedef EigenSolver<T> Base;
    typedef typename Base::Vec Vec;
    typedef typename Base::RVec RVec;
    typedef typename Base::Mat Mat;
    typedef typename Base::SMat SMat;
    typedef typename Base::UdType UdType;
public:
    EigenSolverOMLSVD(const SMat& S, double tolerance);
};


#endif //_EigenSolver_H_
