// File: LASolver/Internal/Imp/LASolverImp.C  Linear algebra for Lowden orthogonalization and eigne solutions.
module;
#include <cassert>
#include <iostream>
module qchem.LASolver.Internal.Common;
import qchem.Math;
import qchem.Blaze;

using std::cout;
using std::endl;
//-------------------------------------------------------------------------
//
//  Common level.
//
template <class T> typename LASolver<T>::Ud_t LASolverCommon<T>::Solve(const smat_t<T>& Ham) const
{
    assert(!blazem::isnan(Ham));
    mat_t<T> Hprime = Vd * Ham * V;  //Transform to orthogonal coordinates.
    smat_t<T> Hsym=MakeSymmetric(Hprime,"Hamiltonian");
    rvec_t d;
    mat_t<T>  U;
    blazem::eigen(Hsym,d,U);
    U = V * U;                      //Back transform.
    return std::make_tuple(U,d);
}
template <class T> typename LASolver<T>::UUd_t LASolverCommon<T>::SolveOrtho(const smat_t<T>& Hprime) const
{
    assert(!blazem::isnan(Hprime));
    rvec_t d;
    mat_t<T>  Uprime;
    blazem::eigen(Hprime,d,Uprime);
    // auto [Uprime,e]  =itsLapackEigenSolver->SolveAll(HPrime,itsParams.abstol);  //Get eigen solution.
    mat_t<T> U = V * Uprime;                      //Back transform.
    return std::make_tuple(U,Uprime,d);
}


template <class T> rsmat_t  LASolverCommon<T>::Transform(const rsmat_t& M) const
{
    mat_t<T> Mprime=Vd * M * V;  //Transform to orthogonal coordinates.
    return MakeSymmetric(Mprime,"Test matrix");
}
template <class T> mat_t<T>  LASolverCommon<T>::Transform(const mat_t<T>& M) const
{
    mat_t<T> Mprime=Vd * M * V;  //Transform to orthogonal coordinates.
    return Mprime;
}


template <class T>  smat_t<T> LASolverCommon<T>::MakeSymmetric(mat_t<T>& A,std::string name)
{
#ifdef false
    double del=::MakeSymmetric(A); // A=0.5*(A+~A)
    if (fabs(del) > 1e-9)
        std::cerr << "Warning: " << name << " asymmetry = " << del << " is big!" << std::endl;
#else
    A=0.5*(A+blazem::trans(A));
#endif
    return A;
}

 template class LASolverCommon<double>;



