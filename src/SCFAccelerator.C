// FIle: SCFAccelerator.C  Interface for an accelerator alogrithm


#include "Imp/SCFAccelerator.H"
#include <LASolver.H>

SCFIrrepAccelerator_DIIS::SCFIrrepAccelerator_DIIS(const SMat& _S) : S(_S) {};
SCFIrrepAccelerator_DIIS::~SCFIrrepAccelerator_DIIS() {};
void SCFIrrepAccelerator_DIIS::Init(const LASolver<double>* las)
{
    itsLaSolver=las;
}

SCFIrrepAccelerator::SMat SCFIrrepAccelerator_DIIS::Project(const SMat& F, const SMat& DPrime)
{
    SMat FPrime=itsLaSolver->Transform(F); // Fprime = Vd*F*V
    if (Max(fabs(DPrime))==0.0) return FPrime;
    Mat E=FPrime*DPrime-DPrime*FPrime;// Make the communtator
    double En=FrobeniusNorm(E);
    if (En<1.0  && En> 1e-7)
    {
        // std::cout << "||E|| = " << FrobeniusNorm(E) << std:: endl;
        itsEs.push_back(E);
        itsEns.push_back(En);
        itsF_Primes.push_back(FPrime);
        Purge(8);
        return Solve(FPrime);
    }
    else
        return FPrime;
}




template <class T> const SMatrix<T>& operator+=(SMatrix<T>& a, const SMatrix<T>& b)
{
    if (a.size()==0) 
    {
        a.SetLimits(b.GetLimits());
        Fill(a,0.0);
    }
    return ArrayAdd(a,b);
}

SCFIrrepAccelerator::SMat SCFIrrepAccelerator_DIIS::BuildB() const
{
    index_t N=itsEs.size()+1;
    SMat B(N,N);
    for (index_t i:B.rows())
        for (index_t j:B.cols(i))
            if (i<N)
                if (j<N)
                    B(i,j)=Dot(itsEs[i-1],itsEs[j-1]);
                else
                    B(i,j)=1.0;
            else
                B(i,j)=0.0; //i=j=N
    
    return B;
}

#include "oml/numeric/LapackSVDSolver.H"
#include "oml/diagonalmatrix.h"
bool SCFIrrepAccelerator_DIIS::IsSingular(const SMat& B, double SVtol) const
{
    auto svd_solver=new oml::LapackSVDSolver<double>();
    auto [U,s,V]=svd_solver->SolveAll(B);
    delete svd_solver;
    size_t N=s.GetNumRows();
    return s(N,N)<SVtol;
}

#include "oml/vector.h"
#include "oml/numeric/LapackLinearSolver.H"

SCFIrrepAccelerator_DIIS::RVec SCFIrrepAccelerator_DIIS::SolveC(const SMat& B) const
{
    size_t N=B.GetNumRows();
    RVec v(N);
    Fill(v,0.0); 
    v(N)=1.0;
    auto solver=new LapackLinearSolver<double>();
    Vector<double> C=solver->Solve(B,v);
    delete solver;
    // std:: cout << B << C << v << del <<std::endl;
    // std:: cout << "del,[F,D] = " << sqrt(del*del) << " " << C(N) << std::endl;
    return C.SubVector(N-1);   
}

#include "Imp/Containers/stl_io.h"
SCFIrrepAccelerator::SMat SCFIrrepAccelerator_DIIS::Solve(const SMat& Fprime)
{
    SMat B=BuildB();
    if (IsSingular(B,1e-9)) return Fprime;
    RVec c=SolveC(B); //Solve for the projection coefficients
    assert(fabs(Sum(c)-1.0)<1e-13); //Check that the constraint worked.

    // Now do the projection for the Fock matrix.
    SMat Fproj;
    index_t i=1;
    for (const auto& f:itsF_Primes) Fproj+=SMat(c(i++)*f);
    return Fproj;

}

#include <algorithm>
void SCFIrrepAccelerator_DIIS::Purge(int N)
{
    assert(itsEs.size()==itsF_Primes.size());
    while ((int)itsEs.size()>N)
    {
        auto iter=std::max_element(itsEns.begin(),itsEns.end());
        size_t index=std::distance(itsEns.begin(),iter);
        itsEns.erase(iter);
        itsEs.erase(itsEs.begin()+index);
        itsF_Primes.erase(itsF_Primes.begin()+index);
    }
    
}


#include <Irrep_BS.H>

SCFAccelerator_DIIS::SCFAccelerator_DIIS() {};
SCFAccelerator_DIIS::~SCFAccelerator_DIIS() {};
SCFIrrepAccelerator* SCFAccelerator_DIIS::Create(const TOrbital_IBS<double>* bs) const
{
    return new SCFIrrepAccelerator_DIIS(bs->Overlap());
}

