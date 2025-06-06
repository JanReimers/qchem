// FIle: SCFAccelerator.C  Interface for an accelerator alogrithm


#include "Imp/SCFAccelerator.H"
#include <LASolver.H>
#include "oml/numeric/LapackLinearSolver.H"
#include "oml/numeric/LapackSVDSolver.H"


SCFIrrepAccelerator_DIIS::SCFIrrepAccelerator_DIIS(const DIISParams& p) 
    : itsParams(p)
    , itsLastError(0.0)
    , itsLastSVMin(1e20)
    , itsLinearSolver(new LapackLinearSolver<double>())
    , itsSVDSolver(new oml::LapackSVDSolver<double>())
{
    assert(itsLinearSolver);
    assert(itsSVDSolver);
};

SCFIrrepAccelerator_DIIS::~SCFIrrepAccelerator_DIIS() 
{
    delete itsLinearSolver;
    delete itsSVDSolver;
};

void SCFIrrepAccelerator_DIIS::Init(const LASolver<double>* las)
{
    itsLaSolver=las;
}

SCFIrrepAccelerator::SMat SCFIrrepAccelerator_DIIS::Project(const SMat& F, const SMat& DPrime)
{
    SMat FPrime=itsLaSolver->Transform(F); // Fprime = Vd*F*V
    if (Max(fabs(DPrime))==0.0) return FPrime;
    Mat E=FPrime*DPrime-DPrime*FPrime;// Make the communtator
    itsLastError=FrobeniusNorm(E);
    if (itsLastError<itsParams.EMax  && itsLastError> itsParams.EMin)
    {
        // std::cout << "||E|| = " << FrobeniusNorm(E) << std:: endl;
        AppendAndPurge(FPrime,E,itsLastError,itsParams.Nproj);
        return Project(FPrime);
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

#include "oml/diagonalmatrix.h"
bool SCFIrrepAccelerator_DIIS::IsSingular(const SMat& B, double SVtol) 
{
    auto [U,s,V]=itsSVDSolver->SolveAll(B);
    size_t N=s.GetNumRows();
    itsLastSVMin=s(N,N);
    return itsLastSVMin<itsParams.SVTol;
}

#include "oml/vector.h"

SCFIrrepAccelerator_DIIS::RVec SCFIrrepAccelerator_DIIS::SolveC(const SMat& B) const
{
    size_t N=B.GetNumRows();
    RVec v(N);
    Fill(v,0.0); 
    v(N)=1.0;
    Vector<double> C=itsLinearSolver->Solve(B,v);
    // std:: cout << B << C << v << del <<std::endl;
    // std:: cout << "del,[F,D] = " << sqrt(del*del) << " " << C(N) << std::endl;
    return C.SubVector(N-1);   
}

#include "Imp/Containers/stl_io.h"
SCFIrrepAccelerator::SMat SCFIrrepAccelerator_DIIS::Project(const SMat& Fprime)
{
    SMat B=BuildB();
    if (IsSingular(B,itsParams.SVTol)) 
    {
        itsLastSVMin=1e20;  //indicate bailout.
        return Fprime;
    }
    RVec c=SolveC(B); //Solve for the projection coefficients
    assert(fabs(Sum(c)-1.0)<1e-13); //Check that the constraint worked.

    // Now do the projection for the Fock matrix.
    SMat Fproj;
    index_t i=1;
    for (const auto& f:itsFPrimes) Fproj+=SMat(c(i++)*f);
    return Fproj;

}

#include <algorithm>
void SCFIrrepAccelerator_DIIS::AppendAndPurge(const SMat& FPrime, const Mat& E, double En, size_t N)
{
    assert(itsEs.size()==itsFPrimes.size());
    assert(itsEns.size()==itsFPrimes.size());
    // Purge first to avoid purging the latest
    while (N>0 && itsEs.size()>N-1)
    {
        auto iter=std::max_element(itsEns.begin(),itsEns.end()); //Find the maximum Error
        size_t index=std::distance(itsEns.begin(),iter);
        itsEns.erase(iter);
        itsEs.erase(itsEs.begin()+index);
        itsFPrimes.erase(itsFPrimes.begin()+index);
    }
    itsEs.push_back(E);
    itsEns.push_back(En);
    itsFPrimes.push_back(FPrime);
    assert(itsEs.size()<=N);
    assert(itsEns.size()<=N);
    assert(itsFPrimes.size()<=N);
    
}


#include <Irrep_BS.H>

SCFAccelerator_DIIS::SCFAccelerator_DIIS(const DIISParams& p) 
: itsParams(p)
{};

SCFAccelerator_DIIS::~SCFAccelerator_DIIS() {};
SCFIrrepAccelerator* SCFAccelerator_DIIS::Create(const TOrbital_IBS<double>* bs) 
{
    itsIrreps.push_back(new SCFIrrepAccelerator_DIIS(itsParams));
    return itsIrreps.back();;
}

void SCFAccelerator_DIIS::ShowLabels(std::ostream& os) const
{
    os << " [F,D]   SVMin";
}
void SCFAccelerator_DIIS::ShowConvergence(std::ostream& os) const
{
    double EMax=0.0,SVMin=1.0e20;
    for (auto i:itsIrreps)
    {
        if (i->GetError()>EMax) EMax=i->GetError();
        if (i->GetSVMin()<SVMin) SVMin=i->GetSVMin();
    }
    os << std::scientific << std::setw(7) << std::setprecision(1) << EMax << " ";
    if (SVMin<1.0e20)
        os << std::scientific << std::setw(7) << std::setprecision(1) << SVMin << " ";
    else
        os << "        ";
}

double SCFAccelerator_DIIS::GetError() const
{
    double EMax=0.0;
    for (auto i:itsIrreps)
        if (i->GetError()>EMax) EMax=i->GetError();
    return EMax;
}


