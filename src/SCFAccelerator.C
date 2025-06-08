// FIle: SCFAccelerator.C  Interface for an accelerator alogrithm


#include "Imp/SCFAccelerator.H"
#include <LASolver.H>
#include "oml/numeric/LapackLinearSolver.H"
#include "oml/numeric/LapackSVDSolver.H"


#include "oml/diagonalmatrix.h"

SCFIrrepAccelerator_DIIS::SMat SCFIrrepAccelerator_DIIS::BuildB(const std::vector< Matrix<double>>& Es)
{
    index_t N=Es.size()+1;
    SMat B(N,N);
    for (index_t i:B.rows())
        for (index_t j:B.cols(i))
            if (i<N)
                if (j<N)
                    B(i,j)=Dot(Es[i-1],Es[j-1]);
                else
                    B(i,j)=1.0;
            else
                B(i,j)=0.0; //i=j=N
    
    return B;
}
double SCFIrrepAccelerator_DIIS::GetMinSV(const SMat& B)
{
    static oml::LapackSVDSolver<double> solver;
    auto [U,s,V]=solver.SolveAll(B);
    size_t N=s.GetNumRows();
    return s(N,N);
}
SCFIrrepAccelerator_DIIS::RVec SCFIrrepAccelerator_DIIS::SolveC(const SMat& B) 
{
    static LapackLinearSolver<double> solver;
    size_t N=B.GetNumRows();
    RVec v(N);
    Fill(v,0.0); 
    v(N)=1.0;
    RVec C=solver.Solve(B,v);
    // std:: cout << B << C << v << del <<std::endl;
    // std:: cout << "del,[F,D] = " << sqrt(del*del) << " " << C(N) << std::endl;
    return C.SubVector(N-1);   
}



SCFIrrepAccelerator_DIIS::SCFIrrepAccelerator_DIIS(const DIISParams& p) 
    : itsParams(p)
    , itsLastError(0.0)
    , itsLastSVMin(1e20)
    
{
};

SCFIrrepAccelerator_DIIS::~SCFIrrepAccelerator_DIIS() 
{

};

void SCFIrrepAccelerator_DIIS::Init(const LASolver<double>* las, const Irrep_QNs& qns)
{
    itsIrrep=qns;
    itsLaSolver=las;
}
#include "oml/vector.h"


void SCFIrrepAccelerator_DIIS::UseFD(const SMat& F, const SMat& DPrime)
{
    itsFPrime=itsLaSolver->Transform(F); // Fprime = Vd*F*V
    assert(itsFPrime.GetLimits()==DPrime.GetLimits());
    itsDPrime=DPrime;
}

SCFIrrepAccelerator::Mat SCFIrrepAccelerator_DIIS::CalculeteError() const
{
    return itsFPrime*itsDPrime-itsDPrime*itsFPrime;
}


SCFIrrepAccelerator::SMat SCFIrrepAccelerator_DIIS::Project(const SMat& F, const SMat& DPrime)
{
    SMat FPrime=itsLaSolver->Transform(F); // Fprime = Vd*F*V
    if (Max(fabs(DPrime))==0.0) return FPrime;
    assert(FPrime.GetLimits()==DPrime.GetLimits());
    Mat E=FPrime*DPrime-DPrime*FPrime;// Make the communtator
    itsLastError=FrobeniusNorm(E);
    if (itsLastError>itsParams.EMax) return FPrime;
  
    AppendAndPurge(FPrime,E,itsLastError,itsParams.Nproj);
    if (itsLastError< itsParams.EMin) return FPrime;
    RVec c=Solve();
    if (c.size()==0) return FPrime;
    return Project(c);

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


#include "Imp/Containers/stl_io.h"
SCFIrrepAccelerator_DIIS::RVec SCFIrrepAccelerator_DIIS::Solve()
{
    SMat B=BuildB(itsEs);
    size_t N=itsEs.size();
    while (GetMinSV(B)<itsParams.SVTol) 
    {
        if (N<=2)
        {
            itsLastSVMin=1e20;  //indicate bailout.
            return RVec();
        }
        N--;
        // std::cout << "Purging to N=" << N << std::endl;
        Purge(N);
        B=BuildB(itsEs);
    }
    
    return SolveC(B); //Solve for the projection coefficients
}

SCFIrrepAccelerator::SMat SCFIrrepAccelerator_DIIS::Project(const RVec& c)
{
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
        //auto iter=std::max_element(itsEns.begin(),itsEns.end()); //Find the maximum Error
        auto iter=itsEns.begin(); //just pick the oldest element.
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
void SCFIrrepAccelerator_DIIS::Purge(size_t N)
{
    while (itsEs.size()>N)
    {
        //auto iter=std::max_element(itsEns.begin(),itsEns.end()); //Find the maximum Error
        auto iter=itsEns.begin(); //just pick the oldest element.
        size_t index=std::distance(itsEns.begin(),iter);
        itsEns.erase(iter);
        itsEs.erase(itsEs.begin()+index);
        itsFPrimes.erase(itsFPrimes.begin()+index);
    }
    assert(itsEs.size()==N);
    assert(itsEns.size()==N);
    assert(itsFPrimes.size()==N);
    
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




void SCFAccelerator_DIIS::CalculateProjections()
{
    itsBailout=false;
    Mat E;
    for (auto k:itsIrreps) E+=k->CalculeteError();
    double En=FrobeniusNorm(E);
    if (itsBailout=En>itsParams.EMax;itsBailout) return;
    itsEs.push_back(E);
    if (itsBailout=itsEs.size()<2;itsBailout) return;
    SMat B=SCFIrrepAccelerator_DIIS::BuildB(itsEs);
    while (SCFIrrepAccelerator_DIIS::GetMinSV(B)<itsParams.SVTol) 
    {
        if (itsBailout=itsEs.size()<2;itsBailout) return;
        itsEs.erase(itsEs.begin());
        B=SCFIrrepAccelerator_DIIS::BuildB(itsEs);
    }
    itsCs=SCFIrrepAccelerator_DIIS::SolveC(B);
    for (auto k:itsIrreps) k->SetProjection(itsCs);
}

void SCFAccelerator_DIIS::ShowLabels(std::ostream& os) const
{
    os << " [F,D]   Nmin   Nmax    SVMin";
}
void SCFAccelerator_DIIS::ShowConvergence(std::ostream& os) const
{
    double EMax=0.0,SVMin=1.0e20;
    size_t NMin=1000,NMax=0;
    for (auto i:itsIrreps)
    {
        if (i->GetError()>EMax) EMax=i->GetError();
        if (i->GetSVMin()<SVMin) SVMin=i->GetSVMin();
        if (i->GetNproj()<NMin) NMin=i->GetNproj();
        if (i->GetNproj()>NMax) NMax=i->GetNproj();
        
    }
    os << std::scientific << std::setw(7) << std::setprecision(1) << EMax << " ";
    if (SVMin<1.0e20)
    {
        os << std::setw(3) << NMin << "    " << std::setw(3) << NMax << "     ";
        os << std::scientific << std::setw(7) << std::setprecision(1) << SVMin << "  ";
    }
    else
        os << "                        ";
}

double SCFAccelerator_DIIS::GetError() const
{
    double EMax=0.0;
    for (auto i:itsIrreps)
        if (i->GetError()>EMax) EMax=i->GetError();
    return EMax;
}


