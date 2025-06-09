// FIle: SCFAccelerator.C  Interface for an accelerator alogrithm


#include "Imp/SCFAccelerator.H"
#include <LASolver.H>
#include "oml/numeric/LapackLinearSolver.H"
#include "oml/numeric/LapackSVDSolver.H"

using std::cout;
using std::endl;

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
    , itsLastEn(0.0)
    , itsLastSVMin(1e20)
    , itsCs(1)
{
    itsCs(1)=1.0;
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
    itsBailout = Max(fabs(itsDPrime))==0.0;
}

SCFIrrepAccelerator::Mat SCFIrrepAccelerator_DIIS::CalculateError() const
{
    return itsFPrime*itsDPrime-itsDPrime*itsFPrime;
}


SCFIrrepAccelerator::SMat SCFIrrepAccelerator_DIIS::Project(const SMat& F, const SMat& DPrime)
{
    itsBailout=false;
    SMat FPrime=itsLaSolver->Transform(F); // Fprime = Vd*F*V
    if (itsBailout=Max(fabs(DPrime))==0.0;itsBailout) return FPrime;
    assert(FPrime.GetLimits()==DPrime.GetLimits());
    itsLastE=FPrime*DPrime-DPrime*FPrime;// Make the communtator
    itsLastEn=FrobeniusNorm(itsLastE);
    if (itsBailout=itsLastEn>itsParams.EMax;itsBailout) return FPrime;
  
    Append(FPrime,itsLastE,itsLastEn);
    if (itsEs.size()>itsParams.Nproj) Purge1();
    assert(itsEs.size()<=itsParams.Nproj);
    if (itsBailout=itsLastEn< itsParams.EMin;itsBailout) return FPrime;
    RVec c=Solve();
    if (itsBailout=c.size()<2;itsBailout) return FPrime;
    return Project(c);

}

SCFIrrepAccelerator::SMat SCFIrrepAccelerator_DIIS::Project()
{
    if (itsBailout) 
        return itsFPrime;
    else
        return Project(itsCs);
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
template <class T> const Matrix<T>& operator+=(Matrix<T>& a, const Matrix<T>& b)
{
    if (a.size()==0) 
    {
        a.SetLimits(b.GetLimits());
        Fill(a,0.0);
    }
    return ArrayAdd(a,b);
}

SCFIrrepAccelerator_DIIS::md_t SCFIrrepAccelerator_DIIS::BuildPrunedB(const mv_t& Es,double svmin)
{
    SMat B=BuildB(Es);
    double sv=GetMinSV(B);
    while (sv<svmin && Es.size()>=2) 
    {
        Purge1(); //Must be a member function for this.
        B=BuildB(Es);
        sv=GetMinSV(B);
    }
    return std::make_pair(B,sv);
}

#include "Imp/Containers/stl_io.h"
SCFIrrepAccelerator_DIIS::RVec SCFIrrepAccelerator_DIIS::Solve()
{
    SMat B;
    std::tie(B,itsLastSVMin)=BuildPrunedB(itsEs,itsParams.SVTol);
    if (B.GetNumRows()<2)
    {
        itsBailout=true;
        RVec c({1});
        return c;
    }
    
    return SolveC(B); //Solve for the projection coefficients
}

bool SCFIrrepAccelerator_DIIS::CalculateProjections()
{
    // cout << itsIrrep << endl;
    itsBailout=false;
    if (itsBailout=Max(fabs(itsDPrime))==0.0;itsBailout) return itsBailout;
    assert(itsFPrime.GetLimits()==itsDPrime.GetLimits());
    Mat E=CalculateError();
    itsLastEn=FrobeniusNorm(E);
    if (itsBailout=itsLastEn>itsParams.EMax;itsBailout) return itsBailout;
  
    Append(itsFPrime,E,itsLastEn);
    if (itsEs.size()>itsParams.Nproj) Purge1();
    assert(itsEs.size()<=itsParams.Nproj);
    if (itsBailout=itsLastEn< itsParams.EMin;itsBailout) 
        return itsBailout;
    
    SMat B;
    std::tie(B,itsLastSVMin)=BuildPrunedB(itsEs,itsParams.SVTol);
    if (B.GetNumRows()<=2)
    {
        itsBailout=true;
        itsCs=RVec(1);
        itsCs(1)=1.0;
    }
    else
        itsCs=SolveC(B); //Solve for the projection coefficients
    
    return itsBailout;
}


SCFIrrepAccelerator::SMat SCFIrrepAccelerator_DIIS::Project(const RVec& c)
{
    assert(!itsBailout);
    assert(fabs(Sum(c)-1.0)<1e-13); //Check that the constraint worked.
    assert(c.size()>=2);
    assert(c.size()==itsFPrimes.size());
    // Now do the projection for the Fock matrix.
    SMat Fproj;
    index_t i=1;
    for (const auto& f:itsFPrimes) Fproj+=SMat(c(i++)*f);
    return Fproj;

}

#include <algorithm>
void SCFIrrepAccelerator_DIIS::Append(const SMat& FPrime, const Mat& E, double En)
{
    assert(itsEs.size()==itsFPrimes.size());
    assert(itsEns.size()==itsFPrimes.size());
    itsEs.push_back(E);
    itsEns.push_back(En);
    itsFPrimes.push_back(FPrime);   
}

void SCFIrrepAccelerator_DIIS::Purge1()
{
    //auto iter=std::max_element(itsEns.begin(),itsEns.end()); //Find the maximum Error
    auto iter=itsEns.begin(); //just pick the oldest element.
    size_t index=std::distance(itsEns.begin(),iter);
    itsEns.erase(iter);
    itsEs.erase(itsEs.begin()+index);
    itsFPrimes.erase(itsFPrimes.begin()+index);
    
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


SCFAccelerator_DIIS::md_t SCFAccelerator_DIIS::BuildPrunedB(const mv_t& Es,double svmin)
{
    SMat B=SCFIrrepAccelerator_DIIS::BuildB(Es);
    double sv=SCFIrrepAccelerator_DIIS::GetMinSV(B);
    while (sv<svmin && Es.size()>=2) 
    {
        Purge1(); //Must be a member function for this.
        B=SCFIrrepAccelerator_DIIS::BuildB(Es);
        sv=SCFIrrepAccelerator_DIIS::GetMinSV(B);
    }
    return std::make_pair(B,sv);    
}

void SCFAccelerator_DIIS::Purge1()
{
    //auto iter=std::max_element(itsEns.begin(),itsEns.end()); //Find the maximum Error
    auto iter=itsEs.begin(); //just pick the oldest element.
    size_t index=std::distance(itsEs.begin(),iter);
    //itsEns.erase(iter);
    itsEs.erase(itsEs.begin()+index);
    //itsFPrimes.erase(itsFPrimes.begin()+index);
    for (auto k:itsIrreps) k->Purge1();
}

bool SCFAccelerator_DIIS::CalculateProjections()
{
    itsBailout=false;
    switch (itsParams.type)
    {
        case DIISParams::irrep:
        {
            for (auto k:itsIrreps) 
            {
                itsBailout = k->CalculateProjections() || itsBailout;
            }
            break;
        }
        case DIISParams::global:
        {
            itsBailout=false;
            for (auto k:itsIrreps) 
                itsBailout=k->Bailout() || itsBailout;
            if (itsBailout) return itsBailout;

            Mat E;
            for (auto k:itsIrreps) E+=k->CalculateError();
            itsEn=FrobeniusNorm(E);
            if (itsBailout=itsEn>itsParams.EMax;itsBailout) return BailoutChildren();
            itsEs.push_back(E);
            for (auto k:itsIrreps) k->AppendFPrime();
            if (itsBailout=itsEs.size()<2;itsBailout) return BailoutChildren();
            SMat B;
            std::tie(B,itsLastSVMin)=BuildPrunedB(itsEs,itsParams.SVTol);
            if (B.GetNumRows()<=2)
            {
                itsBailout=true;
                itsCs=RVec(1);
                itsCs(1)=1.0;
                return BailoutChildren();
            }
            else
                itsCs=SCFIrrepAccelerator_DIIS::SolveC(B);

            for (auto k:itsIrreps) k->SetProjection(itsCs);
            break;
        }

    }
    return itsBailout;
}
bool SCFAccelerator_DIIS::BailoutChildren()
{
    assert(itsBailout);
        for (auto k:itsIrreps) k->GlobalBailout();
    return itsBailout;
}

void SCFAccelerator_DIIS::ShowLabels(std::ostream& os) const
{
    os << " [F,D]   Nmin   Nmax    SVMin";
}
void SCFAccelerator_DIIS::ShowConvergence(std::ostream& os) const
{
    double EMax=0.0,SVMin=1.0e20;
    size_t NMin=1000,NMax=0;
    bool bailout=false;
    switch (itsParams.type)
    {
        case DIISParams::irrep :
        {
            for (auto i:itsIrreps)
            {
                if (i->GetError()>EMax) EMax=i->GetError();
                if (i->GetSVMin()<SVMin) SVMin=i->GetSVMin();
                if (i->GetNproj()<NMin) NMin=i->GetNproj();
                if (i->GetNproj()>NMax) NMax=i->GetNproj();
                if (i->Bailout()) bailout=true;
            }
            
           
            break;
        }
        case DIISParams::global :
        {
            bailout=itsBailout;
            EMax=itsEn;
            SVMin=itsLastSVMin;
            NMax=NMin=itsEs.size();
            break;
        }
    }
    os << std::scientific << std::setw(7) << std::setprecision(1) << EMax << " ";
    if (!bailout)
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


