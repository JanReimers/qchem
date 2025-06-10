// FIle: SCFAccelerator_DIIS.C  Direct Inversion of the Iterative Subspace (DIIS) algorithm


#include "Imp/SCF/SCFAccelerator_DIIS.H"
#include <LASolver.H>
#include "oml/numeric/LapackLinearSolver.H"
#include "oml/numeric/LapackSVDSolver.H"

using std::cout;
using std::endl;

SCFIrrepAccelerator_DIIS::SCFIrrepAccelerator_DIIS(const DIISParams& p,const Irrep_QNs& qns) 
    : itsParams(p)
    , itsIrrep(qns)
    , itsLastEn(0.0)
    , itsCs(1)
{
    itsCs(1)=1.0;
};
SCFIrrepAccelerator_DIIS::~SCFIrrepAccelerator_DIIS() 
{

};
void SCFIrrepAccelerator_DIIS::Init(const LASolver<double>* las)
{
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

SCFIrrepAccelerator::Mat SCFIrrepAccelerator_DIIS::CalculateError() 
{
    itsLastE=Mat(itsFPrime*itsDPrime-itsDPrime*itsFPrime);
    return itsLastE;
}

template <class T> const SMatrix<T>& operator+=(SMatrix<T>& a, const SMatrix<T>& b)
{
    if (a.size()==0) 
    {
        a.SetLimits(b.GetLimits());
        Fill(a,0.0);
    }
    else
        assert(a.GetLimits()==b.GetLimits());
    return ArrayAdd(a,b);
}
template <class T> const Matrix<T>& operator+=(Matrix<T>& a, const Matrix<T>& b)
{
    if (a.size()==0) 
    {
        a.SetLimits(b.GetLimits());
        Fill(a,0.0);
    }
    else
        assert(a.GetLimits()==b.GetLimits());
    return ArrayAdd(a,b);
}

SCFIrrepAccelerator::SMat SCFIrrepAccelerator_DIIS::Project()
{
    if (itsBailout) 
        return itsFPrime;
    else
    {
        assert(fabs(Sum(itsCs)-1.0)<1e-13); //Check that the constraint worked.
        assert(itsCs.size()>=2);
        assert(itsCs.size()==itsFPrimes.size());
        // Now do the projection for the Fock matrix.
        SMat Fproj;
        index_t i=1;
        for (const auto& f:itsFPrimes) Fproj+=SMat(itsCs(i++)*f);
        return Fproj;
        // return Project(itsCs);
    }
}

#include "Imp/Containers/stl_io.h"

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

//----------------------------------------------------------------------------------------------------------------------------
//
// Non irrep code
//
#include <Irrep_BS.H>
#include "oml/diagonalmatrix.h"


SCFAccelerator_DIIS::SCFAccelerator_DIIS(const DIISParams& p) 
: itsParams(p)
{};

SCFAccelerator_DIIS::~SCFAccelerator_DIIS() {};
SCFIrrepAccelerator* SCFAccelerator_DIIS::Create(const Irrep_QNs& qns) 
{
    itsIrreps.push_back(new SCFIrrepAccelerator_DIIS(itsParams,qns));
    return itsIrreps.back();;
}

size_t SCFAccelerator_DIIS::GetNProj() const
{
    size_t N=itsIrreps[0]->GetNproj();
#ifdef DEBUG
    for (auto k:itsIrreps) assert(N==k->GetNproj());
#endif
    return N;
}

double SCFAccelerator_DIIS::GetMinSV(const SMat& B)
{
    static oml::LapackSVDSolver<double> solver;
    auto [U,s,V]=solver.SolveAll(B);
    size_t N=s.GetNumRows();
    return s(N,N);
}
SCFAccelerator_DIIS::RVec SCFAccelerator_DIIS::SolveC(const SMat& B) 
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


SCFAccelerator_DIIS::SMat SCFAccelerator_DIIS::BuildB() const
{
    index_t N=GetNProj()+1;
    SMat B(N);
    Fill(B,0.0);
    for (index_t i=1;i<N;i++)
    {
        B(i,N)=1.0; //B is symmetric so no need to set B(N,i)=1.0
        for (index_t j=i;j<N;j++)
            for (auto k:itsIrreps) B(i,j)+=k->GetError(i-1,j-1);
    }
    // B(N,N)=0.0;  should already be true
    return B;
}

SCFAccelerator_DIIS::md_t SCFAccelerator_DIIS::BuildPrunedB(double svmin)
{
    SMat B=BuildB();
    double sv=SCFAccelerator_DIIS::GetMinSV(B);
    while (sv<svmin &&GetNProj()>=2) 
    {
        Purge1(); //Must be a member function for this.
        B=BuildB();
        sv=SCFAccelerator_DIIS::GetMinSV(B);
    }
    return std::make_pair(B,sv);    
}

void SCFAccelerator_DIIS::Purge1()
{
    
    for (auto k:itsIrreps) k->Purge1();
}

bool SCFAccelerator_DIIS::CalculateProjections()
{
    itsBailout=false;

    itsEn=0.0;
    for (auto k:itsIrreps) 
    {
        double Enk=FrobeniusNorm(k->CalculateError());
        if (Enk==0.0) return BailoutChildren();
        itsEn+=Enk*Enk;
    }
    itsEn=sqrt(itsEn);
    // cout << "itsEn=" << itsEn << endl;
    if (itsEn>itsParams.EMax) return BailoutChildren();
    
    for (auto k:itsIrreps) k->AppendFPrime();
    if (GetNProj()>itsParams.Nproj) Purge1();
    assert(GetNProj()<=itsParams.Nproj);
    if (GetNProj()<2) return BailoutChildren();
    
    SMat B;
    std::tie(B,itsLastSVMin)=BuildPrunedB(itsParams.SVTol);
    if (itsBailout=B.GetNumRows()<=2;itsBailout) return BailoutChildren();
                
    itsCs=SCFAccelerator_DIIS::SolveC(B);
    for (auto k:itsIrreps) k->SetProjection(itsCs);

    return itsBailout;
}
bool SCFAccelerator_DIIS::BailoutChildren()
{
    itsBailout=true;
        for (auto k:itsIrreps) k->GlobalBailout();
    return itsBailout;
}

void SCFAccelerator_DIIS::ShowLabels(std::ostream& os) const
{
    os << " [F,D]   Nproj    SVMin";
}
void SCFAccelerator_DIIS::ShowConvergence(std::ostream& os) const
{
    os << std::scientific << std::setw(7) << std::setprecision(1) << itsEn << " ";
    if (!itsBailout)
    {
        os << std::setw(3) << GetNProj() << "    ";
        os << std::scientific << std::setw(7) << std::setprecision(1) << itsLastSVMin << "  ";
    }
        else
        os << "                 ";
}

double SCFAccelerator_DIIS::GetError() const
{
      
    return itsEn;
}


