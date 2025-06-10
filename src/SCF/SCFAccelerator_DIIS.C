// FIle: SCFAccelerator_DIIS.C  Direct Inversion of the Iterative Subspace (DIIS) algorithm


#include "Imp/SCF/SCFAccelerator_DIIS.H"
#include <LASolver.H>
#include "oml/numeric/LapackLinearSolver.H"
#include "oml/numeric/LapackSVDSolver.H"

using std::cout;
using std::endl;


#include "oml/diagonalmatrix.h"

SCFIrrepAccelerator_DIIS::SMat SCFIrrepAccelerator_DIIS::BuildRawB(const std::vector< Matrix<double>>& Es)
{
    index_t N=Es.size()+1;
    SMat B(N,N);
    for (index_t i:B.rows())
        if (i<N)
            for (index_t j:B.cols(i))
                if (j<N)
                    B(i,j)=Dot(Es[i-1],Es[j-1]);
                
    return B;
}
SCFIrrepAccelerator_DIIS::SMat& SCFIrrepAccelerator_DIIS::AddBEdges(SMat& B)
{
    index_t N=B.GetNumRows();
    for (index_t i:B.rows())
        if (i<N)
            B(i,N)=1;
        else
            B(N,N)=0.0; 
    return B;
}
SCFIrrepAccelerator_DIIS::SMat SCFIrrepAccelerator_DIIS::BuildB(const std::vector< Matrix<double>>& Es)
{
    SMat B(BuildRawB(Es));
    return AddBEdges(B);
   
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



SCFIrrepAccelerator_DIIS::SCFIrrepAccelerator_DIIS(const DIISParams& p,const Irrep_QNs& qns) 
    : itsParams(p)
    , itsIrrep(qns)
    , itsLastEn(0.0)
    , itsLastSVMin(1e20)
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
    if (itsBailout=itsLastEn< itsParams.EMin;itsBailout) return itsBailout;
    
    SMat B;
    std::tie(B,itsLastSVMin)=BuildPrunedB(itsEs,itsParams.SVTol);
    if (itsBailout=B.GetNumRows()<=2;itsBailout) return itsBailout;
    
    itsCs=SolveC(B); //Solve for the projection coefficients
    
    return itsBailout;
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
SCFIrrepAccelerator* SCFAccelerator_DIIS::Create(const Irrep_QNs& qns) 
{
    itsIrreps.push_back(new SCFIrrepAccelerator_DIIS(itsParams,qns));
    return itsIrreps.back();;
}


SCFAccelerator_DIIS::md_t SCFAccelerator_DIIS::BuildPrunedB(double svmin)
{
    SMat B;
    for (auto k:itsIrreps) B+=k->BuildRawB();
    SCFIrrepAccelerator_DIIS::AddBEdges(B);
    
    double sv=SCFIrrepAccelerator_DIIS::GetMinSV(B);
    while (sv<svmin &&itsN>=2) 
    {
        Purge1(); //Must be a member function for this.
        B.SetLimits(0);
        for (auto k:itsIrreps) B+=k->BuildRawB();
        SCFIrrepAccelerator_DIIS::AddBEdges(B);
        sv=SCFIrrepAccelerator_DIIS::GetMinSV(B);
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
    for (auto k:itsIrreps) 
        itsBailout=k->Bailout() || itsBailout;
    if (itsBailout) return BailoutChildren();

    itsEn=0.0;
    for (auto k:itsIrreps) 
    {
        double Enk=FrobeniusNorm(k->CalculateError());
        itsEn+=Enk*Enk;
    }
    itsEn=sqrt(itsEn);
    // cout << "itsEn=" << itsEn << endl;
    if (itsBailout=itsEn>itsParams.EMax;itsBailout) return BailoutChildren();
    
    for (auto k:itsIrreps) k->AppendFPrime();
    itsN= itsIrreps[0]->GetNproj();
    if (itsN>itsParams.Nproj) Purge1();
    itsN= itsIrreps[0]->GetNproj();
    assert(itsN<=itsParams.Nproj);
    if (itsBailout=itsN<2;itsBailout) return BailoutChildren();
    
    SMat B;
    std::tie(B,itsLastSVMin)=BuildPrunedB(itsParams.SVTol);
    if (itsBailout=B.GetNumRows()<=2;itsBailout) return BailoutChildren();
                
    itsCs=SCFIrrepAccelerator_DIIS::SolveC(B);
    for (auto k:itsIrreps) k->SetProjection(itsCs);

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
    os << " [F,D]   Nproj    SVMin";
}
void SCFAccelerator_DIIS::ShowConvergence(std::ostream& os) const
{
    os << std::scientific << std::setw(7) << std::setprecision(1) << itsEn << " ";
    if (!itsBailout)
    {
        os << std::setw(3) << itsN << "    ";
        os << std::scientific << std::setw(7) << std::setprecision(1) << itsLastSVMin << "  ";
    }
        else
        os << "                 ";
}

double SCFAccelerator_DIIS::GetError() const
{
      
    return itsEn;
}


