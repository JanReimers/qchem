// FIle: SCFAccelerator_DIIS.C  Direct Inversion of the Iterative Subspace (DIIS) algorithm


#include <iostream>
#include <cassert>
#include <cmath>
#include <iomanip>
#include "SCFAccelerator_DIIS.H"
#include <LASolver/LASolver.H>
#include <BasisSet/Irrep_BS.H>
#include "SCFAccelerator_Null.H"

using std::cout;
using std::endl;

SCFIrrepAccelerator_DIIS::SCFIrrepAccelerator_DIIS(const DIISParams& p,const LASolver<double>* las,const Irrep_QNs& qns,const RVec& cs) 
    : itsParams(p)
    , itsIrrep(qns)
    , itsEn(0.0)
    , itsCs(cs)
    , itsLaSolver(las)
{
    assert(itsLaSolver);
};
SCFIrrepAccelerator_DIIS::~SCFIrrepAccelerator_DIIS() 
{

};


void SCFIrrepAccelerator_DIIS::UseFD(const SMat& F, const SMat& DPrime)
{
    itsFPrime=itsLaSolver->Transform(F); // Fprime = Vd*F*V
    assert(itsFPrime.GetLimits()==DPrime.GetLimits());
    itsDPrime=DPrime;
    itsE=Mat(itsFPrime*itsDPrime-itsDPrime*itsFPrime);
    itsEn=FrobeniusNorm(itsE);
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
// template <class T> const Matrix<T>& operator+=(Matrix<T>& a, const Matrix<T>& b)
// {
//     if (a.size()==0) 
//     {
//         a.SetLimits(b.GetLimits());
//         Fill(a,0.0);
//     }
//     else
//         assert(a.GetLimits()==b.GetLimits());
//     return ArrayAdd(a,b);
// }

SCFIrrepAccelerator::SMat SCFIrrepAccelerator_DIIS::Project()
{
    if (itsCs.size()<2) 
        return itsFPrime;
    else
    {
        assert(fabs(Sum(itsCs)-1.0)<1e-13); //Check that the constraint worked.
        assert(itsCs.size()==itsFPrimes.size());
        // Now do the projection for the Fock matrix.
        SMat Fproj;
        index_t i=1;
        for (const auto& f:itsFPrimes) Fproj+=SMat(itsCs(i++)*f);
        return Fproj;
    }
}

void SCFIrrepAccelerator_DIIS::Append1()
{
    assert(itsEs.size()==itsFPrimes.size());
    assert(itsEns.size()==itsFPrimes.size());
    itsEs     .push_back(itsE);
    itsEns    .push_back(itsEn);
    itsFPrimes.push_back(itsFPrime);   
}
void SCFIrrepAccelerator_DIIS::Purge1()
{
    assert(itsEs.size()==itsFPrimes.size());
    assert(itsEns.size()==itsFPrimes.size());
    itsEns    .pop_front();
    itsEs     .pop_front();
    itsFPrimes.pop_front();    
}

//----------------------------------------------------------------------------------------------------------------------------
//
// Non irrep code
//


SCFAccelerator_DIIS::SCFAccelerator_DIIS(const DIISParams& p) 
: itsParams(p)
{};

SCFAccelerator_DIIS::~SCFAccelerator_DIIS() {};
SCFIrrepAccelerator* SCFAccelerator_DIIS::Create(const LASolver<double>* las,const Irrep_QNs& qns, int occ) 
{
    if (occ>0)
    {
        itsIrreps.push_back(new SCFIrrepAccelerator_DIIS(itsParams,las,qns,itsCs));
        return itsIrreps.back();
    }
    else
        return new SCFIrrepAccelerator__Null(las,qns);
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
    static oml::LapackLinearSolver<double> solver;
    size_t N=B.GetNumRows();
    RVec v(N);
    Fill(v,0.0); 
    v(N)=1.0;
    RVec C=solver.Solve(B,v);
    // std:: cout << B << C << v << del <<std::endl;
    // std:: cout << "del,[F,D] = " << sqrt(del*del) << " " << C(N) << std::endl;
    return C.SubVector(N-1);   
}
SCFAccelerator_DIIS::md_t SCFAccelerator_DIIS::BuildB() const
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
    return {B,GetMinSV(B)};    
}
SCFAccelerator_DIIS::SMat SCFAccelerator_DIIS::BuildPrunedB(double svmin)
{
    md_t B=BuildB(); //Returns a SMat,double struct.
    while (B.sv<svmin &&GetNProj()>=2) 
    {
        Purge1(); //Must be a member function for this.
        B=BuildB();
    }
    itsLastSVMin=B.sv;
    return B.B;    
}
size_t SCFAccelerator_DIIS::Purge1()
{
    
    for (auto k:itsIrreps) k->Purge1();
    return GetNProj();
}
size_t SCFAccelerator_DIIS::Append1()
{
    
    for (auto k:itsIrreps) k->Append1();
    return GetNProj();
}

bool SCFAccelerator_DIIS::CalculateProjections()
{
    itsCs.SetLimits(0);
    itsEn=0.0;
    for (auto k:itsIrreps) 
    {
        double Enk=k->GetError();
        if (Enk==0.0) return false;
        itsEn+=Enk*Enk;
    }
    itsEn=sqrt(itsEn);
    // cout << "itsEn=" << itsEn << endl;
    if (itsEn>itsParams.EMax) return false;
    
    if (Append1()>itsParams.Nproj) Purge1();
    assert(GetNProj()<=itsParams.Nproj);
    if (GetNProj()<2) return false;
    
    SMat B=BuildPrunedB(itsParams.SVTol);
    if (B.GetNumRows()<=2) return false;
                
    itsCs=SCFAccelerator_DIIS::SolveC(B); //Irreps have a refeence to this in order to the the projections.
   
    return true;
}


void SCFAccelerator_DIIS::ShowLabels(std::ostream& os) const
{
    os << " [F,D]   Nproj    SVMin";
}
void SCFAccelerator_DIIS::ShowConvergence(std::ostream& os) const
{
    os << std::scientific << std::setw(7) << std::setprecision(1) << itsEn << " ";
    if (HasProjection())
    {
        os << std::setw(3) << GetNProj() << "    ";
        os << std::scientific << std::setw(7) << std::setprecision(1) << itsLastSVMin << "  ";
    }
        else
        os << "                ";
}
double SCFAccelerator_DIIS::GetError() const
{
      
    return itsEn;
}


