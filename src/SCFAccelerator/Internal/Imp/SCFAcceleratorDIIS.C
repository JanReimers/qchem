// FIle: SCFAcceleratorDIIS.C  Direct Inversion of the Iterative Subspace (DIIS) algorithm
module;
#include <iostream>
#include <cassert>
#include <cmath>
#include <iomanip>
module qchem.SCFAccelerator.Internal.SCFAcceleratorDIIS;
import qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull;
import qchem.Irrep_BS;

using std::cout;
using std::endl;

SCFIrrepAcceleratorDIIS::SCFIrrepAcceleratorDIIS(const DIISParams& p,const LASolver<double>* las,const Irrep_QNs& qns,const RVec& cs) 
    : itsParams(p)
    , itsIrrep(qns)
    , itsEn(0.0)
    , itsCs(cs)
    , itsLaSolver(las)
{
    assert(itsLaSolver);
};
SCFIrrepAcceleratorDIIS::~SCFIrrepAcceleratorDIIS() 
{

};


void SCFIrrepAcceleratorDIIS::UseFD(const SMatrix<double>& F, const SMatrix<double>& DPrime)
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

SMatrix<double> SCFIrrepAcceleratorDIIS::Project()
{
    if (itsCs.size()<2) 
        return itsFPrime;
    else
    {
        double err=fabs(Sum(itsCs)-1.0);
        if (err>1e-13)
            cout << "Warning: SCFIrrepAcceleratorDIIS::Project() fabs(Sum(itsCs)-1.0)<>e-13 ." << endl;
        // assert(fabs(Sum(itsCs)-1.0)<1e-13); //Check that the constraint worked.
        assert(itsCs.size()==itsFPrimes.size());
        // Now do the projection for the Fock matrix.
        SMatrix<double> Fproj;
        size_t  i=1;
        for (const auto& f:itsFPrimes) Fproj+=SMatrix<double>(itsCs(i++)*f);
        return Fproj;
    }
}

void SCFIrrepAcceleratorDIIS::Append1()
{
    assert(itsEs.size()==itsFPrimes.size());
    assert(itsEns.size()==itsFPrimes.size());
    itsEs     .push_back(itsE);
    itsEns    .push_back(itsEn);
    itsFPrimes.push_back(itsFPrime);   
}
void SCFIrrepAcceleratorDIIS::Purge1()
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


SCFAcceleratorDIIS::SCFAcceleratorDIIS(const DIISParams& p) 
: itsParams(p)
{};

SCFAcceleratorDIIS::~SCFAcceleratorDIIS() {};
SCFIrrepAccelerator* SCFAcceleratorDIIS::Create(const LASolver<double>* las,const Irrep_QNs& qns, int occ) 
{
    if (occ>0)
    {
        itsIrreps.push_back(new SCFIrrepAcceleratorDIIS(itsParams,las,qns,itsCs));
        return itsIrreps.back();
    }
    else
        return new SCFIrrepAcceleratorNull(las,qns);
}

size_t SCFAcceleratorDIIS::GetNProj() const
{
    size_t N=itsIrreps[0]->GetNproj();
#ifdef DEBUG
    for (auto k:itsIrreps) assert(N==k->GetNproj());
#endif
    return N;
}

double SCFAcceleratorDIIS::GetMinSV(const SMat& B)
{
    static oml::LapackSVDSolver<double> solver;
    auto [U,s,V]=solver.SolveAll(B);
    size_t N=s.GetNumRows();
    return s(N,N);
}

RVec SCFAcceleratorDIIS::SolveC(const SMat& B) 
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
SCFAcceleratorDIIS::md_t SCFAcceleratorDIIS::BuildB() const
{
    size_t  N=GetNProj()+1;
    SMatrix<double> B(N);
    Fill(B,0.0);
    for (size_t  i=1;i<N;i++)
    {
        B(i,N)=1.0; //B is symmetric so no need to set B(N,i)=1.0
        for (size_t  j=i;j<N;j++)
            for (auto k:itsIrreps) B(i,j)+=k->GetError(i-1,j-1);
    }
    // B(N,N)=0.0;  should already be true
    return {B,GetMinSV(B)};    
}
SCFAcceleratorDIIS::SMat SCFAcceleratorDIIS::BuildPrunedB(double svmin)
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
size_t SCFAcceleratorDIIS::Purge1()
{
    
    for (auto k:itsIrreps) k->Purge1();
    return GetNProj();
}
size_t SCFAcceleratorDIIS::Append1()
{
    
    for (auto k:itsIrreps) k->Append1();
    return GetNProj();
}

bool SCFAcceleratorDIIS::CalculateProjections()
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
    
    SMatrix<double> B=BuildPrunedB(itsParams.SVTol);
    if (B.GetNumRows()<=2) return false;
                
    itsCs=SCFAcceleratorDIIS::SolveC(B); //Irreps have a refeence to this in order to the the projections.
   
    return true;
}


void SCFAcceleratorDIIS::ShowLabels(std::ostream& os) const
{
    os << " [F,D]   Nproj    SVMin";
}
void SCFAcceleratorDIIS::ShowConvergence(std::ostream& os) const
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
double SCFAcceleratorDIIS::GetError() const
{
      
    return itsEn;
}


