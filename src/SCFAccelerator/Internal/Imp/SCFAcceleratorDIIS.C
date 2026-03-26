// FIle: SCFAcceleratorDIIS.C  Direct Inversion of the Iterative Subspace (DIIS) algorithm
module;
#include <iostream>
#include <cassert>
#include <cmath>
#include <iomanip>
#include "blaze/Math.h" 
module qchem.SCFAccelerator.Internal.SCFAcceleratorDIIS;
import qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull;
import qchem.IrrepBasisSet;
import qchem.Conversions;

using std::cout;
using std::endl;

SCFIrrepAcceleratorDIIS::SCFIrrepAcceleratorDIIS(const DIISParams& p,const LASolver_blaze<double>* lasb,const Irrep_QNs& qns,const rvec_t& cs) 
    : itsParams(p)
    , itsIrrep(qns)
    , itsEn(0.0)
    , itsCs(cs)
    , itsLaSolver_blaze(lasb)
{
    assert(itsLaSolver_blaze);
};
SCFIrrepAcceleratorDIIS::~SCFIrrepAcceleratorDIIS() 
{

};


void SCFIrrepAcceleratorDIIS::UseFD(const smat_t<double>& F, const smat_t<double>& DPrime)
{
    // SMatrix<double> Fb=LASolver_blaze<double>::convert(itsLaSolver_blaze->Transform(LASolver_blaze<double>::convert(F))); // Fprime = Vd*F*V
    // itsFPrime=itsLaSolver->Transform(F); // Fprime = Vd*F*V
    // cout << "Fb-FPrime=" << Fb-itsFPrime << endl;
    itsFPrime=itsLaSolver_blaze->Transform(F); // Fprime = Vd*F*V
    assert(itsFPrime.rows()==DPrime.rows());
    assert(itsFPrime.columns()==DPrime.columns());
    itsDPrime=DPrime;
    itsE=itsFPrime*itsDPrime-itsDPrime*itsFPrime;
    itsEn=norm(itsE);
}

// template <class T> const smat_t<T>& operator+=(smat_t<T>& a, const smat_t<T>& b)
// {
//     if (a.size()==0) 
//     {
//         a=zero<double>(b.rows());
//     }
//     else
//         assert(a.rows()==b.rows());
//     a=a+b;
//     return a;
// }

smat_t<double> SCFIrrepAcceleratorDIIS::Project()
{
    if (itsCs.size()<2) 
        return itsFPrime;
    else
    {
        double err=fabs(sum(itsCs)-1.0);
        if (err>1e-13)
            cout << "Warning: SCFIrrepAcceleratorDIIS::Project() fabs(Sum(itsCs)-1.0)<>e-13 ." << endl;
        // assert(fabs(Sum(itsCs)-1.0)<1e-13); //Check that the constraint worked.
        assert(itsCs.size()==itsFPrimes.size());
        // Now do the projection for the Fock matrix.
        rsmat_t Fproj=zero<double>(itsFPrime.rows()) ;
        size_t  i=0;
        for (const auto& f:itsFPrimes) Fproj+=itsCs[i++]*f;
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
SCFIrrepAccelerator* SCFAcceleratorDIIS::Create(const LASolver_blaze<double>* lasb,const Irrep_QNs& qns, int occ) 
{
    if (occ>0)
    {
        itsIrreps.push_back(new SCFIrrepAcceleratorDIIS(itsParams,lasb,qns,itsCs));
        return itsIrreps.back();
    }
    else
        return new SCFIrrepAcceleratorNull(lasb,qns);
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
    auto [U,s,V]=solver.SolveAll(convert(B));
    size_t N=s.GetNumRows();
    return s(N,N);
}

rvec_t SCFAcceleratorDIIS::SolveC(const SMat& B) 
{
    static oml::LapackLinearSolver<double> solver;
    size_t N=B.rows();
    RVec v(N);
    Fill(v,0.0); 
    v(N)=1.0;
    RVec C=solver.Solve(convert(B),v);
    // std:: cout << B << C << v << del <<std::endl;
    // std:: cout << "del,[F,D] = " << sqrt(del*del) << " " << C(N) << std::endl;
    return convert(C.SubVector(N-1));   
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
    return {convert(B),GetMinSV(convert(B))};    
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
    blaze::clear(itsCs);
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
    if (B.rows()<=2) return false;
                
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


