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
import qchem.Blaze;

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


void SCFIrrepAcceleratorDIIS::UseFD(const rsmat_t& F, const rsmat_t& DPrime)
{
     itsFPrime=itsLaSolver_blaze->Transform(F); // Fprime = Vd*F*V
    assert(itsFPrime.rows()==DPrime.rows());
    assert(itsFPrime.columns()==DPrime.columns());
    itsDPrime=DPrime;
    itsE=itsFPrime*itsDPrime-itsDPrime*itsFPrime;
    itsEn=norm(itsE);
}

rsmat_t SCFIrrepAcceleratorDIIS::Project()
{
    if (itsCs.size()<2) 
        return itsFPrime;
    else
    {
        double err=fabs(sum(itsCs)-1.0);
        if (err>1e-13)
            cout << "Warning: SCFIrrepAcceleratorDIIS::Project() fabs(Sum(itsCs)-1.0)<>e-13 ." << endl;
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

double SCFAcceleratorDIIS::GetMinSV(const rsmat_t& B)
{
    rvec_t s;
    rmat_t  U,Vt;
    blaze::svd(B,U,s,Vt);
    return s[s.size()-1];
}

rvec_t SCFAcceleratorDIIS::SolveC(const rsmat_t& B) 
{
    size_t N=B.rows();
    rvec_t v(N,0.0);
    v[N-1]=1.0;
    rvec_t C=blaze::solve(B,v);
    return subvector(C,0,N-1);   
}
SCFAcceleratorDIIS::md_t SCFAcceleratorDIIS::BuildB() const
{
    size_t  N=GetNProj();
    rsmat_t B=zero<double>(N+1);
    for (size_t  i=0;i<N;i++)
    {
        B(i,N)=1.0; //B is symmetric so no need to set B(N,i)=1.0
        for (size_t  j=i;j<N;j++)
            for (auto k:itsIrreps) B(i,j)+=k->GetError(i,j);
    }
    // B(N,N)=0.0;  should already be true
    return {B,GetMinSV(B)};    
}
rsmat_t SCFAcceleratorDIIS::BuildPrunedB(double svmin)
{
    md_t B=BuildB(); //Returns a SMat,double struct.
    while (B.sv<svmin && GetNProj()>=2) 
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
    
    rsmat_t B=BuildPrunedB(itsParams.SVTol);
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


