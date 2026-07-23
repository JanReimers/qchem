// FIle: SCFAcceleratorDIIS.C  Direct Inversion of the Iterative Subspace (DIIS) algorithm
module;
#include <iostream>
#include <cassert>
#include <iomanip>
#include <complex>   // std::real
module qchem.SCFAccelerator.Internal.SCFAcceleratorDIIS;
import qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull;
import qchem.Blaze;
import qchem.Math;
import qchem.Math.DIIS;   // the shared Pulay bordered-solve (Bordered/MinSV/Coefficients) -- serves Fock-DIIS too

namespace qchem::SCFAccelerators
{
using std::cout;
using std::endl;

template <class T> tSCFIrrepAcceleratorDIIS<T>::
tSCFIrrepAcceleratorDIIS(const DIISParams& p,const LASolver<T>* lasb,const Irrep& qns,const rvec_t& cs)
    : itsParams(p)
    , itsIrrep(qns)
    , itsEn(0.0)
    , itsCs(cs)
    , itsLASolver(lasb)
{
    assert(itsLASolver);
};
template <class T> tSCFIrrepAcceleratorDIIS<T>::~tSCFIrrepAcceleratorDIIS()
{

};


template <class T> void tSCFIrrepAcceleratorDIIS<T>::UseFD(const hmat_t<T>& F, const hmat_t<T>& DPrime)
{
    itsFPrime=itsLASolver->Transform(F); // Fprime = Vd*F*V
    assert(itsFPrime.rows()==DPrime.rows());
    assert(itsFPrime.columns()==DPrime.columns());
    itsDPrime=DPrime;
    itsE=itsFPrime*itsDPrime-itsDPrime*itsFPrime; // [F',D'] (anti-Hermitian)
    itsEn=std::real(blazem::norm(itsE));
}

template <class T> typename LASolver<T>::UUd_t tSCFIrrepAcceleratorDIIS<T>::NextOrbitals()
{
    // Fock-extrapolation accelerator: extrapolate F' then diagonalize it.
    return itsLASolver->SolveOrtho(Project());
}

template <class T> hmat_t<T> tSCFIrrepAcceleratorDIIS<T>::Project()
{
    if (itsCs.size()<2)
        return itsFPrime;
    else
    {
        double err=fabs(blazem::sum(itsCs)-1.0);
        if (err>1e-13)
            cout << "Warning: tSCFIrrepAcceleratorDIIS::Project() fabs(Sum(itsCs)-1.0)<>e-13 ." << endl;
        assert(itsCs.size()==itsFPrimes.size());
        // Now do the projection for the Fock matrix (real coefficients * Hermitian F' = Hermitian).
        hmat_t<T> Fproj=blazem::zeroH<T>(itsFPrime.rows()) ;
        size_t  i=0;
        for (const auto& f:itsFPrimes) Fproj+=itsCs[i++]*f;
        return Fproj;
    }
}

template <class T> void tSCFIrrepAcceleratorDIIS<T>::Append1()
{
    assert(itsEs.size()==itsFPrimes.size());
    assert(itsEns.size()==itsFPrimes.size());
    itsEs     .push_back(itsE);
    itsEns    .push_back(itsEn);
    itsFPrimes.push_back(itsFPrime);
}
template <class T> void tSCFIrrepAcceleratorDIIS<T>::Purge1()
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


template <class T> tSCFAcceleratorDIIS<T>::tSCFAcceleratorDIIS(const DIISParams& p)
: itsParams(p)
{};

template <class T> tSCFAcceleratorDIIS<T>::~tSCFAcceleratorDIIS() {};
template <class T> tSCFIrrepAccelerator<T>* tSCFAcceleratorDIIS<T>::Create(const LASolver<T>* lasb,const Irrep& qns, int occ)
{
    if (occ>0)
    {
        itsIrreps.push_back(new tSCFIrrepAcceleratorDIIS<T>(itsParams,lasb,qns,itsCs));
        return itsIrreps.back();
    }
    else
        return new tSCFIrrepAcceleratorNull<T>(lasb,qns);
}

template <class T> size_t tSCFAcceleratorDIIS<T>::GetNProj() const
{
    size_t N=itsIrreps[0]->GetNproj();
#ifdef DEBUG
    for (auto k:itsIrreps) assert(N==k->GetNproj());
#endif
    return N;
}

template <class T> typename tSCFAcceleratorDIIS<T>::md_t tSCFAcceleratorDIIS<T>::BuildB() const
{
    size_t  N=GetNProj();
    rsmat_t Braw=blazem::zero<double>(N);           // raw error-overlap Bᵢⱼ = Σ_irreps ⟨Eᵢ,Eⱼ⟩ (upper-tri)
    for (size_t  i=0;i<N;i++)
        for (size_t  j=i;j<N;j++)
            for (auto k:itsIrreps) Braw(i,j)+=k->GetError(i,j);
    rsmat_t B=qchem::Math::DIIS::Bordered(Braw);    // the (N+1) bordered Pulay system
    return {B,qchem::Math::DIIS::MinSV(B)};
}
template <class T> rsmat_t tSCFAcceleratorDIIS<T>::BuildPrunedB(double svmin)
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
template <class T> size_t tSCFAcceleratorDIIS<T>::Purge1()
{

    for (auto k:itsIrreps) k->Purge1();
    return GetNProj();
}
template <class T> size_t tSCFAcceleratorDIIS<T>::Append1()
{

    for (auto k:itsIrreps) k->Append1();
    return GetNProj();
}

template <class T> bool tSCFAcceleratorDIIS<T>::CalculateProjections()
{
    blazem::clear(itsCs);
    itsEn=0.0;
    bailoutReason="            ";
    // A zero error guards two distinct cases: (1) the zero-initial-density first iterations (every
    // channel's [F',D']==0 until it has a density -- atoms rely on this), and (2) a permanently-empty
    // irrep (e.g. A2 for H2O).  We must still BAIL on (1) but SKIP (2), else symmetric molecules
    // never extrapolate.  Discriminator: once any irrep has shown a nonzero error we are "seeded"
    // (past case 1), so from then on an exactly-zero error means a permanently-empty irrep.
    for (auto k:itsIrreps) if (k->GetError()!=0.0) itsSeeded=true;
    for (auto k:itsIrreps)
    {
        double Enk=k->GetError();
        if (Enk==0.0)
        {
            if (itsSeeded) continue;            // permanently-empty irrep: contributes nothing, skip it
            bailoutReason="Enk==0.0    ";       // not yet seeded (zero-density start): must bail
            return false;
        }
        itsEn+=Enk*Enk;
    }
    itsEn=sqrt(itsEn);
    // cout << "itsEn=" << itsEn << endl;
    if (itsEn>itsParams.EMax)
    {
        bailoutReason="En>EMax     ";
        itsStuckCount=0;        //early (pre-DIIS) phase, not "out of steam"
        return false;
    }

    if (Append1()>itsParams.Nproj) Purge1();
    assert(GetNProj()<=itsParams.Nproj);
    if (GetNProj()<2)
    {
        bailoutReason="Nproj<2     ";
        itsStuckCount++;
        return false;
    }

    rsmat_t B=BuildPrunedB(itsParams.SVTol);
    if (B.rows()<=2)
    {
        bailoutReason="B.rows()<=2 ";
        itsStuckCount++;
        return false;
    }

    itsCs=qchem::Math::DIIS::Coefficients(B); //Irreps have a reference to this in order to do the projections.
    itsStuckCount=0;
    return true;
}


template <class T> void tSCFAcceleratorDIIS<T>::ShowLabels(std::ostream& os) const
{
    os << " [F,D]   Nproj    SVMin   Bail   ";
}
template <class T> void tSCFAcceleratorDIIS<T>::ShowConvergence(std::ostream& os) const
{
    os << std::scientific << std::setw(7) << std::setprecision(1) << itsEn << " ";
    if (HasProjection())
    {
        os << std::setw(3) << GetNProj() << "    ";
        os << std::scientific << std::setw(7) << std::setprecision(1) << itsLastSVMin << "  ";
    }
        else
        os << "                ";
    os << bailoutReason;
}
template <class T> double tSCFAcceleratorDIIS<T>::GetError() const
{

    return itsEn;
}

template class tSCFIrrepAcceleratorDIIS<double>;
template class tSCFIrrepAcceleratorDIIS<dcmplx>;
template class tSCFAcceleratorDIIS<double>;
template class tSCFAcceleratorDIIS<dcmplx>;

} //namespace
