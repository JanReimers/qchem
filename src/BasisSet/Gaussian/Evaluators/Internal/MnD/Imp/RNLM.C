// File: BasisSet/Gaussian/Evaluators/Internal/MnD/Imp/RNLM.C  Manage the R_NLM Hermite Coulomb auxiliaries.
module;
#include <iostream>
#include <cassert>
#include <vector>
module qchem.BasisSet.Gaussian.Evaluators.Internal.MnD.RNLM;
import qchem.BasisSet.Gaussian.Evaluators.Internal.Boys.AuxillaryFJ;
import qchem.IntPow;

#define MAX 16

namespace qchem::BasisSet::Gaussian::Evaluators::Internal::MnD
{
using Boys::AuxillaryFJ;


RNLM::RNLM() : itsData() {};

RNLM::RNLM(int Max, double α, const rvec3_t& dR)
    : itsLMax(Max)
    , itsData()
{
    Build(Max, α, dR);
}

void RNLM::Build(int Max, double α, const rvec3_t& dR)
{
    itsLMax=Max;
    // The recursion's layers R_j, j = 0..MAX: a per-THREAD scratch.  It was function-static, shared by every caller, so two RNLMs built
    // concurrently overwrote each other's layers (M_BoysTable.RNLMsBuiltConcurrentlyAreTheSerialOnes: ~half of 4000 differed at 4 workers).
    // No production path builds one in a parallel region today (OpenWork §3); this keeps it that way by construction, not by convention.
    thread_local Triangle3D theRjs[MAX+1];
    thread_local bool       theRjsInitialized = false;
    if(!theRjsInitialized)
    {
        assert(MAX==AuxillaryFJ::thejMax);
        int GlobalMax=AuxillaryFJ::thejMax;
        for (int j=0; j<=GlobalMax; j++) theRjs[j] = Triangle3D(GlobalMax-j);
        theRjsInitialized=true;
    }

    AuxillaryFJ FjCalculator(Max); // throws on Max > thejMax, before theRjs[thejMax-Max] below is indexed
    double T=α*dR*dR;
    thread_local rvec_t Fj;   // scratch: GetFjAt overwrites every entry, so only its size matters
    if (Fj.size()!=size_t(Max+1)) Fj.resize(Max+1,false);
    FjCalculator.GetFjAt(T,Fj);

//  cout << "T=" << T << " Max = " << Max << "  Fj=" << Fj << std::endl;
    Triangle3D* Rjs=&theRjs[AuxillaryFJ::thejMax-Max];
//  for (int j=0;j<=Max;j++) Rjs.push_back(new Triangle(Max-j));

    HermiteCoulomb(Max, α, dR, Fj, [Rjs](int j, int N, int L, int M) -> double& {return Rjs[j](N,L,M);});
    itsData = Rjs[0];
}

void RNLM::Add(const RNLM& theR, double theScale)
{
    itsData.Add(theR.itsData,theScale);
}


} //namespace
