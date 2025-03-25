// File: FittedCDImplementation.C  General implementation using a density matrix.



#include "Imp/ChargeDensity/FittedCD.H"
#include <Mesh.H>
#include <Irrep_BS.H>
#include "oml/smatrix.h"
#include "oml/imp/binio.h"
#include <cmath>
#include <cassert>
#include <stdlib.h>

//------------------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> FittedCDImp<T>::FittedCDImp()
    : IntegralConstrainedFF<double>()
    , itsExactRep(0)
    , itsTotalCharge(0)
{}; // No UT coverage

template <class T> FittedCDImp<T>::FittedCDImp(bs_t& bs, mesh_t& m)
    : IntegralConstrainedFF<double>(bs,m) //Use repulsion overlap for fitting
    , itsExactRep(0)
    , itsTotalCharge(0)
{}; // No UT coverage

template <class T> FittedCDImp<T>::FittedCDImp(bs_t& bs, mesh_t& m, double totalCharge)
    : IntegralConstrainedFF<double>(bs,m) //Use repulsion overlap for fitting
    , itsExactRep(0)
    , itsTotalCharge(totalCharge)
{
    FittedFunctionImp<double>::ReScale(itsTotalCharge);
    assert(itsTotalCharge>0);
    assert(abs(itsTotalCharge-FitGetCharge())<1e-10);
};


template <class T> double FittedCDImp<T>::DoFit(const DensityFFClient& ffc)
{
    const Exact_CD* cd=dynamic_cast<const Exact_CD*>(&ffc);
    if (!cd)
    {
        std::cerr << "FittedCDImplementation<T>::DoFit could not cast to charge density" << std::endl;
        exit(-1);
    }
    itsExactRep=cd;
    if  (itsTotalCharge==0) itsTotalCharge=itsExactRep->GetTotalCharge();
    return ConstrainedFF<double>::DoFit(ffc);
}
template <class T> double FittedCDImp<T>::DoFit(const Exact_CD& cd)
{
    itsExactRep=&cd;
    if  (itsTotalCharge==0) itsTotalCharge=itsExactRep->GetTotalCharge();
    return ConstrainedFF<double>::DoFit(cd);
}
//-----------------------------------------------------------------------------
//
//  Totale energy terms for a charge density.
//

template <class T> Exact_CD::SMat FittedCDImp<T>::GetRepulsion(const TOrbital_DFT_IBS<double>* bs) const
{
    assert(bs);
    const std::vector<SMat>& repulsions=bs->Repulsion3C(*itsBasisSet);
    int n=bs->GetNumFunctions();
    SMat J(n,n);
    Fill(J,0.0);
    size_t i=0;
    for (auto c:itsFitCoeff) J+=SMat(c*repulsions[i++]);
    assert(!isnan(J));
    return J;
}

// template <class T> double FittedCDImp<T>::GetEnergy(const HamiltonianTerm* v) const
// {
//     // No UT coverage
//     assert(itsExactRep);
//     return itsExactRep->GetEnergy(v);
// }

template <class T> double FittedCDImp<T>::GetSelfRepulsion() const
{
    return 0.5 * FittedFunctionImp<T>::FitGetRepulsion(this);
}

//-------------------------------------------------------------------------
//
//  SCF convergence stuff.
//

template <class T> double FittedCDImp<T>::operator()(const RVec3& r) const
{
    // No UT coverage
    return FittedFunctionImp<T>::operator()(r);
}

template <class T> void FittedCDImp<T>::Eval(const Mesh& m, Vector<double>& v) const
{
    // No UT coverage
    FittedFunctionImp<T>::Eval(m,v);
}

template <class T> typename FittedCDImp<T>::Vec3 FittedCDImp<T>::Gradient(const RVec3& r) const
{
    // No UT coverage
    return FittedFunctionImp<T>::Gradient(r);
}

template <class T> FittedCD* FittedCDImp<T>::Clone() const
{
    return new FittedCDImp<T>(*this);
}


template class FittedCDImp<double>;
