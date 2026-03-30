// File: FittedCDImplementation.C  General implementation using a density matrix.
module;
#include <cmath>
#include <cassert>
#include <memory>
#include <vector>
#include "blaze/Math.h"

module qchem.ChargeDensity.Imp.FittedCD;
import qchem.Orbital_DFT_IBS;
import qchem.Mesh;
import qchem.Blaze;

//------------------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> FittedCDImp<T>::FittedCDImp(bs_t& bs, mesh_t& m, double totalCharge)
    : IntegralConstrainedFF<double>(bs,m) //Use repulsion overlap for fitting
{
    FittedFunctionImp<double>::ReScale(totalCharge);
    assert(totalCharge>0);
    assert(fabs(totalCharge-FitGetCharge())<1e-10);
};


//-----------------------------------------------------------------------------
//
//  Totale energy terms for a charge density.
//

template <class T> smat_t<T> FittedCDImp<T>::GetRepulsion(const Orbital_DFT_IBS<double>* bs) const
{
    assert(bs);
    const ERI3<T>& repulsions=bs->Repulsion3C(*itsBasisSet);
    int n=bs->GetNumFunctions();
    smat_t<T> J=zero<T>(n);
    size_t i=0;
    for (auto c:itsFitCoeff) J+=c*repulsions[i++];
    assert(!isnan(J));
    return J;
}

template <class T> double FittedCDImp<T>::GetSelfRepulsion() const
{
    return 0.5 * FittedFunctionImp<T>::FitGetRepulsion(this);
}

//-------------------------------------------------------------------------
//
//  SCF convergence stuff.
//

template <class T> double FittedCDImp<T>::operator()(const rvec3_t& r) const
{
    // No UT coverage
    return FittedFunctionImp<T>::operator()(r);
}

template <class T> rvec3_t FittedCDImp<T>::Gradient(const rvec3_t& r) const
{
    // No UT coverage
    return FittedFunctionImp<T>::Gradient(r);
}

template <class T> FittedCD* FittedCDImp<T>::Clone() const
{
    return new FittedCDImp<T>(*this);
}


template class FittedCDImp<double>;
