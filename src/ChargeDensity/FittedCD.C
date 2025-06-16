// File: FittedCDImplementation.C  General implementation using a density matrix.



#include "Imp/ChargeDensity/FittedCD.H"
#include <ChargeDensity.H>
#include <Mesh.H>
#include <BasisSet/DFT_IBS.H>
#include "oml/smatrix.h"
#include "oml/imp/binio.h"
#include <cmath>
#include <cassert>
#include <stdlib.h>

//------------------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> FittedCDImp<T>::FittedCDImp(bs_t& bs, mesh_t& m, double totalCharge)
    : IntegralConstrainedFF<double>(bs,m) //Use repulsion overlap for fitting
{
    FittedFunctionImp<double>::ReScale(totalCharge);
    assert(totalCharge>0);
    assert(abs(totalCharge-FitGetCharge())<1e-10);
};


//-----------------------------------------------------------------------------
//
//  Totale energy terms for a charge density.
//

template <class T> DM_CD::SMat FittedCDImp<T>::GetRepulsion(const TOrbital_DFT_IBS<double>* bs) const
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
