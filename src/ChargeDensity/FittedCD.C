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
{};

template <class T> FittedCDImp<T>::FittedCDImp(bs_t& bs, mesh_t& m)
    : IntegralConstrainedFF<double>(bs,m) //Use repulsion overlap for fitting
    , itsExactRep(0)
    , itsTotalCharge(0)
{};

template <class T> FittedCDImp<T>::FittedCDImp(bs_t& bs, mesh_t& m, double totalCharge)
    : IntegralConstrainedFF<double>(bs,m) //Use repulsion overlap for fitting
    , itsExactRep(0)
    , itsTotalCharge(totalCharge)
{
    FittedFunctionImp<double>::ReScale(itsTotalCharge);
    assert(itsTotalCharge>0);
    assert(abs(itsTotalCharge-FitGetCharge())<1e-10);
};

template <class T> double FittedCDImp<T>::DoFit(const ScalarFFClient& ffc)
{
    const FittedCD* cd=dynamic_cast<const FittedCD*>(&ffc);
    if (!cd)
    {
        std::cerr << "FittedCDImplementation<T>::DoFit could not cast to charge density" << std::endl;
        exit(-1);
    }
    itsExactRep=cd;
    if  (itsTotalCharge==0) itsTotalCharge=itsExactRep->GetTotalCharge();
    return ConstrainedFF<double>::DoFit(ffc);
}

template <class T> double FittedCDImp<T>::DoFit(const DensityFFClient& ffc)
{
    const ChargeDensity* cd=dynamic_cast<const ChargeDensity*>(&ffc);
    if (!cd)
    {
        std::cerr << "FittedCDImplementation<T>::DoFit could not cast to charge density" << std::endl;
        exit(-1);
    }
    itsExactRep=cd;
    if  (itsTotalCharge==0) itsTotalCharge=itsExactRep->GetTotalCharge();
    return ConstrainedFF<double>::DoFit(ffc);
}
//-----------------------------------------------------------------------------
//
//  Totale energy terms for a charge density.
//

template <class T> ChargeDensity::SMat FittedCDImp<T>::GetRepulsion(const TOrbital_IBS<double>* bs) const
{
    auto bs_dft=dynamic_cast<const TOrbital_DFT_IBS<double>*>(bs);
    assert(bs_dft);
    const std::vector<SMat>& repulsions=bs_dft->Repulsion3C(*itsBasisSet);
    int n=bs->GetNumFunctions();
    SMat J(n,n);
    Fill(J,0.0);
    size_t i=0;
    for (auto c:itsFitCoeff) J+=SMat(c*repulsions[i++]);
    assert(!isnan(J));
    return J;
}

template <class T> ChargeDensity::SMat FittedCDImp<T>::GetExchange(const TOrbital_IBS<double>* bs) const
{
    std::cerr << "FittedCDImplementation<T>::AddExchange: Warning using four center ERIs from a fitted charge density !?!" << std::endl;
    return SMat();
}

template <class T> double FittedCDImp<T>::GetEnergy(const HamiltonianTerm* v) const
{
    assert(itsExactRep);
    return itsExactRep->GetEnergy(v);
}

template <class T> double FittedCDImp<T>::GetSelfRepulsion() const
{
    return 0.5 * FittedFunctionImp<T>::FitGetRepulsion(this);
}

template <class T> double FittedCDImp<T>::GetRepulsion(const FittedFunction* ff) const
{
    const FittedFunctionImp<T>* ffi = dynamic_cast<const FittedFunctionImp<T>*>(ff);
    assert(ffi);
    return FittedFunctionImp<T>::FitGetRepulsion(ffi); //Cross repulsion between to different fitted charge densities!!
}

template <class T> double FittedCDImp<T>::GetOverlap(const FittedFunction* ff) const
{
    const FittedFunctionImp<T>* ffi = dynamic_cast<const FittedFunctionImp<T>*>(ff);
    assert(ffi);
    return FittedFunctionImp<T>::FitGetOverlap(ffi); //Cross repulsion between to different fitted charge densities!!
}

template <class T> double FittedCDImp<T>::GetTotalCharge() const
{
    if (itsTotalCharge==0) itsTotalCharge=FitGetCharge();
    return itsTotalCharge;
}

//------------------------------------------------------------------------------
//
//  Required by fitting routines.
//
template <class T> Vector<double> FittedCDImp<T>::
GetRepulsion3C(const Fit_IBS* fbs) const
{
    return FitGet2CenterRepulsion(fbs);
}

//-------------------------------------------------------------------------
//
//  SCF convergence stuff.
//
template <class T> void FittedCDImp<T>::MixIn(const ChargeDensity& cd,double c)
{
    const FittedCDImp<T>* fcd = dynamic_cast<const FittedCDImp<T>*>(&cd);
    assert(fcd);
    FitMixIn(*fcd,c);
}

template <class T> double FittedCDImp<T>::GetChangeFrom(const ChargeDensity& cd) const
{
    const FittedCDImp<T>* fcd = dynamic_cast<const FittedCDImp<T>*>(&cd);
    assert(fcd);
    return FitGetChangeFrom(*fcd);
}

//-------------------------------------------------------------------------
//
//  Real space function stuff.
//
template <class T> void   FittedCDImp<T>::ReScale (double factor)
{
    itsTotalCharge*=factor;
    FittedFunctionImp<T>::ReScale(factor);
}

template <class T> void FittedCDImp<T>::ShiftOrigin(const RVec3& newCenter)
{
    FittedFunctionImp<T>::ShiftOrigin(newCenter);
}

template <class T> double FittedCDImp<T>::operator()(const RVec3& r) const
{
    return FittedFunctionImp<T>::operator()(r);
}

template <class T> void FittedCDImp<T>::Eval(const Mesh& m, Vector<double>& v) const
{
    FittedFunctionImp<T>::Eval(m,v);
}

template <class T> typename FittedCDImp<T>::Vec3 FittedCDImp<T>::Gradient(const RVec3& r) const
{
    return FittedFunctionImp<T>::Gradient(r);
}

//-----------------------------------------------------------------------
//
//  Streamable stuff.
//
template <class T> std::ostream& FittedCDImp<T>::Write(std::ostream& os) const
{
    ConstrainedFF<T>::Write(os);
    if (StreamableObject::Binary())
        BinaryWrite(itsTotalCharge,os);
    else
        os << itsTotalCharge << " ";
    return os;
}

template <class T> std::istream& FittedCDImp<T>::Read(std::istream&is)
{
    ConstrainedFF<T>::Read(is);
    if (StreamableObject::Binary())
        BinaryRead(itsTotalCharge,is);
    else
        is >> itsTotalCharge;
    return is;
}

template <class T> FittedCD* FittedCDImp<T>::Clone() const
{
    return new FittedCDImp<T>(*this);
}


template class FittedCDImp<double>;
