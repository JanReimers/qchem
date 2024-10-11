// File: FittedCDImplementation.C  General implementation using a density matrix.



#include "ChargeDensityImplementation/FittedCDImplementation.H"
#include "oml/smatrix.h"
#include "oml/imp/binio.h"
#include <cmath>
#include <cassert>
#include <stdlib.h>

//------------------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> FittedCDImplementation<T>::FittedCDImplementation()
    : IntegralConstrainedFF<double>()
    , itsExactRep(0)
    , itsTotalCharge(0)
{};

template <class T> FittedCDImplementation<T>::FittedCDImplementation(const rc_ptr<IrrepBasisSet>& bs, Mesh* m)
    : IntegralConstrainedFF<double>(bs,m,true) //Use repulsion overlap for fitting
    , itsExactRep(0)
    , itsTotalCharge(0)
{};

template <class T> FittedCDImplementation<T>::FittedCDImplementation(const rc_ptr<IrrepBasisSet>& bs, Mesh* m, double totalCharge)
    : IntegralConstrainedFF<double>(bs,m,true) //Use repulsion overlap for fitting
    , itsExactRep(0)
    , itsTotalCharge(totalCharge)
{
    FittedFunctionImplementation<double>::ReScale(itsTotalCharge);
    assert(itsTotalCharge>0);
    assert(abs(itsTotalCharge-FitGetCharge())<1e-10);
};

template <class T> double FittedCDImplementation<T>::DoFit(const FittedFunctionClient& ffc)
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
template <class T> ChargeDensity::SMat FittedCDImplementation<T>::GetOverlap  (const IrrepBasisSet* bs) const
{
    const FittedFunctionImplementation<T>* ffi=dynamic_cast<const FittedFunctionImplementation<T>*>(this);
    assert(ffi);
    const std::vector<SMat>& overlap=bs->GetOverlap3C(ffi->itsBasisSet.get());
    int n=bs->GetNumFunctions();
    SMat J(n,n);
    Fill(J,0.0);
    size_t i=0;
    for (auto c:ffi->itsFitCoeff) J+=SMat(c*overlap[i++]);
    assert(!isnan(J));
    return J;
}

template <class T> ChargeDensity::SMat FittedCDImplementation<T>::GetRepulsion(const IrrepBasisSet* bs) const
{
    const FittedFunctionImplementation<T>* ffi=dynamic_cast<const FittedFunctionImplementation<T>*>(this);
    assert(ffi);
    const std::vector<SMat>& repulsions=bs->GetRepulsion3C(ffi->itsBasisSet.get());
    int n=bs->GetNumFunctions();
    SMat J(n,n);
    Fill(J,0.0);
    size_t i=0;
    for (auto c:ffi->itsFitCoeff) J+=SMat(c*repulsions[i++]);
    assert(!isnan(J));
    return J;
}

template <class T> ChargeDensity::SMat FittedCDImplementation<T>::GetExchange(const IrrepBasisSet* bs) const
{
    std::cerr << "FittedCDImplementation<T>::AddExchange: Warning using four center ERIs from a fitted charge density !?!" << std::endl;
    return SMat();
}

template <class T> double FittedCDImplementation<T>::GetEnergy(const HamiltonianTerm* v) const
{
    assert(itsExactRep);
    return itsExactRep->GetEnergy(v);
}

template <class T> double FittedCDImplementation<T>::GetSelfRepulsion() const
{
    return 0.5 * FittedFunctionImplementation<T>::FitGetRepulsion(this);
}

template <class T> double FittedCDImplementation<T>::GetRepulsion(const FittedFunction* ff) const
{
    const FittedFunctionImplementation<T>* ffi = dynamic_cast<const FittedFunctionImplementation<T>*>(ff);
    assert(ffi);
    return FittedFunctionImplementation<T>::FitGetRepulsion(ffi); //Cross repulsion between to different fitted charge densities!!
}

template <class T> double FittedCDImplementation<T>::GetOverlap(const FittedFunction* ff) const
{
    const FittedFunctionImplementation<T>* ffi = dynamic_cast<const FittedFunctionImplementation<T>*>(ff);
    assert(ffi);
    return FittedFunctionImplementation<T>::FitGetOverlap(ffi); //Cross repulsion between to different fitted charge densities!!
}

template <class T> double FittedCDImplementation<T>::GetTotalCharge() const
{
    if (itsTotalCharge==0) itsTotalCharge=FitGetCharge();
    return itsTotalCharge;
}

//------------------------------------------------------------------------------
//
//  Required by fitting routines.
//
template <class T> void FittedCDImplementation<T>::
InjectOverlaps  (FittedFunction* ff, const IrrepBasisSet* theFitBasisSet) const
{
    FittedFunctionImplementation<T>* ffi=dynamic_cast<FittedFunctionImplementation<T>*>(ff);
    assert(ffi);
    ffi->GetFitCoeff()+=FitGet2CenterOverlap(theFitBasisSet);
}

template <class T> void FittedCDImplementation<T>::
InjectRepulsions(FittedFunction* ff, const IrrepBasisSet* theFitBasisSet) const
{
    FittedFunctionImplementation<T>* ffi=dynamic_cast<FittedFunctionImplementation<T>*>(ff);
    assert(ffi);
    ffi->GetFitCoeff()+=FitGet2CenterRepulsion(theFitBasisSet);
}

//-------------------------------------------------------------------------
//
//  SCF convergence stuff.
//
template <class T> void FittedCDImplementation<T>::MixIn(const ChargeDensity& cd,double c)
{
    const FittedCDImplementation<T>* fcd = dynamic_cast<const FittedCDImplementation<T>*>(&cd);
    assert(fcd);
    FitMixIn(*fcd,c);
}

template <class T> double FittedCDImplementation<T>::GetChangeFrom(const ChargeDensity& cd) const
{
    const FittedCDImplementation<T>* fcd = dynamic_cast<const FittedCDImplementation<T>*>(&cd);
    assert(fcd);
    return FitGetChangeFrom(*fcd);
}

//-------------------------------------------------------------------------
//
//  Real space function stuff.
//
template <class T> void   FittedCDImplementation<T>::ReScale (double factor)
{
    itsTotalCharge*=factor;
    FittedFunctionImplementation<T>::ReScale(factor);
}

template <class T> void FittedCDImplementation<T>::ShiftOrigin(const RVec3& newCenter)
{
    FittedFunctionImplementation<T>::ShiftOrigin(newCenter);
}

template <class T> double FittedCDImplementation<T>::operator()(const RVec3& r) const
{
    return FittedFunctionImplementation<T>::operator()(r);
}

template <class T> void FittedCDImplementation<T>::Eval(const Mesh& m, Vector<double>& v) const
{
    FittedFunctionImplementation<T>::Eval(m,v);
}

template <class T> typename FittedCDImplementation<T>::Vec3 FittedCDImplementation<T>::Gradient(const RVec3& r) const
{
    return FittedFunctionImplementation<T>::Gradient(r);
}

//-----------------------------------------------------------------------
//
//  Streamable stuff.
//
template <class T> std::ostream& FittedCDImplementation<T>::Write(std::ostream& os) const
{
    ConstrainedFF<T>::Write(os);
    if (StreamableObject::Binary())
        BinaryWrite(itsTotalCharge,os);
    else
        os << itsTotalCharge << " ";
    return os;
}

template <class T> std::istream& FittedCDImplementation<T>::Read(std::istream&is)
{
    ConstrainedFF<T>::Read(is);
    if (StreamableObject::Binary())
        BinaryRead(itsTotalCharge,is);
    else
        is >> itsTotalCharge;
    return is;
}

template <class T> FittedCD* FittedCDImplementation<T>::Clone() const
{
    return new FittedCDImplementation<T>(*this);
}


template class FittedCDImplementation<double>;
