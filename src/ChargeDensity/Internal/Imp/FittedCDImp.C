// File: FittedCDImplementation.C  General implementation using a density matrix.
module;
#include <cassert>
#include <memory>
#include <vector>

module qchem.ChargeDensity.Imp.FittedCD;
import qchem.Mesh;
import qchem.Blaze;
import qchem.Math;

namespace qchem::ChargeDensity
{

//------------------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> FittedCDImp<T>::FittedCDImp(bs_t& bs, mesh_t& m, double totalCharge)
    // Charge-CONSTRAINED Coulomb-metric density fit (Dunlap-Connolly-Sabin 1979): every DoFit yields a
    // density of exactly totalCharge, variationally -- no post-hoc rescale needed.
    : itsFitter(std::make_unique<Fitting::IntegralConstrainedFF<double>>(bs,m))
{
    itsFitter->ReScale(totalCharge);   // normalize the initial guess (each DoFit then re-imposes the charge)
    assert(totalCharge>0);
    assert(fabs(totalCharge-itsFitter->FitGetCharge())<1e-10);
};

template <class T> FittedCDImp<T>::FittedCDImp(const FittedCDImp& o)
    // NOTE: copies the base fitter -- slices the IntegralConstrainedFF constraint.  Harmless today (Clone
    // is unused); revisit (polymorphic fitter clone) when Clone is needed, e.g. polarized CD from unpolarized.
    : itsFitter(std::make_unique<Fitting::FunctionFitter<double>>(*o.itsFitter))
{};


//-----------------------------------------------------------------------------
//
//  Total energy terms for a charge density.
//

template <class T> smat_t<T> FittedCDImp<T>::GetRepulsion(const odftbs_t* bs) const
{
    assert(bs);
    const ERI3<T>& repulsions=bs->Repulsion3C(*itsFitter->itsBasisSet);
    int n=bs->GetNumFunctions();
    smat_t<T> J=blazem::zero<T>(n);
    size_t i=0;
    for (auto c:itsFitter->itsFitCoeff) J+=c*repulsions[i++];
    assert(!blazem::isnan(J));
    return J;
}

template <class T> double FittedCDImp<T>::GetSelfRepulsion() const
{
    return 0.5 * itsFitter->FitGetRepulsion(itsFitter.get());
}

//-------------------------------------------------------------------------
//
//  FittedFunction mixing (delegated to the composed fitter; the "other" is another FittedCDImp).
//
template <class T> void FittedCDImp<T>::FitMixIn(const Fitting::FittedFunction& g,double f)
{
    auto* og=dynamic_cast<const FittedCDImp<T>*>(&g);
    assert(og);
    itsFitter->FitMixIn(*og->itsFitter,f);
}

template <class T> double FittedCDImp<T>::FitGetChangeFrom(const Fitting::FittedFunction& g) const
{
    auto* og=dynamic_cast<const FittedCDImp<T>*>(&g);
    assert(og);
    return itsFitter->FitGetChangeFrom(*og->itsFitter);
}

template <class T> FittedCD* FittedCDImp<T>::Clone() const
{
    return new FittedCDImp<T>(*this);
}


template class FittedCDImp<double>;

} //namespace