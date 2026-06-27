// File: FittedCDImplementation.C  General implementation using a density matrix.
module;
#include <cassert>
#include <memory>
#include <vector>

module qchem.ChargeDensity.Imp.FittedCD;
import qchem.Blaze;
import qchem.Math;

namespace qchem::ChargeDensity
{

// typedef std::shared_ptr<const BasisSet::FIT_CD_ABS> fbs_t; 
//------------------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> FittedCDImp<T>::FittedCDImp(fbs_t& bs, double totalCharge)
    // Charge-CONSTRAINED Coulomb-metric density fit (Dunlap-Connolly-Sabin 1979): every DoFit yields a
    // density of exactly totalCharge, variationally -- no post-hoc rescale needed.
    : itsFitter(Fitting::MakeDensityFitter(bs))
{
    itsFitter->ReScale(totalCharge);   // normalize the initial guess (each DoFit then re-imposes the charge)
    assert(totalCharge>0);
    assert(fabs(totalCharge-itsFitter->Integral())<1e-10);
};

//-----------------------------------------------------------------------------
//
//  Total energy terms for a charge density -- the fitter answers the "your repulsion with this basis?".
//
template <class T> smat_t<T> FittedCDImp<T>::GetRepulsion(const odftbs_t* bs) const
{
    return itsFitter->Repulsion(bs);   // Sum_a c_a <Oi|f_a/r12|Oj>
}

template <class T> double FittedCDImp<T>::GetSelfRepulsion() const
{
    return 0.5 * itsFitter->FitGetSelfRepulsion();   // 1/2 <ro|1/r12|ro>
}

template <class T> FittedCD* FittedCDImp<T>::Clone() const
{
    // Unused today.  A correct Clone needs a POLYMORPHIC fitter clone (so the constrained fitter isn't
    // sliced); implement when Clone is actually needed -- e.g. building a polarized CD from unpolarized.
    assert(false && "FittedCDImp::Clone not implemented -- see polarized-CD TODO");
    return nullptr;
}


template class FittedCDImp<double>;

} //namespace