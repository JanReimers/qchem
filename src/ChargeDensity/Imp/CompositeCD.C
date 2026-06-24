// File: Imp/CompositeCD.C  Composite charged density, which is any array of Irrep DM_CDs.
module;
#include <cassert>
#include <vector>
#include <memory>
module qchem.CompositeCD;
import qchem.ChargeDensity.Types;
import qchem.Blaze;

namespace qchem::ChargeDensity
{

//------------------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> tComposite_CD<T>::tComposite_CD()
{};

template <class T> void tComposite_CD<T>::Insert(tDM_CD<T>* cd)
{
    itsCDs.push_back(std::unique_ptr<tDM_CD<T>>(cd));
}

//-----------------------------------------------------------------------------
//
//  Total energy terms for a charge density.
//
template <class T> void tComposite_CD<T>::AccumulateDirect(hmat_t<T>& Jab, const ohfbs_t* bs_ab) const
{
    for (auto& c:itsCDs) c->AccumulateDirect(Jab,bs_ab);
}

template <class T> void tComposite_CD<T>::AccumulateExchange(hmat_t<T>& Kab, const ohfbs_t* bs_ab) const
{
    for (auto& c:itsCDs) c->AccumulateExchange(Kab,bs_ab);
}

template <class T> double tComposite_CD<T>::DM_Contract(const tStatic_CC<T>* v) const
{
    double ret=0.0;
    for (auto& c:itsCDs) ret+=c->DM_Contract(v);
    return ret;
}

template <class T> double tComposite_CD<T>::DM_Contract(const tDynamic_CC<T>* v,const tDM_CD<T>* cd) const
{
    double ret=0.0;
    for (auto& c:itsCDs) ret+=c->DM_Contract(v,cd);
    return ret;
}

template <class T> double tComposite_CD<T>::GetTotalCharge() const
{
    double ret=0.0;
    for (auto& c:itsCDs) ret+=c->GetTotalCharge();
    return ret;
}

//------------------------------------------------------------------------------
//
//  Required by fitting routines.
//
template <class T> rvec_t tComposite_CD<T>::GetRepulsion3C(const fbs_t* fbs) const
{
    rvec_t ret(fbs->GetNumFunctions(),0);
    for (auto& c:itsCDs) ret+=c->GetRepulsion3C(fbs);
    return ret;
}

//-------------------------------------------------------------------------
//
//  SCF convergence stuff.
//
template <class T> void tComposite_CD<T>::ReScale(double factor)
{
    // No UT coverage
    for (auto& c:itsCDs) c->ReScale(factor);
}

template <class T> void tComposite_CD<T>::MixIn(const tDM_CD<T>& cd,double f)
{
    const tComposite_CD* ecd = dynamic_cast<const tComposite_CD*>(&cd);
    assert(ecd);
    auto  b(ecd->itsCDs.begin());
    for (auto& c:itsCDs)
    {
        c->MixIn(**b,f);
        b++;
    }
}

template <class T> double tComposite_CD<T>::GetChangeFrom(const tDM_CD<T>& cd) const
{
    const tComposite_CD* ecd = dynamic_cast<const tComposite_CD*>(&cd);
    assert(ecd);
    assert(itsCDs.size()==ecd->itsCDs.size());
    auto  b(ecd->itsCDs.begin());
    double ret=0;
    for (auto& c:itsCDs)
    {
        ret += c->GetChangeFrom(**b);
        b++;
    }
    return ret;
}

//-------------------------------------------------------------------------
//
//  Real space function stuff.
//
template <class T> double tComposite_CD<T>::operator()(const rvec3_t& r) const
{
    double ret=0.0;
    for (auto& c:itsCDs) ret+=c->operator()(r);
    return ret;
}

template <class T> rvec3_t tComposite_CD<T>::Gradient  (const rvec3_t& r) const
{
    // No UT coverage
    rvec3_t ret(0,0,0);
    for (auto& c:itsCDs) ret+=c->Gradient(r);
    return ret;
}

template class tComposite_CD<double>;
template class tComposite_CD<dcmplx>;

} //namespace