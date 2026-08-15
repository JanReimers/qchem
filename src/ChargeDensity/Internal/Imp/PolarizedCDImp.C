// File: PolarizedCD.C  Interface for the charge density.
module;
#include <memory>
#include <cassert>
module qchem.ChargeDensity.Imp.PolarizedCD;
import qchem.ChargeDensity;
import qchem.Symmetry.Spin;

namespace qchem::ChargeDensity
{

//---------------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> tPolarized_CDImp<T>::tPolarized_CDImp()
    : itsSpinUpCD  (nullptr)
    , itsSpinDownCD(nullptr)
{}; // No UT coverage

template <class T> tPolarized_CDImp<T>::tPolarized_CDImp(std::unique_ptr<tDM_CD<T>> up,
                                                         std::unique_ptr<tDM_CD<T>> down)
    : itsSpinUpCD  (std::move(up  ))
    , itsSpinDownCD(std::move(down))
{
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
};
// No destructor: unique_ptr members release the channels, and their move-only nature deletes the copy
// operations that made the old raw-pointer form a double-delete waiting to happen (V1.25).

//-------------------------------------------------------------------
//
//  Access to individual components.
//
template <class T> tDM_CD<T>* tPolarized_CDImp<T>::GetChargeDensity(const Spin& s)
{
    assert(s!=Spin::None);
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    tDM_CD<T>* ret=nullptr;                       // an OBSERVER: the channels stay owned by this object
    if(s==Spin::Up  ) ret=itsSpinUpCD  .get();
    if(s==Spin::Down) ret=itsSpinDownCD.get();
    return ret;
}

template <class T> const tDM_CD<T>* tPolarized_CDImp<T>::GetChargeDensity(const Spin& s) const
{
    assert(s!=Spin::None);
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    const tDM_CD<T>* ret=nullptr;                 // an OBSERVER: the channels stay owned by this object
    if(s==Spin::Up  ) ret=itsSpinUpCD  .get();
    if(s==Spin::Down) ret=itsSpinDownCD.get();
    return ret;
}

template class tPolarized_CDImp<double>;
template class tPolarized_CDImp<dcmplx>;

} //namespace
