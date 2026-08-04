// File: PolarizedCD.C  Interface for the charge density.
module;
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
    : itsSpinUpCD  (0)
    , itsSpinDownCD(0)
{}; // No UT coverage

template <class T> tPolarized_CDImp<T>::tPolarized_CDImp(tDM_CD<T>* up, tDM_CD<T>* down)
    : itsSpinUpCD  (up  )
    , itsSpinDownCD(down)
{
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
};

template <class T> tPolarized_CDImp<T>::~tPolarized_CDImp()
{
    delete itsSpinUpCD;
    delete itsSpinDownCD;
}

//-------------------------------------------------------------------
//
//  Access to individual components.
//
template <class T> tDM_CD<T>* tPolarized_CDImp<T>::GetChargeDensity(const Spin& s)
{
    assert(s!=Spin::None);
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    tDM_CD<T>* ret=0;
    if(s==Spin::Up  ) ret=itsSpinUpCD  ;
    if(s==Spin::Down) ret=itsSpinDownCD;
    return ret;
}

template <class T> const tDM_CD<T>* tPolarized_CDImp<T>::GetChargeDensity(const Spin& s) const
{
    assert(s!=Spin::None);
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    const tDM_CD<T>* ret=0;
    if(s==Spin::Up  ) ret=itsSpinUpCD  ;
    if(s==Spin::Down) ret=itsSpinDownCD;
    return ret;
}

template class tPolarized_CDImp<double>;
template class tPolarized_CDImp<dcmplx>;

} //namespace
