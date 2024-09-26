// File: BasisSet/TBasisSetBrowser.C  Basis set browser, lets clients loop through the list.



#include "BasisSet/TBasisFunction.H"
#include "BasisSet/TBasisSet.H"
#include "BasisSet/TBasisSetBrowser.H"
#include "BasisSetImplementation/BasisSetImplementation.H"
#include <cassert>
#include <complex>

//#######################################################################
//
//  Browser constructor..
//

template <class T> TBasisSetBrowser<T>::TBasisSetBrowser(const BasisSet& bs)
    : BasisSetBrowser(bs)
{};

template <class T> const TBasisFunction<T>& TBasisSetBrowser<T>::operator*() const
{
    const TBasisFunction<T>* tbf=dynamic_cast<const TBasisFunction<T>*>(BasisSetBrowser::operator&());
    if (!tbf)
    {
        const BasisFunction* bf=BasisSetBrowser::operator&();
        tbf=dynamic_cast<const TBasisFunction<T>*>(bf);
    }
    assert(tbf);
    return *tbf;
}

template <class T> const TBasisFunction<T>* TBasisSetBrowser<T>::operator->() const
{
    const TBasisFunction<T>* tbf=dynamic_cast<const TBasisFunction<T>*>(BasisSetBrowser::operator&());
    if (!tbf)
    {
        const BasisFunction* bf=BasisSetBrowser::operator&();
        tbf=dynamic_cast<const TBasisFunction<T>*>(bf);
    }
    assert(tbf);
    return tbf;
}

template class TBasisSetBrowser<double>;
template class TBasisSetBrowser<std::complex<double> >;
