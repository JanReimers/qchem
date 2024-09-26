// File: TOrbitalGroupBrowser.C  Templated orbital group browser/Iterator.



#include "Orbital/TOrbital.H"
#include "Orbital/TOrbitalGroupBrowser.H"
#include <cassert>

//#######################################################################
//
//  Browser constructor..
//
template <class T> TOrbitalGroupBrowser<T>::TOrbitalGroupBrowser(const OrbitalGroup& og)
    : OrbitalGroupBrowser(og)
{};

template <class T> const TOrbital<T>& TOrbitalGroupBrowser<T>::operator*() const
{
    const TOrbital<T>* to = dynamic_cast<const TOrbital<T>*>(&OrbitalGroupBrowser::operator*());
    assert(to);
    return *to;
}

template <class T> const TOrbital<T>* TOrbitalGroupBrowser<T>::operator->() const
{
    const TOrbital<T>* to = dynamic_cast<const TOrbital<T>*>(&OrbitalGroupBrowser::operator*());
    assert(to);
    return to;
}
//#######################################################################
//
//  Iterator constructor..
//
template <class T> TOrbitalGroupIterator<T>::TOrbitalGroupIterator(OrbitalGroup& og)
    : OrbitalGroupIterator(og)
{};

template <class T> TOrbital<T>& TOrbitalGroupIterator<T>::operator*()
{
    TOrbital<T>* to = dynamic_cast<TOrbital<T>*>(&OrbitalGroupIterator::operator*());
    assert(to);
    return *to;
}

template <class T> TOrbital<T>* TOrbitalGroupIterator<T>::operator->()
{
    TOrbital<T>* to = dynamic_cast<TOrbital<T>*>(&OrbitalGroupIterator::operator*());
    assert(to);
    return to;
}

template class TOrbitalGroupBrowser <double>;
template class TOrbitalGroupIterator<double>;
template class TOrbitalGroupBrowser <std::complex<double> >;
template class TOrbitalGroupIterator<std::complex<double> >;
