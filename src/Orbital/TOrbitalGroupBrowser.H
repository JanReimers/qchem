// File: TOrbitalGroupBrowser.H  Templated orbital group browser/iterator.
#ifndef _TOrbitalGroupBrowser_H_
#define _TOrbitalGroupBrowser_H_



#include "Orbital/OrbitalGroupBrowser.H"

template <class T> class TOrbital;
class OrbitalGroup;

template <class T> class TOrbitalGroupBrowser
    : public OrbitalGroupBrowser
{
public:
    TOrbitalGroupBrowser(const OrbitalGroup&);
    const TOrbital<T>& operator* () const;
    const TOrbital<T>* operator->() const;
};

template <class T> class TOrbitalGroupIterator
    : public OrbitalGroupIterator
{
public:
    TOrbitalGroupIterator(OrbitalGroup&);
    TOrbital<T>& operator* ();
    TOrbital<T>* operator->();
};

#endif //_TOrbitalGroupBrowser_H_
