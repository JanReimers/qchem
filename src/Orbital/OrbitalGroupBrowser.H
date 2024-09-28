// File: OrbitalGroupBrowser.H  Orbital group browser.
#ifndef _OrbitalGroupBrowser_H_
#define _OrbitalGroupBrowser_H_

#include "Misc/ptr_vector1.h"
#include "Orbital/OrbitalGroup.H"

class OrbitalGroupBrowser
{
public:
    OrbitalGroupBrowser(const OrbitalGroup& og);

    const Orbital* operator++(int)
    {
        current++;
        return *current;
    }
    const Orbital& operator* () const
    {
        return **current;
    }
    const Orbital* operator->() const
    {
        return *current;
    }
    const Orbital* operator& () const
    {
        return *current;
    }
//  operator const BasisFunction*  () const;
    operator bool() const
    {
        return current!=end;
    }
private:
    typedef optr_vector1<Orbital*>::const_iterator ITER;
    ITER begin;
    ITER current;
    ITER end;

};

class OrbitalGroupIterator
{
public:
    OrbitalGroupIterator(OrbitalGroup& og);

    Orbital* operator++(int)
    {
        current++;
        return *current;
    }
    Orbital& operator* () const
    {
        return **current;
    }
    Orbital* operator->() const
    {
        return *current;
    }
    Orbital* operator& () const
    {
        return *current;
    }
//  operator const BasisFunction*  () const;
    operator bool() const
    {
        return current!=end;
    }
private:
    typedef optr_vector1<Orbital*>::iterator ITER;
    ITER begin;
    ITER current;
    ITER end;

};


#endif //_OrbitalGroupBrowser_H_
