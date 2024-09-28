// File: BasisSetBrowser.h  Basis set browser, lets clients loop through the list.
#ifndef _BasisSetBrowser_H_
#define _BasisSetBrowser_H_


#include "BasisSet/BasisSet.H"
#include "Misc/ptr_vector1.h"


class BasisFunction;

class BasisSetBrowser
{
public:
    BasisSetBrowser(const BasisSet& bs);

    const BasisFunction* operator++(int)
    {
        current++;
        return *current;
    }
    const BasisFunction& operator* () const
    {
        return **current;
    }
    const BasisFunction* operator->() const
    {
        return *current;
    }
    const BasisFunction* operator& () const
    {
        return *current;
    }
//  operator const BasisFunction*  () const;
    operator bool() const
    {
        return current!=end;
    }
private:
    typedef optr_vector1<BasisFunction*>::const_iterator ITER;
    ITER begin;
    ITER current;
    ITER end;

};

#endif //_BasisSetBrowser_H_
