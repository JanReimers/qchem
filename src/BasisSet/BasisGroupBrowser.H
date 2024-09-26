// File: BasisGroupBrowser.h  Basis set browser, lets clients loop through the list.
#ifndef _BasisGroupBrowser_H_
#define _BasisGroupBrowser_H_


#include "BasisSet/BasisSet.H"
#include "Misc/ptr_vector.h"

class BasisGroup;
class BasisSet;

//TODO: All these browsers and iterators should be template classes.
class BasisGroupBrowser
{
public:
    BasisGroupBrowser(const BasisGroup& bg);

    const BasisSet* operator++(int)
    {
        current++;
        return &current;
    }
    const BasisSet& operator* () const
    {
        return *current;
    }
    const BasisSet* operator->() const
    {
        return &current;
    }
    const BasisSet* operator& () const
    {
        return &current;
    }
//  operator const BasisSet*  () const;
    operator bool() const
    {
        return current!=end;
    }
private:
    typedef optr_vector<BasisSet*>::const_iterator ITER;
    ITER begin;
    ITER current;
    ITER end;
};

class BasisGroupIterator
{
public:
    BasisGroupIterator(BasisGroup& bg);

    BasisSet* operator++(int)
    {
        current++;
        return &current;
    }
    BasisSet& operator* () const
    {
        return *current;
    }
    BasisSet* operator->() const
    {
        return &current;
    }
    BasisSet* operator& () const
    {
        return &current;
    }
//  operator const BasisFunction*  () const;
    operator bool() const
    {
        return current!=end;
    }
private:
    typedef optr_vector<BasisSet*>::iterator ITER;
    ITER begin;
    ITER current;
    ITER end;

};


#endif //_BasisGroupBrowser_H_
