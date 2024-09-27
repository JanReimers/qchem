// File: TBasisSetBrowser.h  Basis set browser, lets clients loop through the list.
#ifndef _TBasisSetBrowser_H_
#define _TBasisSetBrowser_H_

#include "BasisSet/BasisSetBrowser.H"
//-----------------------------------------------------
//
//  This can be contructed from a BasisSet (non-template)
//  but does a dynamic cast of all BasisFunctions to
//  TBasisFunction<T>'s for the operators.
//
template <class T> class TBasisSetBrowser
    : public BasisSetBrowser
{
public:
    TBasisSetBrowser(const BasisSet&);

    const TBasisFunction<T>& operator* () const;
    const TBasisFunction<T>* operator->() const;
//  operator const TBasisFunction<T>*  () const;
};

#endif //_TBasisSetBrowser_H_
