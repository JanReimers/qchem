// File: BasisSet/BasisSetBrowser.C  Basis set browser, lets clients loop through the list.



#include "BasisSet/BasisSetBrowser.H"
#include "BasisSetImplementation/BasisSetImplementation.H"
#include <cassert>
BasisSetBrowser::BasisSetBrowser(const BasisSet& bs)
{
    const BasisSetImplementation* bsi=dynamic_cast<const BasisSetImplementation*>(&bs);
    assert(bsi);
    begin=bsi->begin();
    current=begin;
    end=bsi->end();
}

