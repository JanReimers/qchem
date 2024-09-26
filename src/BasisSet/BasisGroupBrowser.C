// File: BasisGroup/BasisGroupBrowser.C  Basis group browser, lets clients loop through the list.



#include "BasisSet/BasisGroupBrowser.H"
#include "BasisSet/BasisGroup.H"
#include <cassert>
BasisGroupBrowser::BasisGroupBrowser(const BasisGroup& bg)
{
    begin=bg.begin();
    current=begin;
    end=bg.end();
}

BasisGroupIterator::BasisGroupIterator(BasisGroup& bg)
{
    begin=bg.begin();
    current=begin;
    end=bg.end();
}
