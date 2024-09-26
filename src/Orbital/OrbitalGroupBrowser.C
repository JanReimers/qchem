//  file: OrbitalGroupBrowser.C  Orbital group browser.



#include "Orbital/Orbital.H"
#include "Orbital/OrbitalGroup.H"
#include "Orbital/OrbitalGroupBrowser.H"
#include "OrbitalImplementation/OrbitalGroupImplementation.H"
#include <cassert>

//#######################################################################
//
//  Browser constructor..
//

OrbitalGroupBrowser::OrbitalGroupBrowser(const OrbitalGroup& og)
{
    const OrbitalGroupImplementation* ogi=dynamic_cast<const OrbitalGroupImplementation*>(&og);
    assert(ogi);
    begin=ogi->itsOrbitals.begin();
    current=begin;
    end=ogi->itsOrbitals.end();
};


//#######################################################################
//
//  Iterator constructor..
//
OrbitalGroupIterator::OrbitalGroupIterator(OrbitalGroup& og)
{
    OrbitalGroupImplementation* ogi=dynamic_cast<OrbitalGroupImplementation*>(&og);
    assert(ogi);
    begin=ogi->itsOrbitals.begin();
    current=begin;
    end=ogi->itsOrbitals.end();
};
