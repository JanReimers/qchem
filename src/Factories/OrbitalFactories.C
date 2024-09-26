// File: OrbitalFactories.C

#include "OrbitalImplementation/TOrbitalImplementation.H"
#include "OrbitalImplementation/TOrbitalGroupImplementation.H"

#include <string>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>

//##################################################################
//
//  Orbital factory, reads in name or Orbital derived
//  class and makes a new object using the default constructor.
//

Orbital* Orbital::Factory(std::istream& is)
{
    std::string Name=PeekAtName(is);
    if (Name=="OrbitalImplementation") return new TOrbitalImplementation<double>;

    std::cerr << "Unknown orbital type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

//##################################################################
//
//  Orbital group factory, reads in name or Orbital group derived
//  class and makes a new object using the default constructor.
//

OrbitalGroup* OrbitalGroup::Factory(std::istream& is)
{
    std::string Name=PeekAtName(is);
    if (Name==typeid(OrbitalGroupImplementation).name()) return new TOrbitalGroupImplementation<double>;

    std::cout << "Unknown orbital group type :" << Name << std::endl;
    exit(-1);
    return NULL;
}
