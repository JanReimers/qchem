// File: OrbitalFactories.C

#include "Imp/Orbitals/TOrbital.H"
#include "Imp/Orbitals/TOrbitals.H"

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
    if (Name==typeid(TOrbitalImp<double>).name()) return new TOrbitalImp<double>;

    std::cerr << "Unknown orbital type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

//##################################################################
//
//  Orbital group factory, reads in name or Orbital group derived
//  class and makes a new object using the default constructor.
//

Orbitals* Orbitals::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(TOrbitalsImp<double>).name()) return new TOrbitalsImp<double>;

    std::cout << "Unknown orbital group type :" << Name << std::endl;
    exit(-1);
    return NULL;
}
