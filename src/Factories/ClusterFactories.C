// File: ClusterFactories.C



#include "Imp/Cluster/Lattice.H"
#include "Imp/Cluster/Molecule.H"
#include "Imp/Cluster/Atom.H"
#include "Imp/Cluster/UnitCell.H"
#include <string>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>

//##################################################################
//
//  Basis set factory, reads in name or BasisSet derived
//  class and make a new object using the default constructor.
//

Cluster* Cluster::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(Molecule).name()) return new Molecule;
    if (Name==typeid(Lattice ).name()) return new Lattice ;

    std::cout << "Unknown cluster type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

Atom* Atom::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(Atom).name()) return new Atom;

    std::cout << "Unknown atom type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

UnitCell* UnitCell::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(UnitCell).name()) return new UnitCell;

    std::cout << "Unknown UnitCell type :" << Name << std::endl;
    exit(-1);
    return NULL;
}
