//  File: MeshFactories.C



#include "Mesh/AtomMesh.H"
#include "Mesh/MoleculeMesh.H"
#include "Mesh/UnitCellMesh.H"
#include <string>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>

//##################################################################
//
//  mesh factory, reads in name or BasisFunction derived
//  class and make a new object using the default constructor.
//
Mesh* Mesh::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==    typeid(AtomMesh).name()) return new     AtomMesh;
    if (Name==typeid(MoleculeMesh).name()) return new MoleculeMesh;
    if (Name==typeid(UnitCellMesh).name()) return new UnitCellMesh;

    std::cout << "Unknown mesh type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

