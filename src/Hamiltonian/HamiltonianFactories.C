// File: HamiltonainFactories.C  Interface a Hamiltonian operator.

#include "ExchangeFunctional.H"
#include "Kinetic.H"
#include "Vee.H"
#include "Ven.H"
#include "Vnn.H"
#include "Vxc.H"
#include "VxcPol.H"
#include "FittedVee.H"
#include "FittedVxc.H"
#include "FittedVxcPol.H"
#include "LDAVxc.H"
#include "Hamiltonian.H"
#include <BasisSet/BasisSet.H>

#include <string>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>

//##################################################################
//
//  Hamiltonian factory, reads in name of a Hamiltonian derived
//  class and make a new object using the default constructor.
//

Hamiltonian* Hamiltonian::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(HamiltonianImp).name()) return new HamiltonianImp;

    std::cout << "Unknown Hamiltonian type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

//##################################################################
//
//  HamiltonianTerm factory, reads in name of a HamiltonianTerm derived
//  class and make a new object using the default constructor.
//


