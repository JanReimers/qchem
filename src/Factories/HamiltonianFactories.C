// File: HamiltonainFactories.C  Interface a Hamiltonian operator.

#include "Imp/Hamiltonian/ExchangeFunctional.H"
#include "Imp/Hamiltonian/Kinetic.H"
#include "Imp/Hamiltonian/Vee.H"
#include "Imp/Hamiltonian/Ven.H"
#include "Imp/Hamiltonian/Vnn.H"
#include "Imp/Hamiltonian/Vxc.H"
#include "Imp/Hamiltonian/VxcPol.H"
#include "Imp/Hamiltonian/FittedVee.H"
#include "Imp/Hamiltonian/FittedVxc.H"
#include "Imp/Hamiltonian/FittedVxcPol.H"
#include "Imp/Hamiltonian/LDAVxc.H"
#include "Imp/Hamiltonian/Hamiltonian.H"
#include <BasisSet.H>

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

Static_HT* Static_HT::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(Kinetic  ).name()) return new Kinetic ();
    if (Name==typeid(Vee ).name()) return new Vee();
    if (Name==typeid(Ven ).name()) return new Ven();
    if (Name==typeid(Vnn ).name()) return new Vnn();
    if (Name==typeid(LDAVxc   ).name()) return new LDAVxc();
    if (Name==typeid(Vxc).name()) return new Vxc();
    if (Name==typeid(FittedVee).name()) return new  FittedVee();
    if (Name==typeid(FittedVxc).name()) return new FittedVxc();
    if (Name==typeid(FittedVxcPol).name()) return new FittedVxcPol();
    if (Name==typeid(VxcPol).name()) return new VxcPol();

    std::cout << "Unknown HamiltonianTerm type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

FittablePotential* FittablePotential::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);

    std::cout << "Unknown potential type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

