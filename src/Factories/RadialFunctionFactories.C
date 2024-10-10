// File: RadialFunctionFactories.C

#include "Imp/BasisSet/PolarizedGaussian/Radial/GaussianRF.H"
#include "Imp/BasisSet/PolarizedGaussian/Radial/ContractedGaussianRF.H"

#include <string>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>

//##################################################################
//
//  Radial function factory, reads in name or RadialFunction derived
//  class and make a new object using the default constructor.
//
RadialFunction* RadialFunction::Factory(std::istream& is)
{
    std::string Name=PeekAtName(is);
    if (Name==          typeid(GaussianRF).name()) return new           GaussianRF;
    if (Name==typeid(ContractedGaussianRF).name()) return new ContractedGaussianRF;

    std::cout << "Unknown radial function type :" << Name << std::endl;
    exit(-1);
    return NULL;
}


