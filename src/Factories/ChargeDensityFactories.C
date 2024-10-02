// File: ChargeDensityFactories.C  Interface for the charge density category.

#include "ChargeDensity.H"
#include "ChargeDensityImplementation/FittedCDImplementation.H"
#include "ChargeDensityImplementation/PolarizedCDImplementation.H"
#include "ChargeDensityImplementation/CompositeCD/CompositeCD.H"
#include "ChargeDensityImplementation/ExactIrrepCD/ExactIrrepCD.H"

#include <string>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>

//##################################################################
//
//  ChargeDensity factory, reads in name of a ChargeDensity derived
//  class and make a new object using the default constructor.
//

ChargeDensity* ChargeDensity::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(CompositeCD                   ).name()) return new CompositeCD;
    if (Name==typeid(ExactIrrepCD<double>          ).name()) return new ExactIrrepCD<double>;
    if (Name==typeid(FittedCDImplementation<double>).name()) return new FittedCDImplementation<double>;
    if (Name==typeid(FittedPolarizedCD             ).name()) return new FittedPolarizedCD;
    if (Name==typeid(PolarizedCDImplementation     ).name()) return new PolarizedCDImplementation;

    std::cout << "Unknown charge density type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

FittedCD* FittedCD::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(FittedCDImplementation<double>).name()) return new FittedCDImplementation<double>;
    if (Name==typeid(FittedPolarizedCD             ).name()) return new FittedPolarizedCD;

    std::cout << "Unknown Fitted charge density type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

