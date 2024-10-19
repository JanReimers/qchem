// File: ChargeDensityFactories.C  Interface for the charge density category.

#include "ChargeDensity.H"
#include "Imp/ChargeDensity/FittedCD.H"
#include "Imp/ChargeDensity/PolarizedCD.H"
#include "Imp/ChargeDensity/CompositeCD.H"
#include "Imp/ChargeDensity/IrrepCD.H"

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
    if (Name==typeid(IrrepCD<double>          ).name()) return new IrrepCD<double>;
    if (Name==typeid(FittedCDImp<double>).name()) return new FittedCDImp<double>;
    if (Name==typeid(FittedPolarizedCD             ).name()) return new FittedPolarizedCD;
    if (Name==typeid(PolarizedCDImp     ).name()) return new PolarizedCDImp;

    std::cout << "Unknown charge density type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

FittedCD* FittedCD::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(FittedCDImp<double>).name()) return new FittedCDImp<double>;
    if (Name==typeid(FittedPolarizedCD             ).name()) return new FittedPolarizedCD;

    std::cout << "Unknown Fitted charge density type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

