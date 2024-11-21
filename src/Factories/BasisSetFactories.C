// File: BasisSetFactories.C

#include "BasisSet.H"
#include "Imp/BasisSet/SphericalGaussian/BasisFunction.H"
#include "Imp/BasisSet/SphericalGaussian/IrrepBasisSet.H"
#include "Imp/BasisSet/SphericalGaussian/BasisSet.H"
#include "Imp/BasisSet/PolarizedGaussian/BasisFunction.H"
#include "Imp/BasisSet/PolarizedGaussian/IrrepBasisSet.H"
#include "Imp/BasisSet/PolarizedGaussian/BasisSet.H"
#include "Imp/BasisSet/PolarizedGaussian/Block.H"

#include "Imp/Symmetry/YlmQN.H"
#include "Imp/Symmetry/YlQN.H"
#include "Imp/Symmetry/UnitQN.H"
//#include "BasisSetImplementation/PlaneWave/BlochQN.H"

#include <string>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>


//##################################################################
//
//  Basis function factory, reads in name or BasisFunction derived
//  class and make a new object using the default constructor.
//
BasisFunction* BasisFunction::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(SphericalGaussian::BasisFunction).name()) return new SphericalGaussian::BasisFunction;
    if (Name==typeid(PolarizedGaussian::BasisFunction).name()) return new PolarizedGaussian::BasisFunction;
//    if (Name==typeid(PlaneWaveBF        ).name()) return new         PlaneWaveBF;

    std::cout << "Unknown basis function type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

namespace PolarizedGaussian
{

Block* Block::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(Block).name()) return new Block;

    std::cout << "Unknown basis function block type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

} //namespace PolarizedGaussian

//##################################################################
//
//  Basis set factory, reads in name or BasisSet derived
//  class and make a new object using the default constructor.
//

IrrepBasisSet* IrrepBasisSet::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(SphericalGaussian::IrrepBasisSet).name()) return new SphericalGaussian::IrrepBasisSet;
    if (Name==typeid(PolarizedGaussian::IrrepBasisSet).name()) return new PolarizedGaussian::IrrepBasisSet;
//    if (Name==typeid(        PlaneWaveBS).name()) return new         PlaneWaveBS;

    std::cout << "Unknown irrep basis set type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

BasisSet* BasisSet::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(SphericalGaussian::BasisSet).name()) return new SphericalGaussian::BasisSet;
    if (Name==typeid(PolarizedGaussian::BasisSet).name()) return new PolarizedGaussian::BasisSet;

    std::cout << "Unknown basis group type :" << Name << std::endl;
    exit(-1);
    return NULL;
}


//##################################################################
//
//  Quantum Number factory, reads in name or BasisFunction derived
//  class and make a new object using the default constructor.
//
QuantumNumber* QuantumNumber::Factory(std::istream& is)
{
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(YlQN).name()) return new YlQN;
    if (Name==typeid(     UnitQN).name()) return new      UnitQN;
//    if (Name==typeid(            BlochQN).name()) return new             BlochQN;

    std::cout << "Unknown Quantum Number type :" << Name << std::endl;
    exit(-1);
    return NULL;
}

