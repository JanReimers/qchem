// File: WaveFunctionFactories.C  Interface for a wave function.

#include "Imp/WaveFunction/IrrepWaveFunction.H"
#include "Imp/WaveFunction/MasterPolarizedWF.H"
#include "Imp/WaveFunction/MasterUnPolarizedWF.H"

#include <string>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>

//##################################################################
//
//  WaveFunction factory, reads in name of a WaveFunction derived
//  class and make a new object using the default constructor.
//

IrrepWaveFunction* IrrepWaveFunction::Factory(std::istream& is)
{
    std::string Name=PeekAtName(is);
    // if (Name==typeid(IrrepWaveFunction).name()) return new IrrepWaveFunction;

    std::cout << "Unknown IrrepWaveFunction type :" << Name << std::endl;
    exit(-1);
    return NULL;
}
WaveFunction* WaveFunction::Factory(std::istream& is)
{
    std::string Name=PeekAtName(is);
    
    std::cout << "Unknown WaveFunction type :" << Name << std::endl;
    exit(-1);
    return NULL;
}
