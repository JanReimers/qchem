// File: WaveFunctionFactories.C  Interface for a wave function.

#include "WaveFunctionImp/IrrepWaveFunction/IrrepWaveFunction.H"
#include "WaveFunctionImp/WaveFunctionGroup/WaveFunctionGroup.H"

#include <string>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>

//##################################################################
//
//  WaveFunction factory, reads in name of a WaveFunction derived
//  class and make a new object using the default constructor.
//

WaveFunction* WaveFunction::Factory(std::istream& is)
{
    std::string Name=PeekAtName(is);
    if (Name==typeid(IrrepWaveFunction).name()) return new IrrepWaveFunction;
    if (Name==typeid(WaveFunctionGroup).name()) return new WaveFunctionGroup;

    std::cout << "Unknown WaveFunction type :" << Name << std::endl;
    exit(-1);
    return NULL;
}
