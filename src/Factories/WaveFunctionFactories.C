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

