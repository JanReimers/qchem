// File: WaveFunctionFactories.C  Interface for a wave function.

#include "Imp/WaveFunction/Irrep_WF.H"
#include "Imp/WaveFunction/Polarized_WF.H"
#include "Imp/WaveFunction/UnPolarized_WF.H"

#include <string>
#include <iostream>
#include <typeinfo>
#include <stdlib.h>

//##################################################################
//
//  WaveFunction factory, reads in name of a WaveFunction derived
//  class and make a new object using the default constructor.
//

