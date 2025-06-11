// File: SCFParams.H Parameters used for controling SCF iteration and convergence.
#ifndef _SCFParams_H_
#define _SCFParams_H_

#include "oml/imp/index_t.h"

struct SCFParams
{
    size_t NMaxIter;         //Max allowed number of iterations
    double MinDeltaRo;       //Minimum delta in charge density for convergence.
    double MinDelE;          //Minimum delta energy convergence.
    double MinError;         //Minimum error from SCF accelerator.  i.e. [F,D] (orbital gradient).
    double StartingRelaxRo;  //relaxation for mixing Ro.  Dynamically adjusted during iterations.
    bool   Verbose;          //Display iteration details.        
};


#endif // _SCFParams_H_
