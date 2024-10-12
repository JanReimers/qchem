#ifndef ITERATIONPARAMS_H_INCLUDED
#define ITERATIONPARAMS_H_INCLUDED

#include "types.H"



struct SCFIterationParams
{
    size_t NMaxIter;             //Max allowed number of iterations
    double MinDeltaRo;           //Minimum delta in charge density for convergence.
    double StartingRelaxRo;      //relaxation for mixing Ro.  Dynamically adjusted during iterations.
    double kT;                   //Thermal smearing of electron occupations near fermi level.
};


#endif // ITERATIONPARAMS_H_INCLUDED
