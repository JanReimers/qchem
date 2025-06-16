// File: AuxillaryData.C  Static data for the auxilary function.

#include "Molecule/PolarizedGaussian/AuxillaryFJ.H"

const int AuxillaryFJ::thejMax = 16;

double AuxillaryFJ::theLookUp[121][23] =
{
#include "Molecule/PolarizedGaussian/FJLookUp.H"
};
