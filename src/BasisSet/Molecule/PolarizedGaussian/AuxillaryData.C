// File: AuxillaryData.C  Static data for the auxilary function.

#include "Imp/BasisSet/Molecule/PolarizedGaussian/AuxillaryFJ.H"

const int AuxillaryFJ::thejMax = 16;

double AuxillaryFJ::theLookUp[121][23] =
{
#include "Imp/BasisSet/Molecule/PolarizedGaussian/FJLookUp.H"
};
