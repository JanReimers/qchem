// File: AuxillaryData.C  Static data for the auxilary function.

#include "BasisSetImplementation/PolarizedGaussian/Auxillary/AuxillaryFJ.H"

const int AuxillaryFJ::thejMax = 16;

double AuxillaryFJ::theLookUp[121][23] =
{
#include "BasisSetImplementation/PolarizedGaussian/Auxillary/FJLookUp.H"
};
