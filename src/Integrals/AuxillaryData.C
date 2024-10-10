// File: AuxillaryData.C  Static data for the auxilary function.

#include "Imp/Integrals/AuxillaryFJ.H"

const int AuxillaryFJ::thejMax = 16;

double AuxillaryFJ::theLookUp[121][23] =
{
#include "Imp/Integrals/FJLookUp.H"
};
