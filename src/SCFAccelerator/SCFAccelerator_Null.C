// FIle: SCFAccelerator_Null.C  A simple pass through accerlator proxy that does no acceleration.

#include <LASolver/LASolver.H>
#include "SCFAccelerator_Null.H"

void SCFIrrepAccelerator__Null::UseFD(const SMat& F, const SMat& DPrime)
{
    itsFPrime=itsLaSolver->Transform(F);
}

SCFIrrepAccelerator::SMat SCFIrrepAccelerator__Null::Project()
{
    return itsFPrime;
}


