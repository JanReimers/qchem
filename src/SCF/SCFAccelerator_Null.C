// FIle: SCFAccelerator_Null.C  A simple pass through accerlator proxy that does no acceleration.

#include "Imp/SCF/SCFAccelerator_Null.H"
#include <LASolver.H>

void SCFIrrepAccelerator__Null::UseFD(const SMat& F, const SMat& DPrime)
{
    itsFPrime=itsLaSolver->Transform(F);
}

SCFIrrepAccelerator::SMat SCFIrrepAccelerator__Null::Project()
{
    return itsFPrime;
}


