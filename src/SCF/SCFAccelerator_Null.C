// FIle: SCFAccelerator_Null.C  A simple pass through accerlator proxy that does no acceleration.

#include "Imp/SCF/SCFAccelerator_Null.H"
#include <LASolver.H>

void SCFIrrepAccelerator__Null::Init(const LASolver<double>* las) 
{
     itsLaSolver=las;
}; 
void SCFIrrepAccelerator__Null::UseFD(const SMat& F, const SMat& DPrime)
{
    itsFPrime=itsLaSolver->Transform(F);
}

SCFIrrepAccelerator::Mat SCFIrrepAccelerator__Null::CalculateError()
{
    assert(false);
    return Mat();
}
SCFIrrepAccelerator::SMat SCFIrrepAccelerator__Null::Project()
{
    return itsFPrime;
}


