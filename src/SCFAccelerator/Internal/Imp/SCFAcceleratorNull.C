// FIle: SCFAcceleratorNull.C  A simple pass through accerlator proxy that does no acceleration.
module qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull;
void SCFIrrepAcceleratorNull::UseFD(const SMat& F, const SMat& DPrime)
{
    itsFPrime=itsLaSolver->Transform(F);
}

SCFIrrepAccelerator::SMat SCFIrrepAcceleratorNull::Project()
{
    return itsFPrime;
}


