// FIle: SCFAcceleratorNull.C  A simple pass through accerlator proxy that does no acceleration.
module qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull;
void SCFIrrepAcceleratorNull::UseFD(const SMatrix<double>& F, const SMatrix<double>& DPrime)
{
    itsFPrime=itsLaSolver->Transform(F);
}

SMatrix<double> SCFIrrepAcceleratorNull::Project()
{
    return itsFPrime;
}


