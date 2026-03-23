// FIle: SCFAcceleratorNull.C  A simple pass through accerlator proxy that does no acceleration.
module qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull;
void SCFIrrepAcceleratorNull::UseFD(const smat_t<double>& F, const smat_t<double>& DPrime)
{
    itsFPrime=itsLaSolver_blaze->Transform(F);
}

smat_t<double> SCFIrrepAcceleratorNull::Project()
{
    return itsFPrime;
}


