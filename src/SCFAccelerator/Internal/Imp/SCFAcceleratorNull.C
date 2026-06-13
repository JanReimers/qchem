// FIle: SCFAcceleratorNull.C  A simple pass through accerlator proxy that does no acceleration.
module qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull;

namespace qchem::SCFAccelerators
{
void SCFIrrepAcceleratorNull::UseFD(const smat_t<double>& F, const smat_t<double>& DPrime)
{
    itsFPrime=itsLASolver->Transform(F);
}

LASolver<double>::UUd_t SCFIrrepAcceleratorNull::NextOrbitals()
{
    return itsLASolver->SolveOrtho(itsFPrime);
}

} //namespace

