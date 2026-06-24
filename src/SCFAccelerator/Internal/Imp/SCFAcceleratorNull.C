// FIle: SCFAcceleratorNull.C  A simple pass through accerlator proxy that does no acceleration.
module qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull;

namespace qchem::SCFAccelerators
{
template <class T> void tSCFIrrepAcceleratorNull<T>::UseFD(const hmat_t<T>& F, const hmat_t<T>& DPrime)
{
    itsFPrime=itsLASolver->Transform(F);
}

template <class T> typename LASolver<T>::UUd_t tSCFIrrepAcceleratorNull<T>::NextOrbitals()
{
    return itsLASolver->SolveOrtho(itsFPrime);
}

template class tSCFIrrepAcceleratorNull<double>;
template class tSCFIrrepAcceleratorNull<dcmplx>;
template class tSCFAcceleratorNull<double>;
template class tSCFAcceleratorNull<dcmplx>;

} //namespace

