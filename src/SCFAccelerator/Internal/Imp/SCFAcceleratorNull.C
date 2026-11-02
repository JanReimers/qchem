// FIle: SCFAcceleratorNull.C  A simple pass through accerlator proxy that does no acceleration.
module qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull;
import qchem.Reporting;   // the accelerator files what it was BUILT with, one level into its builder's context

namespace qchem::SCFAccelerators
{
SCFAcceleratorNull::SCFAcceleratorNull() {report::Set("type", "Null");}   // nothing else to state: no knobs

template <class T> void tSCFIrrepAcceleratorNull<T>::UseFD(const hmat_t<T>& F, const hmat_t<T>& DPrime)
{
    itsFPrime=itsLASolver->Transform(F);
    itsEn=GradientNorm<T>(OrbitalGradient<T>(itsFPrime, DPrime));
}

template <class T> typename LASolver<T>::UUd_t tSCFIrrepAcceleratorNull<T>::NextOrbitals()
{
    return itsLASolver->SolveOrtho(itsFPrime);
}

template class tSCFIrrepAcceleratorNull<double>;
template class tSCFIrrepAcceleratorNull<dcmplx>;

} //namespace

