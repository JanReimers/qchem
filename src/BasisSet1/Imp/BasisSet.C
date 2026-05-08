// File: BasisSetImp.C Quantum Chemistry basis set expressed as a sequence of Irrep basis sets.
module;
module qchem.BasisSet1;
import qchem.BasisSet1.Orbital_DFT_IBS;

namespace BasisSet1
{
template <class T> Fit_IBS* BasisSet<T>::CreateCDFitBasisSet(const Cluster* cl) const
{   
    auto dft=*Iterate<Orbital_DFT_IBS<double>>().begin();
    return dft->CreateCDFitBasisSet(cl);
}
template <class T> Fit_IBS* BasisSet<T>::CreateVxcFitBasisSet(const Cluster* cl) const
{
    auto dft=*Iterate<Orbital_DFT_IBS<double>>().begin();
    return dft->CreateVxcFitBasisSet(cl);
}

template class BasisSet<double>;

} //namespace