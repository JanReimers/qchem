// File: BasisSetImp.C Quantum Chemistry basis set expressed as a sequence of Irrep basis sets.
module;
#include <cassert>
module qchem.BasisSet;
import qchem.BasisSet.Orbital_DFT_IBS;

namespace BasisSet
{
template <class T> Fit_IBS* BasisSet<T>::CreateCDFitBasisSet(const Structure* cl) const
{
    auto dft=*Iterate<Orbital_DFT_IBS<double>>().begin();
    return dft->CreateCDFitBasisSet(cl);
}
template <class T> Fit_IBS* BasisSet<T>::CreateVxcFitBasisSet(const Structure* cl) const
{
    auto dft=*Iterate<Orbital_DFT_IBS<double>>().begin();
    return dft->CreateVxcFitBasisSet(cl);
}

// Auxiliary DFT-fit basis sets are a real Gaussian-basis facility; the plane-wave (dcmplx) lineage
// does not fit densities (it expands them in the plane-wave basis directly), so these are NA there.
template <> Fit_IBS* BasisSet<dcmplx>::CreateCDFitBasisSet (const Structure*) const { assert(false); return 0; }
template <> Fit_IBS* BasisSet<dcmplx>::CreateVxcFitBasisSet(const Structure*) const { assert(false); return 0; }

template class BasisSet<double>;
template class BasisSet<dcmplx>;

} //namespace