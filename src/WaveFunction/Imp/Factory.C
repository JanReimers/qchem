// File: WaveFunction/Factory.C  Create wavefunctions.
module;
#include <type_traits>
import qchem.SCFAccelerator;

module qchem.WaveFunction.Factory;
import qchem.WaveFunction.Internal.CompositeWF;
import qchem.LASolver;   // qchem::Ortho

namespace qchem::WaveFunction
{

    template <class T> tSCFWaveFunction<T>* Factory(
        const qchem::Hamiltonian::tHamiltonian<T>* h,
        const tbs_t<T>* bs,
        const ElectronConfiguration* ec,
        SCFAccelerators::SCFAccelerator* acc,
        qchem::Ortho basisOrtho, double basisOrthoTol)
    {
        // ONE wave-function class (V1.37); the Hamiltonian's polarization names the IMPOSED SPIN SUBGROUP
        // the composite is built under, for both lineages (SymmetryUpgradePlan §4 tier 4b: a polarized
        // Ham_PW_DFT gets the two-channel Bloch composite just like the molecular path).
        const SpinGroup g = h->IsPolarized() ? SpinGroup::Polarized : SpinGroup::UnPolarized;
        return new tCompositeWF<T>(bs,ec,g,acc,basisOrtho,basisOrthoTol);
    }

    template tSCFWaveFunction<double>* Factory(const qchem::Hamiltonian::tHamiltonian<double>*,
        const tbs_t<double>*, const ElectronConfiguration*, SCFAccelerators::SCFAccelerator*,
        qchem::Ortho, double);
    template tSCFWaveFunction<dcmplx>* Factory(const qchem::Hamiltonian::tHamiltonian<dcmplx>*,
        const tbs_t<dcmplx>*, const ElectronConfiguration*, SCFAccelerators::SCFAccelerator*,
        qchem::Ortho, double);
}


