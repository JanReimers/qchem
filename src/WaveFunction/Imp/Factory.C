// File: WaveFunction/Factory.C  Create wavefunctions.
module;
#include <type_traits>
import qchem.SCFAccelerator;

module qchem.WaveFunction.Factory;
import qchem.WaveFunction.Internal.UnPolarizedWF;
import qchem.WaveFunction.Internal.PolarizedWF;
import qchem.LASolver;   // qchem::Ortho

namespace qchem::WaveFunction
{

    template <class T> tSCFWaveFunction<T>* Factory(
        const qchem::Hamiltonian::tHamiltonian<T>* h,
        const tbs_t<T>* bs,
        const ElectronConfiguration* ec,
        SCFAccelerators::tSCFAccelerator<T>* acc,
        qchem::Ortho basisOrtho, double basisOrthoTol)
    {
        // Both lineages dispatch on the Hamiltonian's polarization (SymmetryUpgradePlan §4 tier 4b):
        // a polarized Ham_PW_DFT builds the two-channel Bloch wavefunction just like the molecular path.
        return h->IsPolarized() ? (tSCFWaveFunction<T>*)new tPolarizedWF<T>(bs,ec,acc,basisOrtho,basisOrthoTol)
                                : (tSCFWaveFunction<T>*)new tUnPolarizedWF<T>(bs,ec,acc,basisOrtho,basisOrthoTol);
    }

    template tSCFWaveFunction<double>* Factory(const qchem::Hamiltonian::tHamiltonian<double>*,
        const tbs_t<double>*, const ElectronConfiguration*, SCFAccelerators::tSCFAccelerator<double>*,
        qchem::Ortho, double);
    template tSCFWaveFunction<dcmplx>* Factory(const qchem::Hamiltonian::tHamiltonian<dcmplx>*,
        const tbs_t<dcmplx>*, const ElectronConfiguration*, SCFAccelerators::tSCFAccelerator<dcmplx>*,
        qchem::Ortho, double);
}


