// File: WaveFunction/Factory.C  Create wavefunctions.
module;
#include <type_traits>
import qchem.SCFAccelerator;

module qchem.WaveFunction.Factory;
import qchem.WaveFunction.Internal.UnPolarizedWF;
import qchem.WaveFunction.Internal.PolarizedWF;

namespace qchem::WaveFunction
{

    template <class T> tSCFWaveFunction<T>* Factory(
        const qchem::Hamiltonian::tHamiltonian<T>* h,
        const tbs_t<T>* bs,
        const ElectronConfiguration* ec,
        SCFAccelerators::tSCFAccelerator<T>* acc)
    {
        // The plane-wave (dcmplx) lineage is unpolarized; only the real path builds a PolarizedWF
        // (spin density is a real Gaussian-basis facility), so tPolarizedWF<dcmplx> is never needed.
        if constexpr (std::is_same_v<T,dcmplx>)
            return (tSCFWaveFunction<T>*)new tUnPolarizedWF<T>(bs,ec,acc);
        else
            return h->IsPolarized() ? (tSCFWaveFunction<T>*)new tPolarizedWF<T>(bs,ec,acc)
                                    : (tSCFWaveFunction<T>*)new tUnPolarizedWF<T>(bs,ec,acc);
    }

    template tSCFWaveFunction<double>* Factory(const qchem::Hamiltonian::tHamiltonian<double>*,
        const tbs_t<double>*, const ElectronConfiguration*, SCFAccelerators::tSCFAccelerator<double>*);
    template tSCFWaveFunction<dcmplx>* Factory(const qchem::Hamiltonian::tHamiltonian<dcmplx>*,
        const tbs_t<dcmplx>*, const ElectronConfiguration*, SCFAccelerators::tSCFAccelerator<dcmplx>*);
}


