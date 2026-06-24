// File: WaveFunction/Factory.H  Create wavefunctions.
module;
export module qchem.WaveFunction.Factory;
import qchem.WaveFunction.Types;
export import qchem.SCFAccelerator;
export import qchem.WaveFunction.SCF;
export import qchem.Hamiltonian;
export import qchem.ElectronConfiguration;

export namespace qchem::WaveFunction
{
    // Templated on the matrix element type T (rX/cX); T is deduced from the arguments, so existing
    // <double> callers are unchanged.  Builds an (un)polarized WF for double; the dcmplx plane-wave
    // path is always unpolarized.
    template <class T> tSCFWaveFunction<T>* Factory(
            const qchem::Hamiltonian::tHamiltonian<T>*,
            const tbs_t<T>* bs,
            const ElectronConfiguration* ec,
            SCFAccelerators::tSCFAccelerator<T>* acc
        );
}

