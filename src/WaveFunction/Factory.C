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
    SCFWaveFunction* Factory(
            const Hamiltonian::Hamiltonian*,
            const bs_t* bs,
            const ElectronConfiguration* ec,
            SCFAccelerators::SCFAccelerator* acc
        );
}

