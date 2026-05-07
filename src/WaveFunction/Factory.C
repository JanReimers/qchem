// File: WaveFunction/Factory.H  Create wavefunctions.
module;
export module qchem.WaveFunction.Factory;
import qchem.WaveFunction.Types;
export import qchem.SCFAccelerator;
export import qchem.WaveFunction;
export import qchem.Hamiltonian;
export import qchem.Symmetry.ElectronConfiguration;

export namespace qchem::WaveFunction
{
    WaveFunction* Factory(
            const Hamiltonian::Hamiltonian*,
            const bs_t* bs,
            const ElectronConfiguration* ec,
            SCFAccelerators::SCFAccelerator* acc
        );
}

