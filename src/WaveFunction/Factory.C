// File: WaveFunction/Factory.H  Create wavefunctions.
module;
export module qchem.WaveFunction.Factory;
export import qchem.SCFAccelerator;
export import qchem.WaveFunction;
export import qchem.Hamiltonian;
export import qchem.BasisSet;
export import qchem.Symmetry.ElectronConfiguration;

export namespace WaveFunctionF
{
    WaveFunction* Factory(const Hamiltonian*, const BasisSet* bs,const ElectronConfiguration* ec,SCFAccelerator* acc);
}

