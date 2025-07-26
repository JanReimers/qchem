// File: WaveFunction/Factory.C  Create wavefunctions.
module;
import qchem.SCFAccelerator;

module qchem.WaveFunction.Factory;
import qchem.WaveFunction.Internal.UnPolarizedWF;
import qchem.WaveFunction.Internal.PolarizedWF;

namespace WaveFunctionF
{
    WaveFunction* Factory(const Hamiltonian* h, const BasisSet* bs,const ElectronConfiguration* ec,SCFAccelerator* acc)
    {
        return h->IsPolarized() ? (WaveFunction*)new PolarizedWF(bs,ec,acc) : (WaveFunction*)new UnPolarizedWF(bs,ec,acc);
    }
}


