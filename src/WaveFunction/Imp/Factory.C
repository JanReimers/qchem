// File: WaveFunction/Factory.C  Create wavefunctions.
module;
import qchem.SCFAccelerator;

module qchem.WaveFunction.Factory;
import qchem.WaveFunction.Internal.UnPolarizedWF;
import qchem.WaveFunction.Internal.PolarizedWF;

namespace qchem::WaveFunction
{

    WaveFunction* Factory(
        const Hamiltonian::Hamiltonian* h,
        const bs_t* bs,
        const ElectronConfiguration* ec,
        SCFAccelerators::SCFAccelerator* acc)
    {
        return h->IsPolarized() ? (WaveFunction*)new PolarizedWF(bs,ec,acc) : (WaveFunction*)new UnPolarizedWF(bs,ec,acc);
    }
}


