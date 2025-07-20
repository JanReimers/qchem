// File: WaveFunction/Factory.C  Create wavefunctions.

#include <vector>
#include <WaveFunction/Factory.H>
#include "UnPolarized_WF.H"
#include "Polarized_WF.H"
import qchem.Hamiltonian;

namespace WaveFunctionF
{
    WaveFunction* Factory(const Hamiltonian* h, const BasisSet* bs,const ElectronConfiguration* ec,SCFAccelerator* acc)
    {
        return h->IsPolarized() ? (WaveFunction*)new Polarized_WF(bs,ec,acc) : (WaveFunction*)new UnPolarized_WF(bs,ec,acc);
    }
}


