// File: WaveFunction/Factory.C  Create wavefunctions.

#include <WaveFunction/Factory.H>
#include <Hamiltonian/Hamiltonian.H>
#include "UnPolarized_WF.H"
#include "Polarized_WF.H"

namespace WaveFunctionF
{
    WaveFunction* Factory(const Hamiltonian* h, const BasisSet* bs,const ElectronConfiguration* ec,SCFAccelerator* acc)
    {
        return h->IsPolarized() ? (WaveFunction*)new Polarized_WF(bs,ec,acc) : (WaveFunction*)new UnPolarized_WF(bs,ec,acc);
    }
}


