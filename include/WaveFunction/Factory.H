// File: WaveFunction/Factory.H  Create wavefunctions.
#ifndef _WaveFunction_Factory_H_
#define _WaveFunction_Factory_H_

#include <Symmetry/fwd.H>
#include <BasisSet/fwd.H>
#include <Hamiltonian/fwd.H>
#include <SCFAccelerator/fwd.H>
#include <WaveFunction/fwd.H>

namespace WaveFunctionF
{
    WaveFunction* Factory(const Hamiltonian*, const BasisSet* bs,const ElectronConfiguration* ec,SCFAccelerator* acc);
}

#endif //_WaveFunction_Factory_H_