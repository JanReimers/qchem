// File: Orbitals/Factory.H  Create orbitals.
#ifndef _Orbitals_Factory_H_
#define _Orbitals_Factory_H_

#include <Symmetry/fwd.H>
#include <BasisSet/fwd.H>
#include <Orbitals/fwd.H>
import qchem.Symmetry.Spin;

namespace OrbitalsF
{
    TOrbitals<double>* Factory(const TOrbital_IBS<double>* bs, Spin ms);
}

#endif //_Orbitals_Factory_H_