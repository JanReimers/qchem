// File: Orbitals/Factory.H  Create orbitals.
#ifndef _Orbitals_Factory_H_
#define _Orbitals_Factory_H_

#include <vector>
#include <Orbitals/fwd.H>
import qchem.Irrep_BS;
import qchem.Symmetry.Spin;
import qchem.Symmetry.ElectronConfiguration;

namespace OrbitalsF
{
    TOrbitals<double>* Factory(const TOrbital_IBS<double>* bs, Spin ms);
}

#endif //_Orbitals_Factory_H_