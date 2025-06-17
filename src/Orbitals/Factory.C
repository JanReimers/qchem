// File: Orbitals/Factory.H  Create orbitals.

#include <Orbitals/Factory.H>
#include "TOrbitals.H"

namespace OrbitalsF 
{
    TOrbitals<double>* Factory(const TOrbital_IBS<double>* bs, Spin ms)
    {
        return new  TOrbitalsImp<double>(bs,ms);
    }
}


