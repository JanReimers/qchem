// File: Imp/Factory.C  Create orbitals.
module qchem.Orbitals.Factory;
import qchem.Orbitals.Internal.OrbitalsImp;

namespace OrbitalsF 
{
    TOrbitals<double>* Factory(const TOrbital_IBS<double>* bs, Spin ms)
    {
        return new  TOrbitalsImp<double>(bs,ms);
    }
}


