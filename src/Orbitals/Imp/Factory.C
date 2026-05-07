// File: Imp/Factory.C  Create orbitals.
module qchem.Orbitals.Factory;
import qchem.Orbitals.Internal.OrbitalsImp;

namespace qchem::Orbitals
{
    TOrbitals<double>* Factory(const obs_t* bs, Spin ms)
    {
        return new  TOrbitalsImp<double>(bs,ms);
    }
}


