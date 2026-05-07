// File: Orbitals/Factory.C  Create orbitals.
export module qchem.Orbitals.Factory;
import qchem.Orbitals.Types;
import qchem.Orbitals.Internal.OrbitalsImp;

export namespace qchem::Orbitals
{
    TOrbitals<double>* Factory(const obs_t* bs, Spin ms);
}

