// File: Orbitals/Factory.C  Create orbitals.
export module qchem.Orbitals.Factory;
import qchem.Orbitals.Internal.OrbitalsImp;

export namespace qchem::Orbitals
{
    TOrbitals<double>* Factory(const Orbital_IBS<double>* bs, Spin ms);
}

