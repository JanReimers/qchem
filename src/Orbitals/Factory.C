// File: Orbitals/Factory.C  Create orbitals.
export module qchem.Orbitals.Factory;
import qchem.Orbitals.Internal.OrbitalsImp;

export namespace OrbitalsF
{
    TOrbitals<double>* Factory(const TOrbital_IBS<double>* bs, Spin ms);
}

