// File: Orbitals/Factory.C  Create orbitals.
export module qchem.Orbitals.Factory;
import qchem.Orbitals.Types;
import qchem.Orbitals.Internal.OrbitalsImp;

export namespace qchem::Orbitals
{
    TOrbitals<double>* Factory(const robs_t* bs, Spin ms);
    TOrbitals<dcmplx>* Factory(const tobs_t<dcmplx>* bs, Spin ms);   // plane-wave (complex) orbitals
}

