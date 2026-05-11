// File: Orbitals/Types.C  Define types used throughout the Orbitals module.
module;

export module qchem.Orbitals.Types;

export import qchem.BasisSet1.Orbital_1E_IBS;

export namespace qchem::Orbitals
{
    template <class T> using tobs_t=BasisSet1::Orbital_1E_IBS<T>;
    using obs_t=tobs_t<double>;
}
