// File: WaveFunction/Types.C  Define types used throughout the WaveFunction module.
module;

export module qchem.WaveFunction.Types;

export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet;

export namespace qchem::WaveFunction
{
    template <class T> using tobs_t=BasisSet::Orbital_1E_IBS<T>;
    using robs_t=tobs_t<double>;
    template <class T> using tbs_t=BasisSet::tBasisSet<T>;
    using rbs_t=tbs_t<double>;
}

