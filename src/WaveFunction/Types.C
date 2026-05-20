// File: WaveFunction/Types.C  Define types used throughout the WaveFunction module.
module;

export module qchem.WaveFunction.Types;

export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet;

export namespace qchem::WaveFunction
{
    using obs_t=BasisSet::Orbital_1E_IBS<double>;
    using bs_t=BasisSet::BasisSet<double>;
}

