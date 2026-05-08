// File: WaveFunction/Types.C  Define types used throughout the WaveFunction module.
module;

export module qchem.WaveFunction.Types;

#ifdef LegacyBasisSet
export import qchem.Orbital_1E_IBS;
export import qchem.BasisSet;

export namespace qchem::WaveFunction
{
    using obs_t=Orbital_IBS<double>;
    using bs_t=BasisSet;
}
#else
export import qchem.BasisSet1.Orbital_1E_IBS;
export import qchem.BasisSet1;

export namespace qchem::WaveFunction
{
    using obs_t=BasisSet1::Orbital_1E_IBS<double>;
    using bs_t=BasisSet1::BasisSet<double>;
}

#endif