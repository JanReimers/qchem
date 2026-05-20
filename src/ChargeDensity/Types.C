// File: ChargeDensity/Types.C  Define types used throughout the ChargeDensity module.
module;

export module qchem.ChargeDensity.Types;

export import qchem.BasisSet.Orbital_HF_IBS;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet.Orbital_DFT_IBS;
export import qchem.BasisSet.Fit_IBS;

export namespace qchem::ChargeDensity
{
    using fbs_t=BasisSet::Fit_IBS;
    using ohfbs_t=BasisSet::Orbital_HF_IBS<double>;
    template <class T> using tobs_t=BasisSet::Orbital_1E_IBS<T>;
    template <class T> using todftbs_t=BasisSet::Orbital_DFT_IBS<T>;
    using obs_t=tobs_t<double>;
    using odftbs_t=todftbs_t<double>;
}
