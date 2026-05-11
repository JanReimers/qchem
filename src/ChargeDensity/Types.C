// File: ChargeDensity/Types.C  Define types used throughout the ChargeDensity module.
module;

export module qchem.ChargeDensity.Types;

export import qchem.BasisSet1.Orbital_HF_IBS;
export import qchem.BasisSet1.Orbital_1E_IBS;
export import qchem.BasisSet1.Orbital_DFT_IBS;
export import qchem.BasisSet1.Fit_IBS;

export namespace qchem::ChargeDensity
{
    using fbs_t=BasisSet1::Fit_IBS;
    using ohfbs_t=BasisSet1::Orbital_HF_IBS<double>;
    template <class T> using tobs_t=BasisSet1::Orbital_1E_IBS<T>;
    template <class T> using todftbs_t=BasisSet1::Orbital_DFT_IBS<T>;
    using obs_t=tobs_t<double>;
    using odftbs_t=todftbs_t<double>;
}
