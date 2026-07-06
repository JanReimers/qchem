// File: Fitting/Types.C  Define types used throughout the Fitting module.
module;

export module qchem.Fitting.Types;

export import qchem.BasisSet.Orbital_DFT_IBS;
export import qchem.BasisSet.Fit_IBS;

export namespace qchem::Fitting
{
    // The orbital basis the fitter contracts against is named by its COMMON base (Orbital_1E_IBS), so both
    // the Gaussian (Orbital_DFT_IBS, 3-centre) and plane-wave (Band_FT_IBS, ΔG_Map) bases qualify; each
    // concrete fitter dynamic_casts down to the capability it needs (abstract->abstract, sanctioned).
    template <class T> using robs_t=BasisSet::Orbital_1E_IBS<T>;
}
