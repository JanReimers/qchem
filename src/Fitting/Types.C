// File: Fitting/Types.C  Define types used throughout the Fitting module.
module;

export module qchem.Fitting.Types;

export import qchem.BasisSet1.Orbital_DFT_IBS;
export import qchem.BasisSet1.Fit_IBS;

export namespace qchem::Fitting
{
    using fbs_t=BasisSet1::Fit_IBS;
    template <class T> using obs_t=BasisSet1::Orbital_DFT_IBS<T>;
}
