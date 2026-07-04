// File: ChargeDensity/Types.C  Define types used throughout the ChargeDensity module.
module;

export module qchem.ChargeDensity.Types;

export import qchem.BasisSet.Orbital_HF_IBS;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet.Orbital_DFT_IBS;
export import qchem.BasisSet.Fit_IBS;
import qchem.Types;   // dcmplx (for the c* instantiations)

export namespace qchem::ChargeDensity
{
    // using fbs_t=BasisSet::Fit_IBS;
    using rohfbs_t=BasisSet::Orbital_HF_IBS<double>;
    template <class T> using tobs_t=BasisSet::Orbital_1E_IBS<T>;
    template <class T> using todftbs_t=BasisSet::Orbital_DFT_IBS<T>;
    // r* = <double>, c* = <dcmplx> (mirrors rsmat_t/chmat_t).
    using robs_t=tobs_t<double>;        using cobs_t=tobs_t<dcmplx>;
    using rodftbs_t=todftbs_t<double>;  using codftbs_t=todftbs_t<dcmplx>;
    using odftbs_t=rodftbs_t;
}
