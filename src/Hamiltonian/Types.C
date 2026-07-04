// File: Hamiltonian/Types.C  Define types used throughout the Hamiltonian module.
module;

export module qchem.Hamiltonian.Types;

export import qchem.BasisSet;
export import qchem.BasisSet.Orbital_HF_IBS;
export import qchem.BasisSet.Orbital_DFT_IBS;
export import qchem.BasisSet.Orbital_DHF_IBS;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet.Fit_IBS;
import qchem.Types;   // dcmplx (for cobs_t)


export namespace qchem::Hamiltonian
{
    template <class T> using tbs_t=BasisSet::tBasisSet<T>;         // whole (composite) basis: Iterate<tobs_t>() yields the per-irrep bases
    using bs_t    =tbs_t<double>;
    using fbs_t   =BasisSet::Fit_IBS;
    template <class T> using tobs_t=BasisSet::Orbital_1E_IBS<T>;  // T-parametric orbital basis
    // r* = <double>, c* = <dcmplx> (mirrors rsmat_t/chmat_t); bare obs_t transitional (= robs_t).
    using robs_t  =tobs_t<double>;  using cobs_t=tobs_t<dcmplx>;
    using obs_t   =robs_t;
    using ohfbs_t =BasisSet::Orbital_HF_IBS<double>;
    using odftbs_t=BasisSet::Orbital_DFT_IBS<double>;
    using orkbbs_t=BasisSet::Orbital_RKB_IBS<double>;
}
