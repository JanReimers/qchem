// File: Hamiltonian/Types.C  Define types used throughout the Hamiltonian module.
module;

export module qchem.Hamiltonian.Types;

export import qchem.BasisSet;
export import qchem.BasisSet.Orbital_HF_IBS;
export import qchem.BasisSet.Orbital_DFT_IBS;
export import qchem.BasisSet.Orbital_DHF_IBS;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet.Fit_IBS;


export namespace qchem::Hamiltonian
{
    using bs_t    =BasisSet::BasisSet<double>;
    using fbs_t   =BasisSet::Fit_IBS;
    using obs_t   =BasisSet::Orbital_1E_IBS<double>;
    using ohfbs_t =BasisSet::Orbital_HF_IBS<double>;
    using odftbs_t=BasisSet::Orbital_DFT_IBS<double>;
    using orkbbs_t=BasisSet::Orbital_RKB_IBS<double>;
}
