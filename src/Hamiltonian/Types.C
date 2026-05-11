// File: Hamiltonian/Types.C  Define types used throughout the Hamiltonian module.
module;

export module qchem.Hamiltonian.Types;

export import qchem.BasisSet1;
export import qchem.BasisSet1.Orbital_HF_IBS;
export import qchem.BasisSet1.Orbital_DFT_IBS;
export import qchem.BasisSet1.Orbital_DHF_IBS;
export import qchem.BasisSet1.Orbital_1E_IBS;
export import qchem.BasisSet1.Fit_IBS;


export namespace qchem::Hamiltonian
{
    using bs_t    =BasisSet1::BasisSet<double>;
    using fbs_t   =BasisSet1::Fit_IBS;
    using obs_t   =BasisSet1::Orbital_1E_IBS<double>;
    using ohfbs_t =BasisSet1::Orbital_HF_IBS<double>;
    using odftbs_t=BasisSet1::Orbital_DFT_IBS<double>;
    using orkbbs_t=BasisSet1::Orbital_RKB_IBS<double>;
}
