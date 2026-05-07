// File: Orbitals/Types.C  Define types used throughout the Orbitals module.
module;

export module qchem.Orbitals.Types;
// export import qchem.Orbital_HF_IBS;
// export import qchem.Orbital_DFT_IBS;
// export import qchem.Orbital_DHF_IBS;
export import qchem.Orbital_1E_IBS;
// export import qchem.Fit_IBS;

// import qchem.BasisSet1.Orbital_DFT_IBS;
// import qchem.BasisSet1.Fit_IBS;

export namespace qchem::Orbitals
{
    // using fbs_t=Fit_IBS;
    template <class T> using tobs_t=Orbital_IBS<T>;
    // template <class T> using todftbs_t=Orbital_DFT_IBS<T>;
    using obs_t=tobs_t<double>;
    // using ohfbs_t=Orbital_HF_IBS<double>;
    // using odftbs_t=Orbital_DFT_IBS<double>;
    // using orkbbs_t=Orbital_RKB_IBS<double>;
}