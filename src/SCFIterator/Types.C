// File: SCFIterator/Types.C  Define types used throughout the SCFIterator module.
module;

export module qchem.SCFIterator.Types;
// export import qchem.Orbital_1E_IBS;
export import qchem.BasisSet;

// import qchem.BasisSet1.Orbital_DFT_IBS;
// import qchem.BasisSet1.Fit_IBS;

export namespace qchem::SCFIterator
{
    // using obs_t=Orbital_IBS<double>;
    using bs_t=BasisSet;
}