// File: WaveFunction/Types.C  Define types used throughout the WaveFunction module.
module;

export module qchem.WaveFunction.Types;
export import qchem.Orbital_1E_IBS;
export import qchem.BasisSet;

// import qchem.BasisSet1.Orbital_DFT_IBS;
// import qchem.BasisSet1.Fit_IBS;

export namespace qchem::WaveFunction
{
    using obs_t=Orbital_IBS<double>;
    using bs_t=BasisSet;
}