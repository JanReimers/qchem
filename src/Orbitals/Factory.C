// File: Orbitals/Factory.C  Create orbitals.
export module qchem.Orbitals.Factory;
import qchem.Orbitals.Types;
export import qchem.ChargeDensity.Factory;   // DensityMaker (its callers hand one in)
import qchem.Orbitals.Internal.OrbitalsImp;

export namespace qchem::Orbitals
{
    //! \a maker: the run's density maker (its setup's resolved rho route): every density these orbitals build is made by it.
    TOrbitals<double>* Factory(const robs_t* bs, Spin ms, const ChargeDensity::DensityMaker& maker);
    TOrbitals<dcmplx>* Factory(const tobs_t<dcmplx>* bs, Spin ms, const ChargeDensity::DensityMaker& maker);   // plane-wave (complex) orbitals
}

