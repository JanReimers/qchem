// File: Imp/Factory.C  Create orbitals.
module qchem.Orbitals.Factory;
import qchem.Orbitals.Internal.OrbitalsImp;

namespace qchem::Orbitals
{
    TOrbitals<double>* Factory(const robs_t* bs, Spin ms, const ChargeDensity::DensityMaker& maker)
    {
        return new  TOrbitalsImp<double>(bs,ms,maker);
    }
    TOrbitals<dcmplx>* Factory(const tobs_t<dcmplx>* bs, Spin ms, const ChargeDensity::DensityMaker& maker)
    {
        return new  TOrbitalsImp<dcmplx>(bs,ms,maker);
    }
}


